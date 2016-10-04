""" fractionalRank.py -- Calculating fractional pixel rank and fractional flux

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-07
    Licesed under the MIT License
"""
from __future__ import print_function, division
import re #for regular expressions!
import warnings

import numpy as np
from scipy import interpolate
from astropy.coordinates import SkyCoord, Latitude, Longitude
#block_reduce needs scikit-image
from astropy.nddata.utils import block_reduce
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS
from sep import mask_ellipse

import ancillary

def get_SN_SDSS_coord(SNID):
    """
    Gets the cooridnates, in SDSS system, of a given SN defined by its SDSS SN 
    identification number. The data is taken from photometric data stored in 
    `data/SDSS - photometry/` from 
    `http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html`.
    #todo(figure out the paper that this is from)

    # Parameters 
    SNID: int
        The identification number of the SDSS Supernova, used in file names.
    
    # Returns
    SNPosition : astopy.SkyCoord
        The position of the SN found in SDSS Photomentric data. 
    """
    data_location='data/SDSS - photometry/'

    # read from SDSS files SN's RA & Dec
    SNID_string = str(SNID).zfill(6)        #pad with zeros to match SMP
    with open(data_location+'SMP_{0}.dat'.format(SNID_string), 'r') as f:
        first_line = f.readline()
    
    # split on white space, convert numpy array so np.where works
    split = np.array( re.split(r'\s+', first_line) )
    RAValueLocation = np.where(np.array(split) == 'RA:')[0][0]+1
    DecValueLocation = np.where(np.array(split) == 'DEC:')[0][0]+1

    # create SkyCoord
    RA = float(split[RAValueLocation])*u.deg
    Dec = float(split[DecValueLocation])*u.deg
    
    SNPosition = SkyCoord(ra = RA, dec = Dec)
    #todo(Determing if this is the correct frame: FK5, icrs, etc.)

    return SNPosition

def get_SN_HST_coord(SNID):
    """
    Gets the cooridnates, in HST system, of a given SN defined by its SDSS SN 
    identification number. The data has a shift applied to change from the 
    SDSS to HST WCS.

    # Parameters
    SNID : int
        The identification number of the SDSS Supernova, used in file names.

    # Returns
    SNPosition : astopy.SkyCoord
        The position of the SN, with SDSS-to-HST shift applied. Shifts are 
        defined in `resources/shift.csv` and reasoning can be found in 
        `README.md`
    """
    # get position from SDSS information
    SDSS_SNPosition = get_SN_SDSS_coord(SNID)
    print('SDSS position: ', SDSS_SNPosition.to_string(style=u'hmsdms'))

    #skip 13038 because I don't have a shift
    if SNID == 13038:
        warnings.warn('SN13038 does not have a shift')
        SNPosition = SDSS_SNPosition
    else:
        # import shift data
        #todo(change this location to be dataset.csv)
        shift = Table.read('resources/shift.csv', format='ascii.commented_header')
        #todo(add units to table, from shift.meta)
        shift['delta RA'].unit, shift['delta Dec'].unit = u.arcsec, u.arcsec

        # apply shift with delta = HST - SDSS or HST = delta + SDSS
        # `[0]` are added to drop unneeded array wrapper. 
        # It also makes an internal standard between hst-coords and sdss-cords
        deltaRA = Longitude(
            shift[shift['SDSS SN Name'] == SNID]['delta RA'].quantity[0]
            )
        deltaDec = Latitude(
            shift[shift['SDSS SN Name'] == SNID]['delta Dec'].quantity[0]
            )
        SNPosition = SkyCoord(ra = deltaRA + SDSS_SNPosition.ra,
                              dec = deltaDec + SDSS_SNPosition.dec)
        print('deltas: ', deltaRA, deltaDec)
        print('SNPosition: ', SNPosition.to_string(style=u'hmsdms'))
    return SNPosition

def get_galaxy_pixels(SNID, hdu, sciData=None, key='', cutoff=0, block_size=None):
    """
    Returns the numerical value of the pixels for the host galaxy of `SNID`.
    Uses the 26 mag/sqr-arcsec limit for a `key` of `hst` and `2-sigma` for a
    `key` of `sdss`.
    #todo(change to a variable sb cut off)

    # Parameters
    SNID : int
        The numerical SDSS ID for the SN.

    hdu : astropy.io.fits.hdu.hdulist.HDUList
        The fits object from the data extension. This needs to cantain a 
        WCS. Assumes its an HST style object with WCS in extention `1`.

    sciData : np.array, optional
        A 2D array of the sciecne data from the CCD Chip. If science data is 
        not given then it will be extracted. Provide pre-extracted data if any 
        manipulation is needed prior to analaysis, ie byteswap.

    key : string
        use the strings `'hst'` or `'sdss'` to designate what galaxy source you
        want to use to calcualte the fpr.

    block_size : int
        The integer block size. Used to downsample a data array by applying a 
        function to local blocks.

    # Returns
    galaxy : np.array
        A 1D array of all the value of the pixels of the galaxy. Unsorted, 
        or not?
    """
    # get science data if needed
    if sciData is None:
        sciData = hdu.data

    # init mask that defines pixels belonging to host
    mask = np.zeros(sciData.shape, dtype=np.bool)

    # get galaxy definition & format its astorpy.table object
    # galaxyData = Table.read('resources/2635-galaxy.csv', format='ascii.commented_header', header_start=1)    #hardcoded for testing
    try:
        if block_size is None:
            galaxyData = Table.read('resources/{0}_hosts_{1}.csv'.format(key, cutoff), format='ascii.commented_header', header_start=2, delimiter=',')
        else:
            galaxyData = Table.read('resources/{0}_hosts_{1}_{2}.csv'.format(key, cutoff, block_size), format='ascii.commented_header', header_start=2, delimiter=',')
    #note IOError is for python2 and FileNotFoundError is for python3
    except IOError as e:
        raise IOError(e+' change key ({0}), cutoff ({1}) or block_size ({2}) to find the correct file.'.format(key, cutoff, block_size))
    # table headers are: SN ID, npix, x, y, a, b, theta
    # add apropriate units to the table.
    galaxyData['theta'].unit = u.radian  
    # make the rows follow the SNID
    # SN ID might be 'SN number' or 'SN name' or 'SN'. I have issues with conistency
    galaxyData.add_index('SN')
    #todo(I might need to not read this in each time, but rather take it as a parmeter?)
    
    hostParams = galaxyData.loc[SNID]
    #todo(2635 needs to be a parameter)
    hostParams = Table(hostParams)     #convert to a Table because Rows suck

    # set scalling factor for host parameters. Silly `SE` and `sep`!
    n = 3 
    print(hostParams['x'].quantity[0], hostParams['y'].quantity[0],
                 hostParams['a'].quantity[0], hostParams['b'].quantity[0], 
                 hostParams['theta'].to(u.radian).value[0], n)
    mask_ellipse(mask, hostParams['x'].quantity[0], hostParams['y'].quantity[0]
                , hostParams['a'].quantity[0], hostParams['b'].quantity[0] 
                , hostParams['theta'].to(u.radian).value[0], n)

    # create a nd.array of the pixel values inside the galaxy.
    sciDataFlattened = (mask*sciData).flatten()
    ranked = np.sort(sciDataFlattened[sciDataFlattened != 0])

    return ranked, mask

def get_sn_value(position, hdu, sciData=None, mask=None, key='', cutoff=0):
    """
    Returns the numerical value of the pixels for the SN

    # Parameters
    position : astropy.SkyCoord
        The position of the SN. It will be used with the WCS found in `hdu`.

    hdu : astropy.io.fits.hdu.hdulist.HDUList
        The fits object from the data extension. This needs to cantain a 
        WCS. Assumes its an HST style object with WCS in extention `1`.

    sciData : np.array, optional
        A 2D array of the sciecne data from the CCD Chip. If science data is 
        not given then it will be extracted. Provide pre-extracted data if any 
        manipulation is needed prior to analaysis, ie byteswap.

    mask : np.array (boolian)
        the reulst of `sep.mask_ellipse`. States if a pixel is inside an
        ellipse (defined by the galaxy) or not.

    key : string
        use the strings `'hst'` or `'sdss'` to designate what galaxy source you
        want to use to calcualte the fpr.

    cutoff : float
        This matches the cutoff value for the galaxy definition. Likely 26 for 
        `'hst'`, and either 1, 2, or 3 for `'sdss'`. 

    # Returns
    sn : float
        The pixel value of the.
    
    inside : boolian
        A flag stating if the SN is inside or outside the galaxy.
    """
    # get science data if needed
    if sciData is None:
        sciData = hdu.data

    # Clean up position
    #The results from HST

    # get world cordinate system
    # can't just do `WCS('filename')` because of HST has multihead fits files
    if key == 'hst':
        w = WCS(hdu[1].header)
    elif key == 'sdss':
        w = WCS(hdu[0].header)
    else:
        raise ValueError('{} is an unknown source flag. Use "hst" or "sdss" to represent the galaxy source you want to use to calcualte the fpr.'.format(key))

    # Get SN's pixle position
    # astopy.wcs.WCS only degrees as floats
    # origin is `0` because we will be working with `sciData`
    #todo(this is doing this wrong. It does not agree with what DS9 says the pixel value should be for that WCS location.)
    print('position: ', position)
    print(position.to_string('hmsdms'))
    #Testing on SN12781 showes that the WCS are 1-based for both HST and SDSS.
    if key == 'hst':
        SNPixels = w.all_world2pix(
        position.ra.to(u.deg).value, 
        position.dec.to(u.deg).value, 
        1
        )
    elif key == 'sdss':
        #stupid SDSS! Astropy issue #4976 should be able to make this better
        SNPixels = w.all_world2pix( 
        position.dec.to(u.deg).value,
        position.ra.to(u.deg).value, 
        1
        )
    else:
        raise ValueError('{} is an unknown source flag. Use "hst" or "sdss" to represent the galaxy source you want to use to calcualte the fpr.'.format(key))
    #todo(add an assert to test if I am pulling from the WCS correctly.)
    # can't use `SNPixels.round()` becuase `w.all_world2pix` returns a list!
    SNPixels = np.round(np.array(SNPixels))        # now a veritcal np.array
    print('pixels: ', SNPixels)

    # Test if SN is outide the galaxy
    #todo(what if mask was not pasted in?)
    inside = mask[SNPixels[1], SNPixels[0]]

    # Get value of SN's pixel
    #take median of 3x3 box. 
    sn = np.array([])
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            # numpy.arrays are accessed row then collumn and `all_pix2wold()`
            # returns (x, y) so we need to call sciData[row #,column #]
            sn = np.append(sn, sciData[SNPixels[1]+i, SNPixels[0]+j])
    sn = np.median(sn)
    #todo(do I want median or mean?)

    return sn, inside


def get_FPR(galaxy, SN):#, positions, sigma=2, box_size=3):
    """
    Fractional Pixel Rank (FPR) is CDF value of a particular number, in this 
    case of the pixel that hosted the SN. It is ok if `SN` is or is not found 
    in `galaxy`.

    # Parameters
    galaxy : np.array
        The float value of the galaxy pixels. 

    SN : float
        The value of the pixel (actual, median of a 3x3 grid or etc) where the 
        SN went off.

    # Returns
    fpr : float
        The fractional pixel rank. 
    """
    #todo(if SN pixel value is less then galaxy edge cut off (look this up), then break and return 0 or -1? This minght need to be done in `get_pixels()` or `get_pixels()` might not be a useful function.)

    # calculate the FPR with range [0,1] both inclusive.
    galaxy.sort()
    
    # You can't simple do the interpelation for the whole CDF. CDFs are too 
    # jumpy/verticle for `scipy.interpolate` to work well. It works fine 
    # between two points but not for the whole range. See notes on 2016-03-14.

    # since galaxy is a 1d array, there is a non-needed tuple wrapper around the result of `np.where()`
    rank = np.where(galaxy == SN)[0]

    #if `SN` is found in `galaxy` calculate its FPR (or the mean of its FPRs) driectly
    if len(rank) > 0:
        # subtact 1 so that the denominator is equal to the largest rank can be.
        # Notes are available from 2015-03-09
        fpr = 1.0*rank/(len(galaxy)-1.0)
        if len(fpr) > 1:
            fpr = fpr.mean()
    # if `SN` is NOT found in `galaxy` calculate its FPR by interpolation
    else:
        # find nearist neighbors in galaxy and coresponding CDF values
        # find last entry position of where galaxy is less then N
        # `np.where()` returns a tuple where we only care about arg-0.
        #breaks if SN value is smaller then galaxy. #todo(fix)
        rank_min = np.where(galaxy[galaxy < SN])[0][-1]
        print(rank_min, len(galaxy))
        # find first entry of where galaxy is greater then SN (or add 1 to `rank_min`)
        rank_max = rank_min + 1
        fpr_min = 1.0*rank_min/(len(galaxy)-1.0)
        fpr_max = 1.0*rank_max/(len(galaxy)-1.0)
        if fpr_min == 1:
            fpr = 1.0
        else:
            f = interpolate.interp1d([galaxy[rank_min], galaxy[rank_max]]
                                    ,[fpr_min, fpr_max])
            fpr = f(SN)

    fractionalFlux = galaxy[galaxy<SN].sum()/galaxy.sum()

    return fpr, fractionalFlux

def saveNaN(SNID, key, cutoff, header, block_size):
    save_location = 'resources/SN{0}/'.format(SNID)
    if block_size is None:
        np.savetxt(save_location+'SN{0}_{1}_{2}_host_pixel_values.csv'.format(SNID, key, cutoff), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_pixel_values.csv'.format(SNID, key, cutoff), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_fpr.csv'.format(SNID, key, cutoff), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_fractionalFlux.csv'.format(SNID, key, cutoff), [np.nan], delimiter=',', header=header)
    else:
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_host_pixel_values.csv'.format(SNID, key, cutoff, block_size), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_pixel_values.csv'.format(SNID, key, cutoff, block_size), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_fpr.csv'.format(SNID, key, cutoff, block_size), [np.nan], delimiter=',', header=header)
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_fractionalFlux.csv'.format(SNID, key, cutoff, block_size), [np.nan], delimiter=',', header=header)

def main(key, SNID, cutoff, block_size=None):
    """
    This is the default method for calculucating fractional pixel rank for a
    single SN and it's SDSS host. This is set up so `map()` can be used to
    iterate through all the SN while calling `main_sdss()`.

    # Parameters
    key : string
        use the strings `'hst'` or `'sdss'` to designate what galaxy source you
        want to use to calcualte the fpr.

    SNID : int
        The numerical SDSS ID for the SN.

    cutoff : float
        This matches the cutoff value for the galaxy definition. Likely 26 for 
        `'hst'`, and either 1, 2, or 3 for `'sdss'`.

    block_size : int
        The integer block size. Used to downsample a data array by applying a 
        function to local blocks.

    # Notes
    This method saves three files to `'resources/SN{0}/'.format(SNID)`. They 
    are `SN*_{key}_host_pixel_values.csv`, `SN*_{key}_pixel_values.csv`, and
    `SN*_{key}_fpr.csv`. Currently it cannot make the folders, `SN*`. This
    needs to be updated.
    """
    #skip over two objects with obstructions
    if SNID in [6491, 15345]:
        header = 'SN{0} has an object obstructing its host galaxy'.format(SNID)
        # display warning that there is an issue. 
        warnings.warn(header)
        saveNaN(SNID, key, cutoff, header, block_size)
        return None
    if key == 'hst':
        # get SN position
        print('getting SN{} position'.format(SNID))
        position = get_SN_HST_coord(SNID)
        
        # get hdu and extract science data
        print('getting HDU and SciData')
        if block_size:
            filePath = 'data/HST - reduced/SN{0}_reduced{1}.fits'.format(SNID,
                                                                  block_size)
        else:
            filePath = 'data/HST - combined/SN{0}_combined.fits'.format(SNID)
        hdu, scidata = ancillary.import_fits(filePath, extention=1)
    elif key == 'sdss':
        # skip over sdss images I do not have
        if SNID in [12928, 15171, 19023]:
            header = "SN{0} does not have an sdss coadd image.".format(SNID)
            saveNaN(SNID, key, cutoff, header, block_size)
            return None

        # get SN position
        print('getting SN{} position'.format(SNID))
        position = get_SN_SDSS_coord(SNID)
        
        # get hdu and extract science data
        print('getting HDU and SciData')
        filePath = 'data/SDSS - coadded/SN{0}-g.fits'
        if SNID in [12928, 15171, 19023]:
            # These dont exists yet. Also in defGalaxy.py:l246
            return None
        else:
            hdu, scidata = ancillary.import_fits(filePath.format(SNID))
    else:
        raise ValueError('{} is an unknown source flag. Use "hst" or "sdss" to represent the galaxy source you want to use to calcualte the fpr.'.format(key))

    # get pixel values of galaxy and SN
    print('getting galaxy pixels values')
    galaxy_pixels, mask = get_galaxy_pixels(SNID, hdu, scidata, key, cutoff, block_size)
    print(galaxy_pixels)
    print('getting SN pixel value')
    print('SN location to pixel is wrong')
    sn_pixel_value, inside = get_sn_value(position, hdu, scidata, mask, key, cutoff)
    print(sn_pixel_value)
    
    # get Fractional Pixel Rank
    print('getting FPR')
    #test if SN is outside galaxy, or other error
    if inside:
        fpr, fractionalFlux = get_FPR(galaxy_pixels, sn_pixel_value)
    else:
        print('SN{} is outside of galaxy'.format(SNID))
        fpr, fractionalFlux = 0, 0
    print(fpr, fractionalFlux)

    # save data
    save_location = 'resources/SN{0}/'.format(SNID)
    #todo(make it so a new folder can be made if need be.)
    if block_size is None:
        np.savetxt(save_location+'SN{0}_{1}_{2}_host_pixel_values.csv'.format(SNID, key, cutoff), galaxy_pixels, delimiter=',', header="The pixel values of SN{0} host galaxy from {1}. The galaxy's edge is defined hard coded currently".format(SNID, key))
        np.savetxt(save_location+'SN{0}_{1}_{2}_pixel_values.csv'.format(SNID, key, cutoff), [sn_pixel_value], delimiter=',', header='The pixel value of SN'+str(SNID)+' as seen by '+key)
        np.savetxt(save_location+'SN{0}_{1}_{2}_fpr.csv'.format(SNID, key, cutoff), [fpr], delimiter=',', header='The fractional pixel rank calculated for SN'+str(SNID)+' calculated with '+key)
        np.savetxt(save_location+'SN{0}_{1}_{2}_fractionalFlux.csv'.format(SNID, key, cutoff), [fractionalFlux], delimiter=',', header='The fractional flux calculated for SN'+str(SNID)+' calculated with '+key)
    else:
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_host_pixel_values.csv'.format(SNID, key, cutoff, block_size), galaxy_pixels, delimiter=',', header="The pixel values of SN{0} host galaxy from {1}. The galaxy's edge is defined hard coded currently".format(SNID, key))
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_pixel_values.csv'.format(SNID, key, cutoff, block_size), [sn_pixel_value], delimiter=',', header='The pixel value of SN'+str(SNID)+' as seen by '+key)
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_fpr.csv'.format(SNID, key, cutoff, block_size), [fpr], delimiter=',', header='The fractional pixel rank calculated for SN'+str(SNID)+' calculated with '+key)
        np.savetxt(save_location+'SN{0}_{1}_{2}_{3}_fractionalFlux.csv'.format(SNID, key, cutoff, block_size), [fractionalFlux], delimiter=',', header='The fractional flux calculated for SN'+str(SNID)+' calculated with '+key)

if __name__ == "__main__":
    #Single SN testing
    #hst block reduced
    # main('hst', 6057, 26, 8)
    #sdss 
    # main('sdss', 12781, 2)
    

    #Set up parameters 
    SN = ancillary.get_sn_names()
    SN = np.array(SN, dtype=np.int)
    flag = ['sdss']*len(SN)
    reduced = [8]*len(SN)

    #run once
    # map(main, flag, SN)
    # map(main, flag, SN, reduced)

    #run all settings
    for i in [('hst', 26), ('hst', 26, 8), ('sdss', 1), ('sdss', 2), ('sdss', 3)]:
        flag = [i[0]]*len(SN)
        cutoff = [i[1]]*len(SN)
        if len(i) == 3:
            #this is redudent to the part in set up parameters
            reduced = [i[2]]*len(SN)
            map(main, flag, SN, cutoff, reduced)
        else:
            map(main, flag, SN, cutoff)