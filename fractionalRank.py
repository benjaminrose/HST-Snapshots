""" fractionalRank.py -- Calculating fractional pixel rank

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-03-07
    Licesed under the MIT License
"""
from __future__ import print_function, division
import re #for regular expressions!

import numpy as np
from scipy import interpolate
from astropy.coordinates import SkyCoord, Latitude, Longitude
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
        import warnings
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

def get_galaxy_pixels(SNID, hdu, sciData=None, key=''):
    """
    Returns the numerical value of the pixels for the host galaxy of `SNID`.
    Uses the 25 mag/sqr-arcsec limit for a `key` of `hst` and `2-sigma` for a
    `key` of `sdss`.

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
    if key == 'hst':
        galaxyData = Table.read('resources/hosts_25.csv', format='ascii.commented_header', header_start=2)
    elif key == 'sdss':
        galaxyData = Table.read('resources/sdss_hosts_2.csv', format='ascii.commented_header', header_start=2)   
        # no idea which is better between these two.
        # galaxyData = Table.read('resources/sdss_hosts_2.csv', format='ascii.csv', header_start=2)
    else:
        raise ValueError('{} is an unknown source flag. Use "hst" or "sdss" to represent the galaxy source you want to use to calcualte the fpr.'.format(key))
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

def get_sn_value(position, hdu, sciData=None, key='', mask=None):
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

    key : string
        use the strings `'hst'` or `'sdss'` to designate what galaxy source you
        want to use to calcualte the fpr.

    mask : np.array (boolian)
        the reulst of `sep.mask_ellipse`. States if a pixel is inside an
        ellipse (defined by the galaxy) or not.

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
    print(position.ra.to(u.deg).value, position.dec.to(u.deg).value)
    if key == 'hst':
        SNPixels = w.all_world2pix(
        position.ra.to(u.deg).value, 
        position.dec.to(u.deg).value, 
        0
        )
    elif key == 'sdss':
        #stupid SDSS! Astropy issue #4976 should be able to make this better
        SNPixels = w.all_world2pix( 
        position.dec.to(u.deg).value,
        position.ra.to(u.deg).value, 
        0
        )
    else:
        raise ValueError('{} is an unknown source flag. Use "hst" or "sdss" to represent the galaxy source you want to use to calcualte the fpr.'.format(key))
    print('pixels: ', SNPixels)
    # can't use `SNPixels.round()` becuase `w.all_world2pix` returns a list!
    SNPixels = np.round(SNPixels)        # now a veritcal np.array

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
    return fpr

def main(key, SNID = 2635):
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

    # Notes
    This method saves three files to `'resources/SN{0}/'.format(SNID)`. They 
    are `SN*_{key}_host_pixel_values.csv`, `SN*_{key}_pixel_values.csv`, and
    `SN*_{key}_fpr.csv`. Currently it cannot make the folders, `SN*`. This
    needs to be updated.
    """
    if key == 'hst':
        # get SN position
        print('getting SN{} position'.format(SNID))
        position = get_SN_HST_coord(SNID)
        
        # get hdu and extract science data
        print('getting HDU and SciData')
        filePath = 'data/HST - combined/SN{0}_combined.fits'
        hdu, scidata = ancillary.import_fits(filePath.format(SNID), extention=1)
    elif key == 'sdss':
        # skip over sdss images I do not have
        if SNID in [12928, 15171, 19023]:
            save_location = 'resources/SN{0}/'.format(SNID)
            np.savetxt(save_location+'SN{0}_{1}_host_pixel_values.csv'.format(SNID, key), [0], delimiter=',', header="SN{0} does not have an sdss coadd image.".format(SNID))
            np.savetxt(save_location+'SN{0}_{1}_pixel_values.csv'.format(SNID, key), [0], delimiter=',', header="SN{0} does not have an sdss coadd image.".format(SNID))
            np.savetxt(save_location+'SN{0}_{1}_fpr.csv'.format(SNID, key), [0], delimiter=',', header="SN{0} does not have an sdss coadd image.".format(SNID))
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
    galaxy_pixels, mask = get_galaxy_pixels(SNID, hdu, scidata, key=key)
    print(galaxy_pixels)
    print('getting SN pixel value')
    print('SN location to pixel is wrong')
    sn_pixel_value, inside = get_sn_value(position, hdu, scidata, key, mask)
    print(sn_pixel_value)

    # get Fractional Pixel Rank
    print('getting FPR')
    #test if SN is outside galaxy, or other error
    if inside:
        fpr = get_FPR(galaxy_pixels, sn_pixel_value)
    else:
        print('SN{} is outside of galaxy'.format(SNID))
        fpr = 0
    print(fpr)

    # save data
    save_location = 'resources/SN{0}/'.format(SNID)
    #todo(make it so a new folder can be made if need be.)
    np.savetxt(save_location+'SN{0}_{1}_host_pixel_values.csv'.format(SNID, key), galaxy_pixels, delimiter=',', header="The pixel values of SN{0} host galaxy from {1}. The galaxy's edge is defined hard coded currently".format(SNID, key))
    np.savetxt(save_location+'SN{0}_{1}_pixel_values.csv'.format(SNID, key), [sn_pixel_value], delimiter=',', header='The pixel value of SN'+str(SNID)+' as seen by '+key)
    np.savetxt(save_location+'SN{0}_{1}_fpr.csv'.format(SNID, key), [fpr], delimiter=',', header='The fractional pixel rank calculated for SN'+str(SNID)+' calcuated with '+key)

if __name__ == "__main__":
    # main('hst', SNID = 12781)
    
    # SN = np.array([2102, 2635], dtype=np.int)
    # flag = ['sdss']*len(SN)
    # map(main, flag, SN)

    SN = ancillary.get_sn_names()
    SN = np.array(SN, dtype=np.int)
    flag = ['hst']*len(SN)
    map(main, flag, SN)

    ###########################################
    ## Testing if pixel value is correct for SN
    ###########################################
    # SNID = 2635
    # print('getting SN position')
    # position = get_SN_HST_coord(SNID)

    # # get hdu and extract science data
    # print('getting HDU and SciData')
    # filePath = 'data/HST - combined/SN{0}_combined.fits'
    # hdu, scidata = ancillary.import_fits(filePath.format(SNID), extention=1)

    # print('getting SN pixel value')
    # print('SN location to pixel is wrong')
    # sn_pixel_value = get_sn_value(position, hdu, scidata)
    # print(sn_pixel_value)



    # SNID = 2635
    # position = get_SN_HST_coord(SNID)

    # # galaxy, SN = get_pixels(SNID, position) #a sup-part of `get_FPR`, for testing
    # # print galaxy, SN

    # fpr = get_FPR(SNID, position)#, position)
    # print fpr