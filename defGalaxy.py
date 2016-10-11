""" defGalaxy.py -- set of functions to define galaxy shape size to 
    `resources/galaxies_{0}.csv` where `{0}` is the edge definition.

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-03
    Licesed under the MIT License
"""
from __future__ import print_function, division
from sys import exit
from datetime import datetime
from os import path, makedirs
import warnings

import numpy as np
import sep

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs #for WCS, and util
from astropy.io import fits     #to read fits files

import ancillary 
import fractionalRank

def smooth_Image(sciData, stDev):
    """
    this takes an image and smooths it via `astropy.convolution`.

    # Parameters
    sciData : np.ndarray
        sci data from an hdu or other image, as a 2D numpy array.

    stDev : int
        The standard deviation used as a parameter in 
        `astropy.convolution.Gaussian2DKernel`.


    # Returens
    smoothedData : np.ndarry
        the reuslt of `astropy.convolution.convolve`

    """
    gauss = Gaussian2DKernel(stDev)
    smoothedData = convolve(sciData, gauss) 
    return smoothedData

def run_sep(sciData, threshold, minarea=50, deblendCont=0.1):
    """
    this takes an image and smooths it

    # Parameters
    sciData : np.ndarray
        sci data from an hdu or other image, as a 2D numpy array.

    threshold : float
        Pixel value used for galaxy edge definion. 

    minarea : int
        The minimum area of the source. Passed directly to `sep`.

    deblendCont : float
        One of the deblending parrameters. Passed directly to `sep`'s `deblend_cont`.

    # Returens
    sources : np.ndarray
        the sources, in a structured array, found from sep.
    """
    bkg = sep.Background(sciData)
    bkg.subfrom(sciData)

    sources = sep.extract(sciData, threshold, minarea=minarea, deblend_cont=deblendCont)
    return sources

def find_host(sources, initialGuess=(2090/2.0, 2108/2.0), searchRadius=200):
    """
    A search to find the host galaxy from a list of possible sources. It retuns the largest source (from `source`) found within the pixel `searchRadius` of the `initialGuess`.

    # Paramenters
    sources : np.ndarray
        an structured array of sources, in the format returned by 
        [`sep`](https://github.com/kbarbary/sep). Need values 
        `['npix', 'x', 'y', 'a', 'b', 'theta']`.

    initialGuess : tuple
        The likely pixel location of the host. Needs to be of length 2. This 
        defaults to the center of an HST image but the SN's location can be 
        used instead.

    searchRadius : int
        The number of pixels from the `initialGUess` to search.

    # Returns
    host : np.void
        A list of all the parameters from sep of the determined host galaxy.
        It still a structured array that contains fields named:
        `['npix', 'x', 'y', 'a', 'b', 'theta']`
    """
    #todo(this fails. It is not selecting the right thing at all.)
    # Where, in array `soucres`, are the objects close to initialGuess?

    #make a holding varriable
    centerIDs = []
    #loop through `sources` and save if its center is inside search area
    for i, x in enumerate(sources['x']):
        if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
            centerIDs.append(i)
    
    #check to see if something is found, or else `np.argmax` fails badly
    if len(centerIDs) == 0:
        warnings.warn("This SN can't be found in initial search")
        centerIDs = []    #just play it safe

        #increase search area till you find somthing. 
        while len(centerIDs) == 0: 
            #increase from 200 -> 500 for hst or propotinally
            searchRadius *= 2.5
            print('searchRadius: ', searchRadius)
            print('initialGuess: ', initialGuess)
            for i, x in enumerate(sources['x']):
                if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
                    centerIDs.append(i)
            if searchRadius > 1000:
                from sys import exit; exit('You fail and got a search radius of {}.'.format(searchRadius))
    
    # Select largest of the center objects, but save the ID
    #this selects the centerID associated with the max (in size) of the central cources
    idx = centerIDs[np.argmax(sources['npix'][centerIDs])]
    #todo(add error for too many found)
    #todo(implement a way to fail greacefully if argmax is empty)

    host = sources[['npix', 'x', 'y', 'a', 'b', 'theta']][idx]
    return host

def find_sdss_host(sources, SNID, hdu):
    """
    This finds the likly host from the souces in the sdss field

    # Parameters
    sources : np.ndarray
        an structured array of sources, in the format returned by 
        [`sep`](https://github.com/kbarbary/sep). Need values 
        `['npix', 'x', 'y', 'a', 'b', 'theta']`.

    SNID : int
        The identification number of the SDSS Supernova, used in file names.

    hdu : astropy.io.fits.hdu.hdulist.HDUList
        The fits object from the data extension. This needs to cantain a 
        WCS. Assumes its an SDSS style object with WCS in extention `0`.

    # Returns
    host : np.void
        A list of all the parameters from sep of the determined host galaxy.
        It still a structured array that contains fields named:
        `['npix', 'x', 'y', 'a', 'b', 'theta']`
    """
    SNCoords = fractionalRank.get_SN_SDSS_coord(SNID)

    # get world cordinate system
    w = wcs.WCS(hdu[0].header)

    # Get SN's pixle position
    # astopy.wcs.WCS only degrees as floats
    # origin is `0` because we will be working with `sciData`
    # SDSS has the stupid WCS backwards
    print(SNCoords.ra.to(u.deg).value, SNCoords.dec.to(u.deg).value)
    SNPixels = w.all_world2pix( 
        SNCoords.dec.to(u.deg).value, 
        SNCoords.ra.to(u.deg).value,
        0
        )

    # can't use `SNPixels.round()` becuase `w.all_world2pix` returns a list!
    SNPixels = np.round(SNPixels)        # now a veritcal np.array

    radius = 25

    host = find_host(sources, SNPixels, radius)
    return host

def saveGalaxies(sources, host, SNNumber, sb, telescope='hst', block_size=1):
    """
    This saves the results of this script using a sytematic convention
    
    # Parameters
    sources : np.ndarray
        The output of `sep`, a structured array of all the objects.

    host : np.ndarray
        The pertinate values of the suspected host. The array should be the 
        result of calling `[['npix', 'x', 'y', 'a', 'b', 'theta']]` on the 
        stuctured array returned by `sep` for this object.

    SNNumber : int, str
        This is the number of the SDSS SN that this data coresponds too.

    sb : float
        The surface brightness threshold that the soucers were found at, in
        mag/sqr-arcsec

    telescope : str
        String that represents the telescope that these hosts were observed
        on. This allows source definintion to be from multiple observations
        without things getting over written.

    block_size : int
        The integer block size. Used to downsample a data array by applying a 
        function to local blocks. Defaults to `1` and is therefore not apart of the output file names. 
    """

    #Save ALL objects
    if block_size == 1:
        allHeader = 'data from sep on SN{} as observed by '.format(SNNumber) + telescope + ' with a SB cutoff of {0} mag/sqr-arcsec on {1}'.format(sb, datetime.now().date()) + '\n' + str(sources.dtype.names)[1:-1].replace("'","")
        allLocationFIle = telescope+'_sources_{}.csv'.format(sb)
    else:
        allHeader = 'data from sep on SN{} as observed by '.format(SNNumber) + telescope + ' with a SB cutoff of {0} mag/sqr-arcsec and a block size of {1} on {2}'.format(sb, block_size, datetime.now().date()) + '\n' + str(sources.dtype.names)[1:-1].replace("'","")
        allLocationFIle = telescope+'_sources_{0}_{1}.csv'.format(sb, block_size)
    allLocationFolder = 'resources/SN{}'.format(SNNumber)
    allLocation = allLocationFolder+'/'+allLocationFIle
    if not path.exists(allLocationFolder): makedirs(allLocationFolder)

    np.savetxt(allLocation, sources, delimiter=',', header=allHeader)

    
    #Save host
    #todo(don't assume this file existes, but rather check.)
    if block_size == 1:
        hostLocation = 'resources/' + telescope + '_hosts_{}.csv'.format(sb)
    else:
        hostLocation = 'resources/' + telescope + '_hosts_{0}_{1}.csv'.format(sb, block_size)
    #combine `SNNumber` to the front of `host` data, but structured arrays are stupid. So is np.savetxt()
    dataToSave = str(SNNumber)
    for i in host:
        dataToSave += ', ' + str(i)
    dataToSave += '\n'

    with open(hostLocation, 'a') as f:
        f.write(dataToSave)

def main_hst(SNNumber = 2635, block_size=1, surfaceBrightness = 26, minArea=50, deblendCont=0.1):
    """
    This is the default method of defining an hst galaxy

    # Parameters
    SNNumber : int
        The sdss identification number of the object you want to find the host 
        for.

    block_size : int
        The integer block size. Used to downsample a data array by applying a 
        function to local blocks.

    surfaceBrightness : int, float
        the disired cut off in mag/sqr-arcsec

    minarea : int
        The minimum area of the source. Passed directly to `sep`. Will be devided by block_size**2.

    deblendCont : float
        One of the deblending parrameters. Passed directly to `sep`'s `deblend_cont`.
    """
    print("running SN{}".format(SNNumber))
    # imageFile = 'data/HST - combined/SN{}_combined.fits'.format(SNNumber)
    if block_size == 1:
        imageFile = 'data/HST - combined/SN{}_combined_flux.fits'.format(SNNumber)
    else:
        imageFile = 'data/HST - reduced/SN{}_reduced{}.fits'.format(SNNumber, block_size)

    #import image data
    data = ancillary.import_fits(imageFile, extention=1)[1]

    #smooth image
    if block_size == 1:
        stDev = 3
        data = smooth_Image(data, stDev)

    # find threshold
    #the magniutde per pixel for the same SB, wikipedia for equation
    magPerPixel = surfaceBrightness - 2.5*np.log10((0.05*block_size)**2)
    #convert to AB-mag to flux
    #more at http://www.stsci.edu/hst/acs/analysis/zeropoints
    fNu = 3.63e-20*(u.erg / u.cm**2 / u.s / u.Hz)     #this is when AB-mag = 0
    #convert from f_nu to f_lambda
    f475 = fNu.to(u.erg / u.cm**2 / u.s / u.nm, equivalencies=u.spectral_density(475*u.nm))    #todo(maybe use the pivot wavelenght instead)
    f625 = fNu.to(u.erg / u.cm**2 / u.s / u.nm, equivalencies=u.spectral_density(625*u.nm))
    #this is the zero-mag flux combined like the data
    fluxAB = (f475 + f625)/2.0
    fluxThresh = fluxAB * 10**(magPerPixel/-2.5)

    #run sep
    #fix `minArea`
    minArea = minArea/block_size**2
    #do not let the value go to less then 5 (default) because of block_size
    if minArea < 5:
        minArea = 5
    sources = run_sep(data, fluxThresh.value, minArea, deblendCont)
    # print('sources: ', sources[['npix', 'x', 'y', 'a', 'b', 'theta']])

    #get "best" object from sep - 
    #or select from hardcoded objects that do not work in the test.
    hardcodeFindHost = [8297, 13038, 13354, 13411, 14113, 14284, 18415, 19282]
    #todo(why am i repeateing the obove line and then swiching through everything below?)
    if SNNumber in hardcodeFindHost:
        # if SN is hardcoded do the search, but with a very small radius and a guess that is what is determined between the R_25 and R_26 runs that were visually inspected.
        warnings.warn("SN{}'s host is searched via a hard-coded method".format(SNNumber))
        radius = 20/block_size     #pixels
        if SNNumber == 8297:
            #todo(fix this crazyness. Is this needed?)
            # Does not work for R_25, and R_26 is just fine.
            # SNPixels = (1111, 1025)
            # host = find_host(sources, SNPixels, radius)
            host = find_host(sources, (1045.0/block_size, 1054.0/block_size), 200.0/block_size)
        elif SNNumber == 13354:
            SNPixels = (978/block_size, 1043/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 13411:
            SNPixels = (1010/block_size, 1067/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 14113:
            SNPixels = (1041/block_size, 1049/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 14284:
            SNPixels = (1063/block_size, 1056/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 18415:
            SNPixels = (1095/block_size, 1066/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 19282:
            SNPixels = (1040/block_size, 1059/block_size)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 13038:
            SNPixels = (1059/block_size, 1033/block_size)
            host = find_host(sources, SNPixels, radius)
        else:
            raise NotImplementedError('Somehow SN{} is designated for HST hardcoding but is not implemented'.format(SNNumber))
    else:
        host = find_host(sources, (1045.0/block_size, 1054.0/block_size), 200.0/block_size)
    print('host: ', host)

    #save resutls to file
    saveGalaxies(sources, host, SNNumber, surfaceBrightness, block_size=block_size)

def main_sdss(SNNumber=2635, sigma=2, fltr='g', minarea=5, deblendCont=0.005):
    """
    This is the default method of defining an sdss galaxy

    # Parameters
    SNNumber : int
        The sdss identification number of the object you want to find the host 
        for.

    fltr : string, char
        The single character of an SDSS filter: u, g, r, i, or z (sill python 
        with `filter` being a built in function!)

    minarea : int
        The minimum area of the source. Passed directly to `sep`. The default 
        is the same as `sep`.

    deblendCont : float
        One of the deblending parrameters. Passed directly to `sep`'s `deblend_cont`. The default is the same as `sep`.

    sigma : float
        The value mutiplied to the global background RMS to calcuate the 
        threshold for galaxy edge.
    """
    #break out if I don't have these sdss files
    if SNNumber in [12928, 15171, 19023]:
        # These dont exists yet. Also in fractionalRank.py:l322
        return None
    print("running SN{}".format(SNNumber))

    imageFile = 'data/SDSS - coadded/SN{0}-{1}.fits'.format(SNNumber, fltr)

    #import image data
    hdu, data = ancillary.import_fits(imageFile, extention=0)

    # find threshold
    bkg = sep.Background(data)
    thresh = sigma*bkg.globalrms
    #run sep
    sources = run_sep(data, thresh, minarea, deblendCont)
    #todo(what is happening here?)
    np.savetxt('temp-sn2635-sdss-sources.csv', sources, delimiter=',', header='temp sn2635 sdss sources')
    print('sources: ', sources[['npix', 'x', 'y', 'a', 'b', 'theta']])

    #get "best" object from sep - 
    #or select from hardcoded objects that do not work in the test.
    hardcoded = [13354, 14113, 14284, 14437, 19282, 20350]
    if SNNumber in hardcoded:
        radius = 10    # the normal sdss search is 25 pixel search
        # hardcoed in the initial guess and search radius, the sdss method did not work on these.
        if SNNumber == 13354:
            #We can run just to find_host() because find_sdss_host() just searches for SNPixels and uses a larger radius.
            #todo(maybe all of these can be cut if the radius in find_sdss_host was smaller?)
            # SNPixels = (1539, 1283)    #SN's position is not helpful
            SNPixels = (1544, 1284)      #need to shift towards galaxy
            radius = 8
        elif SNNumber == 14113:
            SNPixels = (259, 1115)
        elif SNNumber == 14284:
            SNPixels = (328, 612)
        elif SNNumber == 14437:
            SNPixels = (637, 454)
        elif SNNumber == 18415:
            SNPixels = (159, 526)
        elif SNNumber == 19282:
            SNPixels = (1193, 849)
        elif SNNumber == 20350:
            SNPixels = (917, 796)
            radius = 15    #if I let it auto exapnd, it picks the neighbor.
        else:
            raise NotImplementedError('Somehow SN{} is designated for SDSS hardcoding but is not implemented'.format(SNNumber))
        host = find_host(sources, SNPixels, radius)
    else:
        host = find_sdss_host(sources, SNNumber, hdu)
    print('host: ', host)

    #save resutls to file
    saveGalaxies(sources, host, SNNumber, sigma, 'sdss')

def runSEPIndiviually():
    """     
    Run a few SN (through `main_hst`) with specific settings. These need
    to be saved, but this concept does not work well. You need to make sure 
    that these SN are saved properly. Currently they are just appened.

    This should be a data oject in the main-scope and `main_hst` and `main_sdss`
    should search this object to see if they need to adjust any settings. But
    since nothing is actually helpful yet (we only have new settings for SN13038
    and that changes the SB-cutoff) I will not implement that yet.
    """
    # LSB galaxy, Nothing is found with these settings!
    # SN19023 = {
    #     'number' : 19023,
    #     'sb' : 26.5,
    #     'minArea' : 5, 
    #     'deblendCont' : 0.1
    # }

    # LSB galaxy, Nothing is found with these settings!
    # SN15850 = {
    #     'number' : 15850,
    #     'sb' : 26.5,
    #     'minArea' : 5, 
    #     'deblendCont' : 0.1
    # }

    # LSB galaxy, Found with these settings!
    SN13038 = {
        'number' : 13038,
        'sb' : 26.5,
        'minArea' : 50, 
        'deblendCont' : 0.1
    }

    main_hst(SN13038['number'], SN13038['sb'], minArea=SN13038['minArea'], deblendCont=SN13038['deblendCont'])

if __name__ == "__main__":
    # runSEPIndiviually()
    # main_hst(13038)
    # map(main_hst, [13354, 13354], [25, 26])
    # main_hst(13038, 8)

    #test to try and find SN19282 as a 2x2 in sdss
    # main_sdss(19282, minarea=3)

    # get integers of the SN numbers/names
    names = np.array(ancillary.get_sn_names(), dtype=int)
    # names = np.array([13354], dtype=int)
    # map(main, names)

    # main_sdss(14284)
    map(main_hst, names)
    # block_size = [8]*len(names)
    # map(main_hst, names, block_size)
    # sigma = [3]*len(names)
    # map(main_sdss, names, sigma)