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
    this takes an image and smooths it

    # Parameters
    sciData : np.ndarray

    stDev : int
        The standard deviation used as a parameter in 
        `astropy.convolution.Gaussian2DKernel`.


    # Returens
    smoothedData : 

    """
    gauss = Gaussian2DKernel(stDev)
    smoothedData = convolve(sciData, gauss) 
    return smoothedData

def run_sep(sciData, threshold, minarea=50):
    """
    this takes an image and smooths it

    # Parameters
    sciData : 

    threshold : float
        Pixel value used for galaxy edge definion. 

    minarea : int
        The minimum area of the source. Passed directly to sep.

    # Returens
    sources : np.ndarray
        the sources, in a structured array, found from sep.
    """
    bkg = sep.Background(sciData)
    bkg.subfrom(sciData)

    sources = sep.extract(sciData, threshold, minarea=minarea, deblend_cont=0.1)
    return sources

def find_host(sources, initialGuess=(2090/2.0, 2108/2.0), searchRadius=200):
    """
    a search to find the host galaxy from a list of possible sources

    # Issues
        This is failing with 2635, it is getting an object at (1108.1448309503585, 37.95088945611436) instead of at (1027.0670854882237, 1049.7834967223419)

    # Paramenters
    sources : np.ndarray
        an structured array of sources, in the format returned by 
        [`sep`](https://github.com/kbarbary/sep). Need values 
        `['npix', 'x', 'y', 'a', 'b', 'theta']`.

    initialGuess : tuple
        The likely pixel location of the host. Needs to be of length 2. This 
        defaults to the center of an HST image but the SN's location can be 
        used instead.

    # Returns
    host : np.void
        A list of all the parameters from sep of the determined host galaxy.
        It still a structured array that contains fields named:
        `['npix', 'x', 'y', 'a', 'b', 'theta']`
    """
    #todo(this fails. It is not selecting the right thing at all.)
    # Where, in array, are the objects close to initialGuess

    #todo(change so it searches to 200 always but does a while loop till it finds something.)
    centerIDs = []

    for i, x in enumerate(sources['x']):
        if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
            centerIDs.append(i)
    #todo(add error for nothing found)
    if len(centerIDs) == 0:
        warnings.warn("This SN can't be found in initial search")
        searchRadius = 500
        centerIDs = []    #just play it safe

        for i, x in enumerate(sources['x']):
            if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
                centerIDs.append(i)
    
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

def saveGalaxies(sources, host, SNNumber, sb, telescope='hst'):
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
    """

    #Save ALL objects
    allHeader = 'data from sep on SN{} as observed by '.format(SNNumber) + telescope + ' with a SB cutoff of {0} mag/sqr-arcsec on {1}'.format(sb, datetime.now().date()) + '\n' + str(sources.dtype.names)[1:-1].replace("'","")
    allLocationFolder = 'resources/SN{}'.format(SNNumber)
    allLocationFIle = telescope+'_sources_{}.csv'.format( sb)
    allLocation = allLocationFolder+'/'+allLocationFIle
    if not path.exists(allLocationFolder): makedirs(allLocationFolder)

    np.savetxt(allLocation, sources, delimiter=',', header=allHeader)

    
    #Save host
    #todo(don't assume it existes, but rather check.)
    hostLocation = 'resources/' + telescope + '_hosts_{}.csv'.format(sb)
    #combine `SNNumber` to the front of `host` data, but structured arrays are stupid. So is np.savetxt()
    dataToSave = str(SNNumber)
    for i in host:
        dataToSave += ', ' + str(i)
    dataToSave += '\n'

    with open(hostLocation, 'a') as f:
        f.write(dataToSave)

def main_hst(SNNumber = 2635, surfaceBrightness = 25):
    """
    This is the default method of defining an hst galaxy

    # Parameters
    SNNumber : int
        The sdss identification number of the object you want to find the host 
        for.

    surfaceBrightness : int, float
        the disired cut off in mag/sqr-arcsec
    """
    print("running SN{}".format(SNNumber))
    # imageFile = 'data/HST - combined/SN{}_combined.fits'.format(SNNumber)
    imageFile = 'data/HST - combined/SN{}_combined_flux.fits'.format(SNNumber)

    #import image data
    data = ancillary.import_fits(imageFile, extention=1)[1]

    #smooth image
    stDev = 3
    data = smooth_Image(data, stDev)

    # find threshold
    #the magniutde per pixel for the same SB, wikipedia for equation
    magPerPixel = surfaceBrightness - 2.5*np.log10(0.04**2)
    #convert to AB-mag to flux
    #more at http://www.stsci.edu/hst/acs/analysis/zeropoints
    fNu = 3.63e-20*(u.erg / u.cm**2 / u.s / u.Hz)     #this is when AB-mag = 0
    #convert from f_nu to f_lambda
    f475 = fNu.to(u.erg / u.cm**2 / u.s / u.nm, equivalencies=u.spectral_density(475*u.nm))    #todo(maybe use the pivot wavelenght instead)
    f625 = fNu.to(u.erg / u.cm**2 / u.s / u.nm, equivalencies=u.spectral_density(625*u.nm))
    #this is the zero-mag flux combined like the data
    fluxAB = f475 + f625
    fluxThresh = fluxAB * 10**(magPerPixel/-2.5)

    #run sep
    sources = run_sep(data, fluxThresh.value)
    # print('sources: ', sources[['npix', 'x', 'y', 'a', 'b', 'theta']])

    #get "best" object from sep - 
    #or select from hardcoded objects that do not work in the test.
    hardcoded = [8297, 13354, 13411, 14113, 14284, 18415]
    if SNNumber in hardcoded:
        # if SN is hardcoded do the search, but with a very small radius and a guess that is what is determined between the R_25 and R_26 runs that were visually inspected.
        warnings.warn("SN{}'s host is searched via a hard-coded method".format(SNNumber))
        radius = 20     #pixels
        if SNNumber == 8297:
            # Does not work for R_25, and R_26 is just fine.
            # SNPixels = (1111, 1025)
            # host = find_host(sources, SNPixels, radius)
            host = find_host(sources)
        elif SNNumber == 13354:
            SNPixels = (978, 1043)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 13411:
            SNPixels = (1010, 1067)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 14113:
            SNPixels = (1041, 1049)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 14284:
            SNPixels = (1063, 1056)
            host = find_host(sources, SNPixels, radius)
        elif SNNumber == 18415:
            SNPixels = (1095, 1066)
            host = find_host(sources, SNPixels, radius)
        else:
            raise NotImplementedError('Somehow SN{} is designated for HST hardcoding but is not implemented'.format(SNNumber))
    else:
        host = find_host(sources)
    # print('host: ', host)

    #save resutls to file
    saveGalaxies(sources, host, SNNumber, surfaceBrightness)

def main_sdss(SNNumber = 2635, fltr='g'):
    """
    This is the default method of defining an sdss galaxy

    # Parameters
    SNNumber : int
        The sdss identification number of the object you want to find the host 
        for.

    fltr : string, char
        The single character of an SDSS filter: u, g, r, i, or z (sill python 
        with `filter` being a built in function!)
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
    sigma = 2
    bkg = sep.Background(data)
    thresh = sigma*bkg.globalrms
    #run sep
    sources = run_sep(data, thresh, 10)
    np.savetxt('temp-sn2635-sdss-sources.csv', sources, delimiter=',', header='temp sn2635 sdss sources')
    print('sources: ', sources[['npix', 'x', 'y', 'a', 'b', 'theta']])

    #get "best" object from sep - 
    #or select from hardcoded objects that do not work in the test.
    hardcoded = [0]
    if SNNumber in hardcoded:
        # hardcoed in the initial guess and search radius, the sdss method did not work on these.
        if SNNumber == 0:
            host = find_host(sources, SNPixels, radius)
        else:
            raise NotImplementedError('Somehow SN{} is designated for SDSS hardcoding but is not implemented'.format(SNNumber))
    else:
        host = find_sdss_host(sources, SNNumber, hdu)
    print('host: ', host)

    #save resutls to file
    saveGalaxies(sources, host, SNNumber, sigma, 'sdss')


if __name__ == "__main__":
    # main_hst(18415, 25)
    map(main_hst, [13354, 13354], [25, 26])

    # get integers of the SN numbers/names
    # names = np.array(ancillary.get_sn_names(), dtype=int)
    # map(main, names)

    # main_sdss()
    # map(main_hst, names)

    '''
    ## Get SN number
    SNNumber = 2635
    ### get SDSS and HST image file names
    sdssFile = 'data/SDSS - coadded/SN{}-g.fits'.format(SNNumber)
    #todo(this should propbably be combined rather then just g or r)
    hstFile = 'data/HST - combined/SN{}_combined.fits'.format(SNNumber)
    ### Do I need the location of the SN to figure out SDSS oject?


    ## Import science data
    #`ancillary.import_fits() returns (hdu, scidata)
    sdssData = ancillary.import_fits(sdssFile, extention=0)[1]
    hstData = ancillary.import_fits(hstFile, extention=1)[1]


    ## run sep on SDSS
    ### does this need any treatment?
    sdssSources = run_sep(sdssData, 0.006)


    ## select the most likely hosts


    ## save out full (OR just most likey?) sep results from SDSS
    saveAS = 'SN{}_sdss_sources_thresh0.0065_cont0.1.csv'.format(SNNumber)
    np.savetxt(saveAs, sources, delimiter=',')


    ## run sep on HST
    ### still smoothing?
    stDev = 3
    hstData = smooth_Image(hstData, stDev)
    hstSources = run_sep(hstData, 0.006)


    ## select the most likely hosts


    ## save out full (OR just most likey?) sep results from HST
    saveAS = 'SN{}_hst)_sources_thresh0.0065_cont0.1.csv'.format(SNNumber)
    np.savetxt(saveAs, sources, delimiter=',')


    ## Combine and get final galaxy
    ### center from HST
    ### a, b from SDSS
    ### theta from SDSS rotated to HST's wcs
    #### Check to see if `wcs.crval` is corect. It seems that for SN2635 they were too close together.
    '''

