""" defGalaxy.py -- set of functions to define galaxy shape size to 
    `resources/galaxies_{0}.csv` where `{0}` is the edge definition.

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-03
    Licesed under the MIT License
"""
from sys import exit
from datetime import datetime
from os import path, makedirs

import numpy as np
import sep

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits     #to read fits files

import ancillary 

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

def run_sep(sciData, threshold):
    """
    this takes an image and smooths it

    # Parameters
    sciData : 

    threshold : float
        Pixel value used for galaxy edge definion. 

    # Returens
    sources : np.ndarray
        the sources, in a structured array, found from sep.
    """
    bkg = sep.Background(sciData)
    bkg.subfrom(sciData)

    sources = sep.extract(sciData, threshold, minarea=50, deblend_cont=0.1)
    return sources

def find_host(sources, initialGuess = (2090/2.0, 2108/2.0)):
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
    searchRadius = 200
    centerIDs = []

    for i, x in enumerate(sources['x']):
        if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
            centerIDs.append(i)
    #todo(add error for nothing found)
    if len(centerIDs) == 0:
            searchRadius = 500
            centerIDs = []    #just play it safe

            for i, x in enumerate(sources['x']):
                if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
                    centerIDs.append(i)
    
     # Select largest of the center objects, but save the ID
    #this selects the centerID associated with the max (in size) of the central cources
    idx = centerIDs[np.argmax(sources['npix'][centerIDs])]
    #todo(add error for too many found)

    host = sources[['npix', 'x', 'y', 'a', 'b', 'theta']][idx]
    return host

def saveGalaxies(sources, host, SNNumber, sb):
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
        The surface brightness threshold that the soucers were found at, in mag/sqr-arcsec
    """

    #Save ALL objects
    allHeader = 'data from sep on SN{0} with a SB cutoff of {1} mag/sqr-arcsec on {2}'.format(SNNumber, sb, datetime.now().date()) + '\n' + str(sources.dtype.names)[1:-1].replace("'","")
    allLocationFolder = 'resources/SN{}'.format(SNNumber)
    allLocationFIle = 'sources_{}.csv'.format(sb)
    allLocation = allLocationFolder+'/'+allLocationFIle
    if not path.exists(allLocationFolder): makedirs(allLocationFolder)

    np.savetxt(allLocation, sources, delimiter=',', header=allHeader)

    
    #Save host
    #todo(don't assume it existes, but rather check.)
    hostLocation = 'resources/hosts_{1}.csv'.format(SNNumber, sb)
    #combine `SNNumber` to the front of `host` data, but structured arrays are stupid. So is np.savetxt()
    dataToSave = str(SNNumber)
    for i in host:
        dataToSave += ', ' + str(i)
    dataToSave += '\n'

    with open(hostLocation, 'a') as f:
        f.write(dataToSave)

def main(SNNumber = 2635):
    """
    This is the default method of defining a galaxy
    """
    print "running SN{}".format(SNNumber)
    # imageFile = 'data/HST - combined/SN{}_combined.fits'.format(SNNumber)
    imageFile = 'data/HST - combined/SN{}_combined_flux.fits'.format(SNNumber)

    #import image data
    data = ancillary.import_fits(imageFile, extention=1)[1]

    #smooth image
    stDev = 3
    data = smooth_Image(data, stDev)

    # find threshold
    surfaceBrightness = 25    #mag/sqr-arcsec
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
    # print 'sources: ', sources[['npix', 'x', 'y', 'a', 'b', 'theta']]

    #get "best" object from sep
    host = find_host(sources)
    # print host

    #save resutls to file
    saveGalaxies(sources, host, SNNumber, surfaceBrightness)

if __name__ == "__main__":
    # get integers of the SN numbers/names
    names = np.array(ancillary.get_sn_names(), dtype=int)
    map(main, names)
    # main()

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

