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

def run_sep(sciData, edge=0.0016, withSigma=False):
    """
    this takes an image and smooths it

    # Parameters
    sciData : 

    edge : float
        Pixel value used for galaxy edge definion. If `withSigma` is `True`
        then edge gets treated as a multiple of the background noise. The 
        default of `0.0016` was from the attempt of calculating HST's 26 mag/
        sqr arcsec value. But HST has a zeropoint of ~26 mag so this is a bad
        limit.

    withSigma : bool
        determeins the meaning of edge, either the defulat (`False`) as a pure
        value or when set to `True` as a sigma of the background noise.

    # Returens
    sources : np.ndarray
        the sources, in a structured array, found from sep.
    """
    if withSigma:
        bkg = sep.Background(sciData)
        threshold = edge*bkg.globalrms
    else:
        threshold = edge

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
    # Where, in array, are the objects close to initialGuess
    searchRadius = 200
    centerIDs = []
    for i, x in enumerate(sources['x']):
        if ((x-initialGuess[0])**2 + (sources['y'][i]-initialGuess[1])**2) < searchRadius**2:
            centerIDs.append(i)
    #todo(add error for nothing found)

    # Select largest of the close objects
    # np.argmax() is much better then np.where(a = max(a))
    idx = np.argmax(sources['npix'][centerIDs])
    #todo(add error for too many found)

    host = sources[['npix', 'x', 'y', 'a', 'b', 'theta']][idx]
    return host

def main():
    """
    This is the default method of defining a galaxy
    """
    SNNumber = 2635
    imageFile = 'data/HST - combined/SN{}_combined.fits'.format(SNNumber)

    #import image data
    data = ancillary.import_fits(imageFile, extention=1)[1]

    #smooth image
    stDev = 3
    data = smooth_Image(data, stDev)

    #run sep
    sources = run_sep(data, 0.006)
    np.savetxt('SN2635_sources_thresh0.0065_cont0.1.csv', sources, delimiter=',')
    print sources

    #get "best" object from sep
    host = find_host(sources)
    print host

    #save resutls to file

if __name__ == "__main__":
    main()