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

    sources = sep.extract(sciData, threshold, minarea=50)
    return sources

def find_host(sources):
    """
    a search to find the host galaxy from a list of possible sources

    # Paramenters
    sources: np.ndarray
        an structured array of sources, in the format returned by 
        [`sep`](https://github.com/kbarbary/sep).

    # Returns
    host: np.ndarray
        a list of all the parameters from sep of the determined host galaxy
    """
    host = []
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
    sources = run_sep(data, 0.04)

    #get "best" object from sep
    find_host(sources)

    #save resutls to file

if __name__ == "__main__":
    main()