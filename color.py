""" color.py -- Calculating the color at a specific location

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-07-20
    Licesed under the MIT License
"""
import numpy as np
import sep

import astropy.stats
from astropy.wcs import WCS

import ancillary
import fractionalRank

def getData(snid):
    """
    """
    #F475W is like sdss-g (the blue-er filter)
    #F625W is like sdss-r (the red-er filter)

    # get data's location on disk
    blueLocation = 'data/HST - renamed/SN{}_F475W_drz.fits'.format(snid)
    redLocation = 'data/HST - renamed/SN{}_F625W_drz.fits'.format(snid)

    # get blue hdu and data
    blueHDU, blueData = ancillary.import_fits(blueLocation, 1)

    # get red hdu and data
    redHDU, redData = ancillary.import_fits(redLocation, 1)

    # get sn location on the sky
    position = fractionalRank.get_SN_HST_coord(snid)

    # get sn location in pixles
    w = WCS(blueHDU[1].header)
    snPixels = w.all_world2pix(position.ra.to(u.deg).value, 
                               position.dec.to(u.deg).value, 0)

    return blueHDU, blueData, redHDU, redData, snPixels

def calculateSNR(blueData, redData, snPixels):
    """
    Blah blah does preperation to simply run [`astropy.stas.signal_to_noise_oir_ccd`](http://docs.astropy.org/en/stable/api/astropy.stats.signal_to_noise_oir_ccd.html)
    # Parameters

    signal : float 
        fed into `source_eps` and should be electrons per second. 

    # Returns
    """
    #Calcualte varriables
    #same for all
    #from [Wide Field Camera 3 Instrument Handbook for Cycle 24]
    #(http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/wfc3_ihb.pdf)
    dark_eps = 0.001944444    #e-/s
    rd = 3.1                  #e-
    npix = 1.0

    #blue variables
    blueT = blueHDU[0].header['EXPTIME']        #s
    blueSource = blueData[snPixels] #or something like this
    blueBkg = sep.Background(blueData)
    blueSky  = blueBkg.back()[snPixels]

    #red variables
    redT = redHDU[0].header['EXPTIME']        #s
    redSource = redData[snPixels] #or something like this
    redBkg = sep.Background(redData)
    redSky  = redBkg.back()[snPixels]

    #use astropy stats to calculate snr                                    
    blueSNR = astropy.stats.signal_to_noise_oir_ccd(blueT, blueSource, blueSky, 
                                                    dark_eps, rd, npix)
    redSNR = astropy.stats.signal_to_noise_oir_ccd(redT, redSource, redSky, 
                                                   dark_eps, rd, npix)

    #do I want to return the souces value?
    return blueSNR, blueSource, redSNR, redSource

def calcuateColor(blueCountRate, redCountRate, location):
    """
    takes red-blue
    # Parameters
    blue

    # Returns
    """
    # get data

    #get SN location

    #mark says I need to take into account the ccd sensitity
    #todo(above)
    color = 2.5*np.log10(blueCountRate/redCountRate)

    return color


def main(snid):
    """
    # Parameters
    snid : int
        The SDSS-II Candidate ID number of the object you want to find the 
        host for.
 
    # Returns
    """
    #skip objects that do not work
    #objects with forground stars nearby

    #import 
    #import data files
    #F475W is like sdss-g (the blue-er filter)
    imageLocation = 'data/HST - renamed/SN{}_F475W_drz.fits'.format(snid)
    hdu_blue, data_blue = ancillary.import_fits(imageLocation, 1)
    #F625W is like sdss-r (the red-er filter)
    imageLocation = 'data/HST - renamed/SN{}_F625W_drz.fits'.format(snid)
    hdu_red, data_red = ancillary.import_fits(imageLocation, 1)

    #Get SN pixel locations
    #todo(test if we get the same thing for both images.)
    position = fractionalRank.get_SN_HST_coord(snid)
    w = WCS(hdu_blue[1].header)
    SNPixels = w.all_world2pix(position.ra.to(u.deg).value, 
                               position.dec.to(u.deg).value, 0)

    #for each filter, calculate S/N
    snr_red = getSNR(snid)
    snr_blue = getSNR(snid)

    #gnerate save settings/locations

    #save S/N

    #if S/N (in each filter) is great enough, calculate color

    #save color
    
if __name__ == __main__:
