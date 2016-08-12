""" color.py -- Calculating the color at a specific location

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-07-20
    Licesed under the MIT License
"""
import re

import numpy as np
import sep

import astropy.stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

import ancillary
import fractionalRank

def getData(snid):
    """
    given a SDSS Transient ID, this gets the HST data for that object. Returns
    both the HDU and the data as well as the pixles of the SN.

    # Parameters
    snid : int 
        The SDSS Transient ID given as an int. (To be honest a string would be 
        ok.)

    # Returns
    blueHDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the bluer, F475W, filter) returned from the
        import feature of `astropy`.

    blueData : ndarray
        The science data of the bluer, F475W, filter. From 
        `blueHDU[extention].data`.

    redHDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the redder, F625W, filter) returned from 
        the import feature of `astropy`.

    redData : ndarray
        The science data of the redder, F625W, filter. From 
        `redHDU[extention].data`.

    snPixels : tuple, two ints
        The result of `all_world2pix` on the `WCS` of the `HDU` objects. A 
        `tuple` of two `int`'s of the order $(x,y)$ and should be useable like
        `snValue = data[snPixels]`
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
    yPixel, xPixel = w.all_world2pix(position.ra.to(u.deg).value, 
                               position.dec.to(u.deg).value, 0)
    #all_world2pix is stupid with just one item to search for. I want a tuple of
    # two float and this is the simple way to do that. With out this I get an
    # array of arrays, but the iner arrays seem to not have a length? They are
    # stupid and now floats!
    #todo(these should be `round()` not `float()`, correct?)
    snPixels = (int(np.round(xPixel)), int(np.round(yPixel)))    
    #python2 round produces float, python3 produces int. I like python3!
    #somehow in python3 we need to use numpy's version of round.

    return blueHDU, blueData, redHDU, redData, snPixels

def calculateSNR(blueHDU, blueData, redHDU, redData, snPixels, size=0):
    """
    Does the signal and noise calculations to simply run 
    [`astropy.stas.signal_to_noise_oir_ccd`](http://docs.astropy.org/en/stable/api/astropy.stats.signal_to_noise_oir_ccd.html).
    Returns the SNR and the signal for both the blue and red filters.
    
    # Parameters
    blueHDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the bluer, F475W, filter) returned from the
        import feature of `astropy`.

    blueData : ndarray
        The science data of the bluer, F475W, filter. From 
        `blueHDU[extention].data`. Should be in [electrons/s].

    redHDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the redder, F625W, filter) returned from 
        the import feature of `astropy`.

    redData : ndarray
        The science data of the redder, F625W, filter. From 
        `redHDU[extention].data`. Should be in [electrons/s].

    snPixels : tuple, two ints
        The result of `all_world2pix` on the `WCS` of the `HDU` objects. A 
        `tuple` of two `int`'s of the order $(x,y)$ and should be useable like 
        `snValue = data[snPixels]`

    size : int
        The number of pixels to around `snPixels` that should be accounted for
        in the size of the appeture. By delfault `size = 0`, so only the 
        `snPixel` is used. If `size = 1` then an extra pixel is added in each 
        direction making a 3x3 square with `snPixel` at the center. With `size 
        = 2`, a 5x5 square is used.

    # Returns
    blueSNR : float
        The result of `signal_to_noise_oir_ccd` for `blueData` at `snPixles`.

    blueSource : float
        The value of `blueData` at `snPixles`. Should be in [electrons/s].

    redSNR : float
        The result of `signal_to_noise_oir_ccd` for `redData` at `snPixles`.

    redSource : float
        The value of `redData` at `snPixles`. Should be in [electrons/s].
    """
    #Calcualte varriables
    #same for all
    #from [Wide Field Camera 3 Instrument Handbook for Cycle 24]
    #(http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/wfc3_ihb.pdf)
    #todo(update to ACS/WFC)
    dark_eps = 0.001944444    #e-/s
    rd = 3.1                  #e-
    npix = (2.0*size + 1.0)**2

    #blue variables
    blueT = blueHDU[0].header['EXPTIME']        #s
    blueSource = blueData[snPixels[0]-size:snPixels[0]+size+1,
                          snPixels[1]-size:snPixels[1]+size+1].sum()
    blueBkg = sep.Background(blueData)
    #sky needs to be for one pixel signal_to_noise_oir_ccd multiplies by npix
    blueSky  = blueBkg.back()[snPixels[0]-size:snPixels[0]+size+1,
                              snPixels[1]-size:snPixels[1]+size+1].mean()

    #red variables
    redT = redHDU[0].header['EXPTIME']        #s
    redSource = redData[snPixels[0]-size:snPixels[0]+size+1,
                        snPixels[1]-size:snPixels[1]+size+1].sum()
    redBkg = sep.Background(redData)
    #sky needs to be for one pixel signal_to_noise_oir_ccd multiplies by npix
    redSky  = redBkg.back()[snPixels[0]-size:snPixels[0]+size+1,
                            snPixels[1]-size:snPixels[1]+size+1].mean()
    # redSky  = redBkg.back()[snPixels]                            

    #use astropy stats to calculate snr                                    
    blueSNR = astropy.stats.signal_to_noise_oir_ccd(blueT, blueSource, blueSky, 
                                                    dark_eps, rd, npix)
    redSNR = astropy.stats.signal_to_noise_oir_ccd(redT, redSource, redSky, 
                                                   dark_eps, rd, npix)

    return blueSNR, blueSource, redSNR, redSource

def calcuateColor(blueCountRate, redCountRate):
    """
    Calcuates the color for two count rate objects. Takes blue-red or more
    exactly: $2.5 \times \log{\frac{red}{blue}}$. This is from the definition of
    [magnitude](https://en.wikipedia.org/wiki/Magnitude_%28astronomy%29).

    # Parameters
    blueCountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s]. This is assumed to be
        the F475W HST filter.

    redCountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s]. This is assumed to be
        the F625W HST filter.

    # Returns
    color : float
        The blue mag minus red mag from the count rates given. In [mag]. From
        the AB system.
    """
    # Content and idea from http://www.stsci.edu/hst/wfc3/phot_zp_lbn
    #todp(are these a r=10 radius apature)

    # blue header info
    inverseSensitivity475W = 1.8267533E-19     #ergs/cm2/Ang/electron                    
    PhotPlam475W = 4.7456221E+03               #Pivot wavelength (Angstroms) 
    ABMagZpt475W = -2.5*np.log10(inverseSensitivity475W) - 5*np.log10(PhotPlam475W)-2.408
    #from http://www.stsci.edu/hst/acs/analysis/zeropoints

    # red header info
    inverseSensitivity625W = 1.1931234E-19     #ergs/cm2/Ang/electron 
    PhotPlam625W = 6.3120513E+03               #Pivot wavelength (Angstroms)                      
    ABMagZpt625W = -2.5*np.log10(inverseSensitivity625W) - 5*np.log10(PhotPlam625W)-2.408
    #from http://www.stsci.edu/hst/acs/analysis/zeropoints

    #Calcualate AB mag
    #get initial converstion 
    blueMag = -2.5*np.log10(blueCountRate) + ABMagZpt475W
    redMag = -2.5*np.log10(redCountRate) + ABMagZpt625W
    #account for resoved source needing to be mag/sqr-arcsec
    pixelScale = 0.05     #arcsec/pixel
    scaling = 1.0/pixelScale**2
    blueMag = blueMag - 2.5*np.log10(scaling)
    redMag = redMag - 2.5*np.log10(scaling)

    #calculate color
    color = blueMag - redMag

    #######
    #Not needed, alternative calcuation methods
    #for SN1415, these calcuate the same thing. 
    # PhotZpt475W = -2.1100000E+01               #ST magnitude zero point 
    # PhotZpt625W = -2.1100000E+01               #ST magnitude zero point 
    # blueMag = -2.5*np.log10(blueCountRate) + PhotZpt475W
    # redMag = -2.5*np.log10(redCountRate) + PhotZpt625W

    # colorCount = 2.5*np.log10(redCountRate/blueCountRate)
    # colorST = blueMag - redMag
    #######
    return blueMag, redMag, color

def getSDSSColor(snID):
    """
    """
    #get snID in form of sdss data
    sn = str(snID).zfill(6)
    
    #read second line from sdss photometry data
    with open('data/SDSS - photometry/SMP_{0}.dat'.format(sn), 'r') as f:
        f.readline()
        f.readline()
        third_line = f.readline()

    #clean up second line and extract desired info
    split = np.array( re.split(r'\s+', third_line) )     #split on whitespace
    #the index of the g (r) values are 3 (4) spaces over from the first '='
    gIndex = np.where(split == '=')[0]+3
    rIndex = gIndex + 1
    
    #calcualte color
    #todo(accoutn for units of asinh-mag/square-arcsec)
    # S = m +2.5log(Area), m = S - 2.5log(Area)
    gSB, rSB = float(split[gIndex][0]), float(split[rIndex][0])
    # print('split: ', split)
    # print('index: ', gIndex, rIndex)
    gmag = -2.5*np.log10(gSB*1e-6/3631)
    rmag = -2.5*np.log10(rSB*1e-6/3631)
    # print('mag: ', gmag, rmag)
    # from sys import exit; exit()
    #g-r is -2.5log(f_g/f_r)
    color = -2.5*np.log10(gSB/rSB)

    return gmag, rmag, color

# def saveData(snid, blueSNR, blueSource, redSNR, redSource, color, sdssColor):
def saveData(*args):
    """
    Saves the input data as an appeneded line to `'resources/hst_color.csv'`. 
    The file needs to exists. Puting header info would be good, such as:

    ```
    #The SNR and color analysis of local environment around SDSS SNIa viewed by HST                 
    #snid, blueSNR, blueSource, redSNR, redSource, color, sdss color
    #    ,        , [counts/s],     , [counts/s], F475W-F625W [mag], g-r [mag]
    ```

    # Parameters
    snid : int 
        The SDSS Transient ID given as an int. (To be honest a string would be 
        ok.) 

    blueSNR : float
        The result of `signal_to_noise_oir_ccd` for `blueData` at `snPixles`.

    blueSource : float
        The value of `blueData` at `snPixles`. Should be in [electrons/s].

    redSNR : float
        The result of `signal_to_noise_oir_ccd` for `redData` at `snPixles`.

    redSource : float
        The value of `redData` at `snPixles`. Should be in [electrons/s].

    color : float
        The blue mag minus red mag from the count rates given. In [mag].

    sdssColor : float
        The g-r (in mag) of from the SDSS data
    """
    #Create saving location
    #`defGalaxy.saveGalaxies()` has a more robust way to doing this.
    saveFileName = 'resources/hst_color.csv'

    #combine the data to make it writeable
    #this chould be done with `np.array()` but this also works. 
    #need double list because `np.savetxt` removes one. This makes `toWrite`
    #a single row as desired.
    #using `args` might be better?
    # toWrite = np.stack([[snid, blueSNR, blueSource, redSNR, redSource, 
                         # color, sdssColor]])
    toWrite = [args]
    
    #thanks to stackoverflow.com/questions/27786868/
    #python3-numpy-appending-to-a-file-using-numpy-savetxt
    #open the file to `append` and in `binary` (for python3) mode.
    with open(saveFileName,'ab') as f:
        np.savetxt(f, toWrite, delimiter=',')

def main(snid, snr=2):
    """
    # Parameters
    snid : int
        The SDSS-II Candidate ID number of the object you want to find the 
        host for.

    snr : float
        The SNR cut off for calcualting the HST color.
    """

    #Get Data
    blueHDU, blueData, redHDU, redData, snPixels = getData(snid)

    #Calculate SNR
    size = 4     #this is a 2*4+1 or 9 per side. 
                 #this gets us a bit larger then sdss (0.4 arcsec/pixel) form
                 #hst's (0.05 arcsec/pixel)
    blueSNR, blueSource, redSNR, redSource = calculateSNR(blueHDU, blueData, 
                                                          redHDU, redData, 
                                                          snPixels, size)

    #Calculate Color
    if blueSNR > snr and redSNR > snr:
        blueMag, redMag, color = calcuateColor(blueSource, redSource)
    else:
        blueMag, redMag, color = np.nan, np.nan, np.nan

    sdssG, sdssR, sdssColor = getSDSSColor(snid)

    #Save results 
    saveData(snid, blueSNR, blueSource, blueMag, redSNR, redSource, redMag,
             color, sdssG, sdssR, sdssColor)
    
if __name__ == '__main__':
    # main(20874)
    # main(1415)
    # main(14279, 5.0)

    names = np.array(ancillary.get_sn_names(), dtype=int)
    snr = np.ones(names.shape)*18.0
    list(map(main, names, snr))