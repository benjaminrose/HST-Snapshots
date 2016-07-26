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
    snPixels = (np.round(xPixel), np.round(yPixel))    
    #python2 round produces float, python3 produces int. I like python3!
    #somehow in python3 we need to use numpy's version of round.

    return blueHDU, blueData, redHDU, redData, snPixels

def calculateSNR(blueHDU, blueData, redHDU, redData, snPixels):
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
    dark_eps = 0.001944444    #e-/s
    rd = 3.1                  #e-
    npix = 1.0

    #blue variables
    blueT = blueHDU[0].header['EXPTIME']        #s
    blueSource = blueData[snPixels] #or something like this
    blueBkg = sep.Background(blueData)
    blueSky  = blueBkg.back()[snPixels]
    
    print('blue ', blueT, blueSource, blueSky)

    #red variables
    redT = redHDU[0].header['EXPTIME']        #s
    redSource = redData[snPixels] #or something like this
    redBkg = sep.Background(redData)
    redSky  = redBkg.back()[snPixels]

    print('red ', redT, redSource, redSky)

    #use astropy stats to calculate snr                                    
    blueSNR = astropy.stats.signal_to_noise_oir_ccd(blueT, blueSource, blueSky, 
                                                    dark_eps, rd, npix)
    redSNR = astropy.stats.signal_to_noise_oir_ccd(redT, redSource, redSky, 
                                                   dark_eps, rd, npix)

    print('snr ', blueSNR, redSNR)

    return blueSNR, blueSource, redSNR, redSource

def calcuateColor(blueCountRate, redCountRate):
    """
    Calcuates the color for two count rate objects. Takes red-blue or more
    exactly: $2.5 \times \log{\frac{blue}{red}}$. This is from the definition of
    [magnitude](https://en.wikipedia.org/wiki/Magnitude_%28astronomy%29).

    # Parameters
    blueCountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s].

    redCountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s].

    # Returns
    color : float
        The red mag minus blue mag from the count rates given. In [mag].
    """
    #mark says I need to take into account the ccd sensitity
    #todo(above)
    color = 2.5*np.log10(blueCountRate/redCountRate)

    return color

def saveData(snid, blueSNR, blueSource, redSNR, redSource, color):
    """
    Saves the input data as an appeneded line to `'resources/hst_color.csv'`. 
    The file needs to exists. Puting header info would be good, such as:

    ```
    #The SNR and color analysis of local environment around SDSS SNIa viewed by HST                 
    #snid, blueSNR, blueSource, redSNR, redSource, color
    #    ,        , [counts/s],       , [counts/s], F625W-F475W [mag]
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
        The red mag minus blue mag from the count rates given. In [mag].
    """
    #Create saving location
    #`defGalaxy.saveGalaxies()` has a more robust way to doing this.
    saveFileName = 'resources/hst_color.csv'

    #combine the data to make it writeable
    #this chould be done with `np.array()` but this also works. 
    #need double list because `np.savetxt` removes one. This makes `toWrite`
    #a single row as desired.
    toWrite = np.stack([[snid, blueSNR, blueSource, redSNR, redSource, color]])
    
    #thanks to stackoverflow.com/questions/27786868/
    #python3-numpy-appending-to-a-file-using-numpy-savetxt
    #open the file to `append` and in `binary` (for python3) mode.
    with open(saveFileName,'ab') as f:
        np.savetxt(f, toWrite, delimiter=',')

def main(snid):
    """
    # Parameters
    snid : int
        The SDSS-II Candidate ID number of the object you want to find the 
        host for.
    """

    #Get Data
    blueHDU, blueData, redHDU, redData, snPixels = getData(snid)
    print(type(tuple(snPixels)))
    #Calculate SNR
    blueSNR, blueSource, redSNR, redSource = calculateSNR(blueHDU, blueData, 
                                                          redHDU, redData, 
                                                          snPixels)

    #Calculate Color
    if blueSNR > 1.0 and redSNR > 1.0:
        color = calcuateColor(blueSource, redSource)
    else:
        color = np.nan
    print(color)

    #Save results 
    saveData(snid, blueSNR, blueSource, redSNR, redSource, color)
    
if __name__ == '__main__':
    # main(20874)

    names = np.array(ancillary.get_sn_names(), dtype=int)
    list(map(main, names))