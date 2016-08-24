""" color.py -- Calculating the color at a specific location

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-07-20
    Licesed under the MIT License
"""
import re
import warnings

import numpy as np
import sep
import pandas as pd

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
    F475HDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the bluer, F475W, filter) returned from the
        import feature of `astropy`.

    F475Data : ndarray
        The science data of the bluer, F475W, filter. From 
        `blueHDU[extention].data`.

    F625HDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the redder, F625W, filter) returned from 
        the import feature of `astropy`.

    F625Data : ndarray
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
    F475Location = 'data/HST - renamed/SN{}_F475W_drz.fits'.format(snid)
    F625Location = 'data/HST - renamed/SN{}_F625W_drz.fits'.format(snid)

    # get blue hdu and data
    F475HDU, F475Data = ancillary.import_fits(F475Location, 1)

    # get red hdu and data
    F625HDU, F625Data = ancillary.import_fits(F625Location, 1)

    # get sn location on the sky
    position = fractionalRank.get_SN_HST_coord(snid)

    # get sn location in pixles
    w = WCS(F475HDU[1].header)
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

    return F475HDU, F475Data, F625HDU, F625Data, snPixels

def calculateSNR(F475HDU, F475Data, F625HDU, F625Data, snPixels, size=0):
    """
    Does the signal and noise calculations to simply run 
    [`astropy.stas.signal_to_noise_oir_ccd`](http://docs.astropy.org/en/stable/api/astropy.stats.signal_to_noise_oir_ccd.html).
    Returns the SNR and the signal for both the blue and red filters.
    
    # Parameters
    F475HDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the bluer, F475W, filter) returned from the
        import feature of `astropy`.

    F475Data : ndarray
        The science data of the bluer, F475W, filter. From 
        `F475HDU[extention].data`. Should be in [electrons/s].

    F625HDU : astropy.io.fits.HDUList
        The resulting `HDUList` (for the redder, F625W, filter) returned from 
        the import feature of `astropy`.

    F625Data : ndarray
        The science data of the redder, F625W, filter. From 
        `F625HDU[extention].data`. Should be in [electrons/s].

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
    F475SNR : float
        The result of `signal_to_noise_oir_ccd` for `F475HDU` at `snPixles`.

    F475Source : float
        The value of `F475HDU` at `snPixles`. Should be in [electrons/s].

    F625SNR : float
        The result of `signal_to_noise_oir_ccd` for `F625Data` at `snPixles`.

    F625ource : float
        The value of `F625Data` at `snPixles`. Should be in [electrons/s].
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
    F475T = F475HDU[0].header['EXPTIME']        #s
    F475Source = F475Data[snPixels[0]-size:snPixels[0]+size+1,
                          snPixels[1]-size:snPixels[1]+size+1].sum()
    F475Bkg = sep.Background(F475Data)
    #sky needs to be for one pixel signal_to_noise_oir_ccd multiplies by npix
    F475Sky  = F475Bkg.back()[snPixels[0]-size:snPixels[0]+size+1,
                              snPixels[1]-size:snPixels[1]+size+1].mean()

    #red variables
    F625T = F625HDU[0].header['EXPTIME']        #s
    F625Source = F625Data[snPixels[0]-size:snPixels[0]+size+1,
                        snPixels[1]-size:snPixels[1]+size+1].sum()
    F625Bkg = sep.Background(F625Data)
    #sky needs to be for one pixel signal_to_noise_oir_ccd multiplies by npix
    F625Sky  = F625Bkg.back()[snPixels[0]-size:snPixels[0]+size+1,
                            snPixels[1]-size:snPixels[1]+size+1].mean()
    # F625Sky  = F625Bkg.back()[snPixels]                            

    #use astropy stats to calculate snr                                    
    F475SNR = astropy.stats.signal_to_noise_oir_ccd(F475T, F475Source, F475Sky, 
                                                    dark_eps, rd, npix)
    F625SNR = astropy.stats.signal_to_noise_oir_ccd(F625T, F625Source, F625Sky, 
                                                   dark_eps, rd, npix)

    return F475SNR, F475Source, F625SNR, F625Source

def calcuateColor(F475CountRate, F625CountRate, scale):
    """
    Calcuates the color for two count rate objects. Takes blue-red or more
    exactly: $2.5 \times \log{\frac{red}{blue}}$. This is from the definition of
    [magnitude](https://en.wikipedia.org/wiki/Magnitude_%28astronomy%29).

    # Parameters
    F475CountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s]. This is assumed to be
        the F475W HST filter, and the bluer of the two inputs.

    F625CountRate : float
        The value, from the science data, at the location where the color
        caluation is desired. Should be in [electrons/s]. This is assumed to be
        the F625W HST filterand the reder of the two inputs.

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
    F475Mag = -2.5*np.log10(F475CountRate) + ABMagZpt475W
    F625Mag = -2.5*np.log10(F625CountRate) + ABMagZpt625W
    #account for resoved source needing to be mag/sqr-arcsec & scaling number of pixels
    pixelScale = 0.05                #arcsec/pixel
    npix = (2.0*scale + 1.0)**2      #pixels**2
    area = npix/pixelScale**2        #arcsec**2
    F475Mag = F475Mag - 2.5*np.log10(area)
    F625Mag = F625Mag - 2.5*np.log10(area)

    #calculate color
    color = F475Mag - F625Mag

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

    #calucate uncertatinty
    return F475Mag, F625Mag, color

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
    if gSB < 0 or rSB < 0:
        #Abort, we got a negative flux!
        warnings.warn(r'SN{} has a negative flux: g = {} μJy/sqr-arcsec, r = {} μJy/sqr-arcsec,'.format(snID, gSB, rSB))
        gMag, gMagUncert, rMag, rMagUncert, color, colorUncert = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan 
    else:
        #todo(do this correct)
        gUncertIndex = gIndex + 8
        rUncertIndex = rIndex + 8
        uncertG, uncertR = float(split[gUncertIndex][0]), float(split[rUncertIndex][0])
        
        gMag = -2.5*np.log10(gSB*1e-6/3631)
        rMag = -2.5*np.log10(rSB*1e-6/3631)
        gMagUncert = uncertG*np.abs(2.5 / (gSB*np.log(10)))
        rMagUncert = uncertR*np.abs(2.5 / (rSB*np.log(10)))
        # print('mag: ', gMag, rMag)
        # from sys import exit; exit()
        #g-r is -2.5log(f_g/f_r)
        color = -2.5*np.log10(gSB/rSB)
        #magnitude error is like a percet error, so add in quadriture
        colorUncert = np.sqrt(gMagUncert**2 + rMagUncert**2)

    return gMag, gMagUncert, rMag, rMagUncert, color, colorUncert

# def saveData(snid, blueSNR, blueSource, redSNR, redSource, color, sdssColor):
def saveData(*args):
    """
    Saves the input data as an appeneded line to `'resources/hst_color.csv'`. 
    The file needs to exists. Puting header info would be good, such as:

    ```
    #The SNR and color analysis of local environment around SDSS SNIa viewed by HST                 
    #snid, sdss g mag, sdss g uncert, sdss r mag, sdss r uncert, sdss color, sdss c olor uncert,   hst size, F475W mag, F475W SNR, F625W mag, F625W SNR, hst color, hst color uncert, use, color, color uncert
    #    ,      [mag],         [mag],      [mag],         [mag],      [mag],              [mag], [0 1 2 or 3],     [mag],          ,     [mag],          ,     [mag],            [mag],    , [mag],        [mag]
    ```

    Explantions and details of these parameters can be seen below.

    # Parameters
    snid : int 
        The SDSS Transient ID given as an int. (To be honest a string would be 
        ok.) 

    sdssG : 

    sdssGUncert : 

    sdssR :
    
    sdssrUncert :

    sdssColor :
        The g-r (in mag) of from the SDSS data

    sdssColorUncert :

    size :

    F475Mag :

    F475SNR :
        The result of `signal_to_noise_oir_ccd` for `F475Data` at `snPixles`.

    F625Mag :

    F625SNR :
        The result of `signal_to_noise_oir_ccd` for `F625Data` at `snPixles`.

    hstColor :
        The F475W (blue) mag minus F625W (red) mag from the count rates given. 
        In [mag].

    hstColorUncert :

    use :

    color :

    colorUncert :
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
    #convert to dataframe, becasue `np.savetxt` appears to not be able to take strings
    toWrite = pd.DataFrame([args])

    
    #thanks to stackoverflow.com/questions/27786868/
    #python3-numpy-appending-to-a-file-using-numpy-savetxt
    #open the file to `append` and in `binary` (for python3) mode.
    # with open(saveFileName,'ab') as f:
    #     np.savetxt(f, toWrite, delimiter=',')
    # `binary` might not like strings.
    toWrite.to_csv(saveFileName, mode='a', header=False, index=False, na_rep='nan')

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
    F475HDU, F475Data, F625HDU, F625Data, snPixels = getData(snid)

    #get SDSS color for referance quality
    sdssG, sdssGUncert, sdssR, sdssRUncert, sdssColor, sdssColorUncert = getSDSSColor(snid)

    #Calculate SNR
    #set defualts to be sdss, will be overwritten if hst is better
    use = 'sdss'
    color, colorUncert = sdssColor, sdssColorUncert
    F475Mag, F625Mag, hstColor, hstColorUncert = np.nan, np.nan, np.nan, np.nan
    #size is a 2*n+1 or 9x9 for n=4. 
    #n=4 overscales hst (0.05 arcsec/pixel) more than sdss (0.4 arcsec/pixel)
    #perfect scaling is 8x8. 
    #Don't search for HST scaled by 4, that is not going to be better then SDSS
    
    for size in np.arange(4):
        (F475SNR, F475Source, F625SNR, 
            F625Source) = calculateSNR(F475HDU, F475Data, F625HDU, F625Data, 
                                      snPixels, size)
        if F475Source < 0 or F625Source <0:
            #continue if source is negative in either filter
            continue
        #flag what color to use, exit if quality if better
        #todo(or we can have the error just be better then the SDSS error?)
        #todo(do not let negative SNR rates through.)
        #1.08 is the difference between mag and fractional uncertanties.
        #note that F475W seems to always have the lowest SNR, just cause?
        if 1.08/F475SNR < 0.1 and F475SNR > 0:
            #Calculate HST Color
            F475Mag, F625Mag, hstColor = calcuateColor(F475Source, F625Source, 
                                                      size)
            # I should do this, http://spiff.rit.edu/classes/phys445/lectures/signal/signal_illus.html
            # it says that "the uncertainty in magnitudes will be 1.08 times the
            # fractional uncertainty in brightness" with the fractional = 1/SNR.
            #1.08 is the difference between mag and fractional uncertanties.
            hstColorUncert = np.sqrt((1.08/F475SNR)**2 + (1.08/F625SNR)**2)
            use = 'hst'
            color, colorUncert = hstColor, hstColorUncert
            break


    #Save results, Outlined in 2016-08-19 lab notebook.
    saveData(snid, sdssG, sdssGUncert, sdssR, sdssRUncert, sdssColor, 
             sdssColorUncert, size, F475Mag, F475SNR, F625Mag, F625SNR, 
             hstColor, hstColorUncert, use, color, colorUncert)
    # #OLD DONT USE# saveData(snid, F475SNR, F475Source, F475Mag, F625SNR, F625Source, F625Mag,
             # color, sdssG, sdssR, sdssColor, sdssColorUncert, use, size)
    # print(snid, sdssColor, sdssColorUncert, size, hstColor, hstColorUncert, use, color, colorUncert)
    # print(snid, sdssG, sdssGUncert, sdssR, sdssRUncert, sdssColor, 
    #          sdssColorUncert, size, F475Mag, F475SNR, F625Mag, F625SNR, 
    #          hstColor, hstColorUncert, use, color, colorUncert)
    
if __name__ == '__main__':
    # main(20874)
    # main(14284)
    # main(14279, 5.0)

    names = np.array(ancillary.get_sn_names(), dtype=int)
    snr = np.ones(names.shape)*18.0
    list(map(main, names, snr))