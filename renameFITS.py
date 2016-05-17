""" renameFITS.py -- rename the HST Fits files. Also get added and 
    color versions. Orignal images can be found relative to this file 
    'data/HST - original/*.fits'. Took 54 mins with 122 images on 2015-05-29.
    Took 18 mins to run just coadd on 2016-05-17.

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-05-16
    Licesed under the MIT License
"""
from __future__ import print_function, division
from copy import deepcopy
from sys import exit
import glob
from datetime import datetime

from numpy import shape
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.constants import c, h

import ancillary

def rename(doc):
    """
    * doc -> file name of the image you want to rename
    * takes fits files (*_drz.fits & *_flt.fits) and remnames them to <<SN>>_<<filter>>_<<img>>_<<count>>.fits
        * <<SN>> is the SDSS tranient number
        * <<filter>> is the HST filter (by number)
        * <<img>> is wether its drz or flt
        * <<count>> optional and only for flt, since there are two images per filter.
    * saves renamed as 'data/HST - renamed/*.fits'

    """
    #from file name, is it drz or flt
    img = doc[-8:-5] #assumes file name is *_drz.fits or *_flt.fits

    #open fits file
    hdu = fits.open(doc)                        # open file
    
    # get names
    sn = hdu[0].header['targname']              #might have trailing white space
    filter1 = hdu[0].header['filter1']          #F475W or F625W
    filter2 = hdu[0].header['filter2']          #this should be CLEAR2L
    
    # check that all is ok - could do more
    if filter2 != "CLEAR2L":                            #if filter2 is not clear
        sys.exit('Not using clear filter in ' + doc)    # abort with error message

    # save renamed file
    if img == 'flt':                #If they are flt there should be 2 files, needs more work
        attempt = 1                     #Is this image 1 or 2?
        while True:                                     #make infinate loop
            name = sn+'_'+filter1+'_'+str(attempt)+doc[-9:]     #define name to save images 
            try:
                hdu.writeto('data/HST - renamed/'+name)         #can I save with this attemt?
                break                                           #if so break
            except IOError:                 #if error
                attempt += 1                #iterate up on attempt (then loop again)
                if attempt > 4:             #if attempt is too high, exit with error message. NO INFINITE LOOPS!
                    sys.exit('Over 3 images of the same SN/filter check '+ doc +' should be ' + name)
    else:                           #for drz images
        name = sn+'_'+filter1+doc[-9:]      #define name to save images 
        hdu.writeto('data/HST - renamed/'+name)     #save fits file
    print('saved '+name)

    # close file to clear memory
    hdu.close()
    return None

def coadd(img1, img2):
    """
    This takes two images (oringally desinged for HST images), and adds them together. The data from the input images are assumed to be in [eletrons/second] and the resulting file will be in [AB-Mag/pixel], maybe? ALso they need to be able to be straight added, no convolution, seeing corrections, or anything fancy like that. comments in new header assume img1 == F475W and img2 == F625W. Should be `drz` images. Saved as `data/<<SN>>_combined.fits`.

    # Parameters
    img1 : string
        The file name of the fits file of the first filter
    img2 : np.ndarray
        The file name of the fits file of the second filter
    """

    ######################
    # open images
    ######################

    #todo(should I use try?)
    #todo(change this whole function to use ancillary)
    if img1[-5:] == '.fits':
        if img2[-5:] == '.fits':
            hdu_1 = fits.open(img1)
            hdu_2 = fits.open(img2)
        else:
         exit('given none .fits files. Files given: ' + img1 + ' ' + img2)
    else:
        exit('given none .fits files. Files given: ' + img1 + ' ' + img2)

    ######################
    # test if images match
    ######################

    # Do they have the same wcs at (1000,1000)
    w1 = WCS(hdu_1[1].header)                           #get the WCS of the images
    w2 = WCS(hdu_2[1].header)
    position1 = w1.all_pix2world(1000., 1000., 1)       #what is the sky possition of pixel 1000, 1000
    position2 = w2.all_pix2world(1000., 1000., 1)
    if (abs(position1[0] - position2[0] ) > 1e-5) or (abs(position1[1] - position2[1] ) > 1e-5):
        #if either RA or Dec's differeance is too large
        #1e-5 is good. 1e-6 does not show up in 1/1000th of a second?
        badWCS.append(hdu_1[0].header['targname']) #@todo(should we be using a list like this? I think it should be an warning.)

    ######################
    # Combine data
    ######################

    # make final hdu
    #   copy from hdu_1 into final image to transfer all the information
    #   note that WCS is the same for both images. All is corrected (looked at SN1415 as an example)
    #   make a unique second copy. I want it in two places in memory
    hdu_combined = hdu_1

    # Get inverse sensitiveity, the convertion factor
    #inverse seinsitivity [ergs/cm2/Ang/electron]
    invSensitivity_1 = hdu_1[1].header['PHOTFLAM'] 
    invSensitivity_2 = hdu_2[1].header['PHOTFLAM']

    # Combine, with a dance if the sizes are not perfectly the same.
    if shape(hdu_1[1].data) != shape(hdu_2[1].data):
        #if not the same size
        #hard code the miss size issue in SN13038, SN15461 and SN3488. F474W_drz is (2076, 2126) and F625W_drz is (2075, 2126) 
        hdu_combined[1].data = invSensitivity_1*hdu_1[1].data[:2075] +invSensitivity_2*hdu_2[1].data
    else:
        #change the data in hdu_combined to be the addition of hdu_1 and hdu_2
        hdu_combined[1].data = invSensitivity_1*hdu_1[1].data + invSensitivity_2*hdu_2[1].data

    # Correct header
    hdu_combined[0].header['comment'] = ''
    hdu_combined[0].header['comment'] = 'modified at ' + datetime.now().isoformat() 
    hdu_combined[0].header['comment'] = 'This is a combined images from F475W and F625W filters by brose3@nd.edu'
    hdu_combined[0].header['comment'] = 'Headers are from F475W images with an update to "BUNIT".'
    hdu_combined[1].header['BUNIT'] = 'ergs/cm2/Ang/sec'
    ##todo(change maybe a few other parts of the header.)

    # Save images
    target = hdu_1[0].header['targname']
    name = target+'_combined_flux.fits'
    hdu_combined.writeto('data/HST - combined/'+name)
    print('saved '+name)

    # Close files to clear memory
    hdu_1.close()
    hdu_2.close()
    hdu_combined.close()
    return None

def main():
    # Rename files
    # note: drz are final images, flt are pre-cosmic ray subtracted
    drz = glob.glob('data/HST - original/*_drz.fits')
    flt = glob.glob('data/HST - original/*_flt.fits')
    data = drz + flt
    map(rename, data)

    # Combine files
    global badWCS                                       #set up global varriable to store info if there is poor allignment of images
    badWCS = [] 
    F475W = glob.glob('data/HST - renamed/SN*F475W*drz.fits')           
    F625W = glob.glob('data/HST - renamed/SN*F625W*drz.fits')           
    map(coadd, F475W, F625W)       #merge images and save to combineddata folter
    print("bad WCS", badWCS)

    return None

if __name__ == '__main__':
    # main()
    # Combine files
    global badWCS                                       #set up global varriable to store info if there is poor allignment of images
    badWCS = [] 
    F475W = glob.glob('data/HST - renamed/SN*F475W*drz.fits')           
    F625W = glob.glob('data/HST - renamed/SN*F625W*drz.fits')           
    map(coadd, F475W, F625W)       #merge images and save to combineddata folter
    print("bad WCS", badWCS)


    #Outline

    ## Select and rename drz and flt images to be accesable via the SN number

    ## coadd

    ## color