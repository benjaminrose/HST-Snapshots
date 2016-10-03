""" blockReduce.py -- reduces the size of fits images 

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-09-16
    Licesed under the MIT License
"""

from copy import deepcopy
from datetime import datetime

import numpy as np
from astropy.io import fits
from astropy import wcs
#block_reduce needs scikit-image
from astropy.nddata.utils import block_reduce

import ancillary

def main(SNNumber):
    #########read in orginal hdu
    imageFile = 'data/HST - combined/SN{}_combined_flux.fits'.format(SNNumber)

    ##########import original image data
    hduOrig, dataOrig = ancillary.import_fits(imageFile, extention=1)
    header0 = hduOrig[0].header
    header1 = hduOrig[1].header

    ##############reduce science data
    block_size = 8
    data = block_reduce(dataOrig, block_size)

    ###############set up WCS
    #HST uses CDi_j nomincature (Section 2.1.2). Each pixel step is not 8
    # times larger, but no new roations are needed. Paper describing WCS at
    # http://adsabs.harvard.edu/abs/2002A%26A...395.1061G
    # also http://astropy.readthedocs.io/en/stable/wcs/index.html

    w = wcs.WCS(naxis=2)

    #note CD2_1 and CD1_2 are the same from SN1415, so I guessed where they go.
    w.wcs.cd = [[ header1['CD1_1']*block_size,header1['CD1_2']*block_size ],
                [ header1['CD2_1']*block_size,header1['CD2_2']*block_size ]]
    w.wcs.crval = [header1['CRVAL1'], header1['CRVAL2']]
    w.wcs.crpix = [header1['CRPIX1']/block_size, 
                   header1['CRPIX2']/block_size]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN' ]

    headerWCS = w.to_header()

    ################### check if wcs's are the same
    #get wcs for original data
    wOrig = wcs.WCS(header1)
    #define pixels to test, not perfect corners, but is divisible by 8!
    pixOrig = np.array([[0, 0], [0, 2072], [1060, 1036], [2120, 0], 
                        [2120, 2072]], np.float_)
    pix = pixOrig/block_size
    #get world values
    worldOrig = wOrig.wcs_pix2world(pixOrig, 1)
    world = w.wcs_pix2world(pix, 1)
    #throw error if ra or dec of any is off by more then 1e-7, 1"==2.7e-4˚
    assert np.max(np.abs(worldOrig - world)) < 1e-7

    ################### copy other header and add new a few chagnes
    #update main header
    headerMain = header0
    headerMain['comment'] = ''
    headerMain['comment'] = 'modified at ' + datetime.now().isoformat() + ' by brose3@nd.edu'
    headerMain['comment'] = 'used astorpy.block_reduce to reduce image by {}'.format(block_size)
    headerMain['comment'] = 'scalled up CDi_j by {}'.format(block_size)

    ###################### update science header with bits other then WCS
    # could pull from `header1`
    headerWCS['EXTNDAME'] = 'SCI'
    headerWCS['BUNIT'] = 'ergs/cm2/Ang/sec'
    headerWCS['PHOTMODE'] = 'ACS WFC1 F475W MJD#55355.1483'
    headerWCS['PHOTFLAM'] = 1.8267533E-19
    headerWCS['PHOTZPT'] = -2.1100000E+01
    headerWCS['PHOTPLAM'] = 4.7456221E+03
    headerWCS['PHOTBW'] = 4.2024533E+02

    #################create new HDU's
    hduMain = fits.PrimaryHDU(header=headerMain)
    hduSci = fits.ImageHDU(data=data, header=headerWCS)

    #####################save new hdu
    target = header0['targname']
    name = target+'_reduced8.fits'
    hdulist = fits.HDUList([hduMain, hduSci])
    hdulist.writeto('data/HST - reduced/'+name)
    print('saved '+name)


if __name__ == '__main__':
    # main('1415')

    #get names
    SN = ancillary.get_sn_names()
    SN = np.array(SN, dtype=np.int)
    list(map(main, SN))