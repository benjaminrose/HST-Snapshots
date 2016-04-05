""" fractionalRank.py -- Calculating fractional pixel rank

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-07
    Licesed under the MIT License
"""
from __future__ import print_function, division
import re #for regular expressions!

import numpy as np
from scipy import interpolate
from astropy.coordinates import SkyCoord, Latitude, Longitude
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS


#to be removed when this call gets moved to seperate file
from astropy.io import fits		#to read fits files
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import ancillary

def get_SN_SDSS_coord(SNID):
	"""
	Gets the cooridnates, in SDSS system, of a given SN defined by its SDSS SN 
	identification number. The data is taken from photometric data stored in 
	`data/SDSS - photometry/` from 
	`http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html`.
	#todo(figure out the paper that this is from)

	# Parameters 
	SNID: int
		The identification number of the SDSS Supernova, used in file names.
	
	# Returns
	SNPosition : astopy.SkyCoord
		The position of the SN found in SDSS Photomentric data. 
	"""
	data_location='data/SDSS - photometry/'

	# read from SDSS files SN's RA & Dec
	SNID_string = str(SNID).zfill(6)        #pad with zeros to match SMP
	with open(data_location+'SMP_{0}.dat'.format(SNID_string), 'r') as f:
		first_line = f.readline()
	
	# split on white space, convert numpy array so np.where works
	split = np.array( re.split(r'\s+', first_line) )
	RAValueLocation = np.where(np.array(split) == 'RA:')[0][0]+1
	DecValueLocation = np.where(np.array(split) == 'DEC:')[0][0]+1

	# create SkyCoord
	RA = float(split[RAValueLocation])*u.deg
	Dec = float(split[DecValueLocation])*u.deg
	
	SNPosition = SkyCoord(ra = RA, dec = Dec)
	#todo(Determing if this is the correct frame: FK5, icrs, etc.)

	return SNPosition

def get_SN_HST_coord(SNID):
	"""
	Gets the cooridnates, in HST system, of a given SN defined by its SDSS SN 
	identification number. The data has a shift applied to change from the 
	SDSS to HST WCS.

	# Parameters
	SNID : int
		The identification number of the SDSS Supernova, used in file names.

	# Returns
	SNPosition : astopy.SkyCoord
		The position of the SN, with SDSS-to-HST shift applied. Shifts are 
		defined in `resources/shift.csv` and reasoning can be found in 
		`README.md`
	"""
	# get position from SDSS information
	SDSS_SNPosition = get_SN_SDSS_coord(SNID)

	#skip 13038 because I don't have a shift
	if SNID == 13038:
		import warnings
		warnings.warn('SN13038 does not have a shift')
		SNPosition = SDSS_SNPosition
	else:
		# import shift data
		#todo(change this location to be dataset.csv)
		shift = Table.read('resources/shift.csv', format='ascii.commented_header')
		#todo(add units to table, from shift.meta)
		shift['delta RA'].unit, shift['delta Dec'].unit = u.arcsec, u.arcsec

		# apply shift with delta = HST - SDSS or HST = delta + SDSS
		deltaRA = Longitude(
			shift[shift['SDSS SN Name'] == SNID]['delta RA'].quantity
			)
		deltaDec = Latitude(
			shift[shift['SDSS SN Name'] == SNID]['delta Dec'].quantity
			)
		SNPosition = SkyCoord(ra = deltaRA + SDSS_SNPosition.ra,
							  dec = deltaDec + SDSS_SNPosition.dec)

	return SNPosition

def get_galaxy_pixels(hdu, sciData=None):
	"""
	Returns the numerical value of the pixels for the host galaxy of `SNID`

	# Parameters
	hdu : astropy.io.fits.hdu.hdulist.HDUList
		The fits object from the data extension. This needs to cantain a 
		WCS. Assumes its an HST style object with WCS in extention `1`.

	sciData : np.array, optional
		A 2D array of the sciecne data from the CCD Chip. If science data is 
		not given then it will be extracted. Provide pre-extracted data if any 
		manipulation is needed prior to analaysis, ie byteswap.

	# Returns
	galaxy : np.array
		A 1D array of all the value of the pixels of the galaxy. Unsorted, 
		or not?
	"""
	# get science data if needed
	if sciData is None:
		sciData = hdu.data

	# init mask that defines pixels belonging to host
	mask = np.zeros(sciData.shape)

	# get galaxy definition & format its astorpy.table object
	galaxyData = Table.read('resources/2635-galaxy.csv', format='ascii.commented_header', header_start=1)
	# table headers are: SN ID, npix, x, y, a, b, theta
	# add apropriate units to the table.
	galaxyData['theta'].unit = u.radian  
	# make the rows follow the SNID
	# SN ID might be 'SN number' or 'SN name'. I have issues with conistency
	galaxyData.add_index('SN ID')
	#todo(I might need to not read this in each time, but rather take it as a parmeter?)
	
	hostParams = galaxyData.loc[2635]
	#todo(2635 needs to be a parameter)
	hostParams = Table(hostParams)     #convert to a Table because Rows suck

	# sn, position, x, y, a, b, theta in zip(SN_num, positions, galaxies['x'], galaxies['y'], galaxies['a'], galaxies['b'], galaxies['theta'].quantity):

	# init variables for ellipse equation
	ctheta = np.cos(hostParams['theta'].to(u.radian).value)
	stheta = np.sin(hostParams['theta'].to(u.radian).value)

	# search a section of mask, and update part inside to be 1.

	# define a search radius as the largest of the two axis + 5 pixels
	r = np.ceil(max(hostParams['a'].quantity, hostParams['b'].quantity))

	for x_index in np.arange(-r, r+1)+hostParams['x'].quantity:
		for y_index in np.arange(-r, r+1)+hostParams['y'].quantity:
			#defing the canonical part of the equation
			x_can = (x_index - hostParams['x'].quantity)*ctheta + (y_index - hostParams['y'].quantity)*stheta
			y_can = -(x_index - hostParams['x'].quantity)*stheta + (y_index - hostParams['y'].quantity)*ctheta
			if (x_can**2)/(hostParams['a'].quantity/2)**2 + (y_can**2)/(hostParams['b'].quantity/2)**2 <= 1: 
				mask[int(y_index), int(x_index)] = 1.0 
	            #todo(does this work correcty. make a test to plot the resulting mask and ellipse. Should x_index and y_index be ints?)

	# create a nd.array of the pixel values inside the galaxy.
	sciDataFlattened = (mask*sciData).flatten()
	ranked = np.sort(sciDataFlattened[sciDataFlattened != 0])

	return ranked

def get_sn_value(position, hdu, sciData=None):
	"""
	Returns the numerical value of the pixels for the SN

	# Parameters
	position : astropy.SkyCoord
		The position of the SN. It will be used with the WCS found in `hdu`.

	hdu : astropy.io.fits.hdu.hdulist.HDUList
		The fits object from the data extension. This needs to cantain a 
		WCS. Assumes its an HST style object with WCS in extention `1`.

	sciData : np.array, optional
		A 2D array of the sciecne data from the CCD Chip. If science data is 
		not given then it will be extracted. Provide pre-extracted data if any 
		manipulation is needed prior to analaysis, ie byteswap.

	# Returns
	sn : float
		The pixel value of the
	"""
	# get science data if needed
	if sciData is None:
		sciData = hdu.data

	# get world cordinate system
	# can't just do `WCS('filename')` because of HST has multihead fits files
	w = WCS(hdu[1].header)

	# Get SN's pixle position
	# astopy.wcs.WCS only degrees as floats
	# origin is `0` because we will be working with `sciData`
	SNPixels = w.all_world2pix(
		position.ra.to(u.deg).value, 
		position.dec.to(u.deg).value, 
		0
		)
	# can't use `SNPixels.round()` becuase `w.all_world2pix` returns a list!
	SNPixels = np.round(SNPixels)        # now a veritcal np.array

	# Get value of SN's pixel
	#take median of 3x3 box. 
	sn = np.array([])
	for i in [-1,0,1]:
		for j in [-1,0,1]:
			# numpy.arrays are accessed row then collumn and `all_pix2wold()`
			# returns (x, y) so we need to call sciData[row #,column #]
			sn = np.append(sn, sciData[SNPixels[1,0]+i, SNPixels[0,0]+j])
	sn = np.median(sn)
	#todo(do I want median or mean?)

	return sn


def get_FPR(galaxy, SN):#, positions, sigma=2, box_size=3):
	"""
	Fractional Pixel Rank (FPR) is CDF value of a particular number, in this 
	case of the pixel that hosted the SN.
	"""
	#todo(if SN pixel value is less then galaxy edge cut off (look this up), then break and return 0 or -1? This minght need to be done in `get_pixels()` or `get_pixels()` might not be a useful function.)

	# calculate the FPR with range [0,1] both inclusive.
	galaxy.sort()
	
	# You can't simple do the interpelation for the whole CDF. CDFs are too 
	# jumpy/verticle for `scipy.interpolate` to work well. It works fine 
	# between two points but not for the whole range. See notes on 2016-03-14.

	# since galaxy is a 1d array, there is a non-needed tuple wrapper around the result of `np.where()`
	rank = np.where(galaxy == SN)[0]

	#if `SN` is found in `galaxy` calculate its FPR (or the mean of its FPRs) driectly
	if len(rank) > 0:
		# subtact 1 so that the denominator is equal to the largest rank can be.
		# Notes are available from 2015-03-09
		fpr = 1.0*rank/(len(galaxy)-1.0)
		if len(fpr) > 1:
			fpr = fpr.mean()
	# if `SN` is NOT found in `galaxy` calculate its FPR by interpolation
	else:
		# find nearist neighbors in galaxy and coresponding CDF values
		# find last entry position of where galaxy is less then N
		# `np.where()` returns a tuple where we only care about arg-0.
		#breaks if SN value is smaller then galaxy. #todo(fix)
		rank_min = np.where(galaxy[galaxy < SN])[0][-1]
		# find first entry of where galaxy is greater then SN (or add 1 to `rank_min`)
		rank_max = rank_min + 1
		fpr_min = 1.0*rank_min/(len(galaxy)-1.0)
		fpr_max = 1.0*rank_max/(len(galaxy)-1.0)
		f = interpolate.interp1d([galaxy[rank_min], galaxy[rank_max]]
								,[fpr_min, fpr_max])

		fpr = f(SN)
	return fpr

def save_rank(galactic, sn, inside):
	''' what do I want to save?
	* galactic list, sn value, if inside?
	* galactic list - for cdf plot, sn fraction, sn cdf value, 
	I should be able to easily change what to do if they are inside or not
	'''
	return None #line can be deleted.


def main(SNID = 2635):
	"""
	This is the default method for calculucating fractional pixel rank.
	"""
	# get SN position
	print('getting SN position')
	position = get_SN_HST_coord(SNID)

	# get hdu and extract science data
	print('getting HDU and SciData')
	filePath = 'data/HST - combined/SN{0}_combined.fits'
	hdu, scidata = ancillary.import_fits(filePath.format(SNID), extention=1)

	# get 
	print('getting galaxy pixels')
	galaxy_pixels = get_galaxy_pixels(hdu, scidata)
	print(galaxy_pixels)
	print('getting SN pixels')
	sn_pixel_value = get_sn_value(position, hdu, scidata)
	print(sn_pixel_value)

	# get Fractional Pixel Rank
	print('getting FPR')
	fpr = get_FPR(galaxy_pixels, sn_pixel_value)
	print(fpr)

	# save data
	save_location = 'resources/SN{0}/'.format(SNID)
	#todo(make `save_location` directory)
	#todo(what in the world am I doing? How do I want this saved? Updating a csv? Creating a new one? -- a new one for at least galaxy.)
	np.savetxt(save_location+'SN{0}_host_pixel_values.csv'.format(SNID), galaxy_pixels, delimiter=',', header="The pixel values of SN{} host galaxy. The galaxy's edge is defined hard coded currently".format(SNID))
	np.savetxt(save_location+'SN{0}_pixel_values.csv'.format(SNID), [sn_pixel_value], delimiter=',', header='The pixel value of SN'+str(SNID))
	np.savetxt(save_location+'SN{0}_fpr.csv'.format(SNID), [fpr], delimiter=',', header='The fractional pixel rank calculated for SN'+str(SNID))

	return None

if __name__ == "__main__":
	main()



	# SNID = 2635
	# position = get_SN_HST_coord(SNID)

	# # galaxy, SN = get_pixels(SNID, position) #a sup-part of `get_FPR`, for testing
	# # print galaxy, SN

	# fpr = get_FPR(SNID, position)#, position)
	# print fpr