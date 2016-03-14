""" fractionalRank.py -- Calculating fractional pixel rank

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-07
    Licesed under the MIT License
"""
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

def get_pixels(SNID, position):
	"""
	Returns the numerical value of the pixels for the host galaxy and the SN

	# Parameters
	SNID : int
		The identification number of the SDSS Supernova, used in file names.

	position : astropy.SkyCoord
		The position of the SN. It should be compatible with the WCS found in 
		the fits file that `SNID` will point to.

	# Returns
	galaxy : np.array
		A 1D array of all the value of the pixels of the galaxy

	sn : float
		The pixel value of the
	"""
	# import HUD
	filePath = 'data/HST - combined/SN{0}_combined.fits'
	hdu, scidata = ancillary.import_fits(filePath.format(SNID), extention=1)

	return get_galaxy_pixels(hdu, scidata), get_sn_value(position, hdu, scidata)

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
		A 1D array of all the value of the pixels of the galaxy
	"""
	# get science data if needed
	if sciData is None:
		sciData = hdu.data

	# get galaxy definition
	# galaxyShape =[95, 1108.1448309503585, 37.95088945611436, 2.43723464012146, 2.371455192565918, -1.4546968936920166]
	

	galaxy = np.array([0,0.1,0.2,1,2,0.2,0.4])
	return galaxy

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
	# numpy.arrays are accessed row then collumn
	# `all_pix2wold() returns (x, y) so we need to call sciData[row #,column #]
	sn = sciData[SNPixels[1,0], SNPixels[0,0]]

	#take median of 3x3 box. 

	sn = 0.3
	return sn


def get_FPR(SNID, position):#, positions, sigma=2, box_size=3):
	"""
	Fractional Pixel Rank (FPR) is CDF value of a particular number, in this 
	case of the pixel that hosted the SN.
	"""
	# get galaxtic pixels as a 1D-array
	galaxy, SN = get_pixels(SNID, position)

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
		# find last entry of where galaxy is less then SN
		rank_min = np.where(galaxy < SN)[0][-1]
		# find first entry of where galaxy is greater then SN (or add 1 to `rank_min`)
		rank_max = rank_min + 1
		fpr_min = 1.0*rank_min/(len(galaxy)-1.0)
		fpr_max = 1.0*rank_max/(len(galaxy)-1.0)
		f = interpolate.interp1d([galaxy[rank_min], galaxy[rank_max]]
								,[fpr_min, fpr_max])

		fpr = f(SN)

	"""
	# import galaxie shape information
	galaxies = Table.read('resources/galaxies_{0}.csv'.format(sigma), format='ascii.commented_header', header_start=1)
	#@todo(fix it so we can change the location of resources)
	#table is: x, y, a, b, theta	
	galaxies['theta'].unit = u.radian

	# for each supernova
	# for i, sn in enumerate(SN_num): #need both SN_num and position to itterate. currently this is a bad method
	for sn, position, x, y, a, b, theta in zip(SN_num, positions, galaxies['x'], galaxies['y'], galaxies['a'], galaxies['b'], galaxies['theta'].quantity):
		print sn, position, x, y, a, b, theta
		# get that host galaxy
		# host = galaxies[galaxies['SN number'] == sn] #still an astropy table
		#'SN number' was first called 'SN name'. A few old files might still have this issue

		# import HST image 
		hdu, data = ancillary.import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn), hdu_return=True)
		
		# init mask that defines pixels belonging to host
		# mask = np.zeros((2075,2126)) #currently what is in SN1415, but I could make this dynamic If I first import fits image
		mask = np.zeros(data.shape)
		#why all the [0], currently host produces arrays and we don't want that?


		# init variables for ellipse equation
		# ctheta = np.cos(host['theta'].quantity.to(u.radian).value[0]) #should already be radians, but lets just be sure.
		ctheta = np.cos(theta.to(u.radian).value) #should already be radians, but lets just be sure.
		stheta = np.sin(theta.to(u.radian).value) #should already be radians, but lets just be sure.
		center = (x, y)
		

		# search a section of mask, and update part inside to be 1.
		r = np.ceil(max(a, b))
		#can't search for all of mask but we can assume its near a circle defined at center & r = max(a,b)
		for x_index in np.arange(-r, r+1)+center[0]: #check indexes, these are not whole numbers!
		    for y_index in np.arange(-r, r+1)+center[1]:
		        x_can = (x_index - center[0])*ctheta + (y_index - center[1])*stheta #defing the canonical part of the equation
		        y_can = -(x_index - center[0])*stheta + (y_index - center[1])*ctheta
		        if (x_can**2)/(a/2)**2 + (y_can**2)/(b/2)**2 <= 1: #maybe do 1.01? or something so that more pixels are countend. It currently can look funny.
		            # mask[center[1]-r+l, center[0]-r+m] = 1 #it is in a (y,x) coordiante system
		            mask[y_index, x_index] = 1.0 #indices floor everything that are given to it.

		# get data that is inside galaxy, & rank it
		ranked_soon = mask*data
		ranked_flat = ranked_soon.flatten()
		ranked =  np.sort(ranked_flat[ranked_flat != 0])

		# get pixel value of SN postion (as a skycoord)
		w = WCS(hdu[1].header) #can't just do WCS('filename') because of HST has multihead fits files

		x_sn, y_sn = w.all_world2pix(position.ra.to(u.deg).value, position.dec.to(u.deg).value, 1) #this should be better, do it all at once?
		#all_world2pix gives non-whole numbers
		x_sn, y_sn = round(x_sn), round(y_sn) #can do np.round (returns lists) or round (returns value)
		# print x, y
		# print data[y,x]
		# print data[x,y]
		sn_value = data[y_sn,x_sn]

		# I want: xy x-1y-1, x-1y, x-1y+1
		sn_values = np.array([])
		 #could make this an input varriable
		for i in range(box_size):
			for j in range(box_size):
				sn_values = np.append(sn_values, data[y_sn-1+i,x_sn-1+j])	
		sn_value = np.mean(sn_values)
		# print sn_values, sn_value, data[y,x]

		# convert SN value to rank
		inside = True
		if mask[y_sn,x_sn] == 0:
			inside = False

		'''
		'''
		print x_sn, y_sn
		#elipse over data
		x_plot = np.arange(center[0]-r, center[0]+r)
		y_plot = np.arange(center[1]-r, center[1]+r)
		plt.figure(2)
		plt.pcolormesh( x_plot, y_plot, data[y_plot[0]:y_plot[-1], x_plot[0]:x_plot[-1]])
		plt.plot(x_sn, y_sn, marker='*', markersize=15, markerfacecolor='r', markeredgecolor='w') #the location of SN1415
		plt.colorbar()
		e = Ellipse(center, a, b, theta.to(u.deg).value)
		e.set_facecolor([1,.733333333,0])
		e.set_alpha(0.25)
		plt.gca().add_patch(e)

		plt.show()
	return ranked, sn_value, inside
	"""
	return fpr

def save_rank(galactic, sn, inside):
	''' what do I want to save?
	* galactic list, sn value, if inside?
	* galactic list - for cdf plot, sn fraction, sn cdf value, 
	I should be able to easily change what to do if they are inside or not
	'''
	return None #line can be deleted.


def main():
	"""
	This is the default method for calculucating fractional pixel rank.
	"""
	SNID = 2635
	position = get_SN_HST_coord(names)

	galaxyShape = [95, 1108.1448309503585, 37.95088945611436, 2.43723464012146, 2.371455192565918, -1.4546968936920166]

	# rank = rank_galactic_pixels([1415], 3)	
	# sn = get_SN_pixel([1415], position) #this is odd. Why am I reimporting the from fits files?

	# galactic, sn, inside  = rank_supernova([1415], position, 1.5)
	galactic, sn, inside  = rank_supernova(names, position, 1.5)

	sigma = 1.5
	sigma_iterate = np.ones(len(names))*sigma
	# stuff = map(rank_supernova, names, position, sigma_iterate)
	# print 	galactic, sn, inside



	print stuff
	from sys import exit
	exit()



	print np.where(sn < galactic)[0][0]

	SN_fractional = sn/galactic[-1] 
	print SN_fractional

	cdf_value = np.where(sn < galactic)[0][0]/float(len(galactic))
	#cdf is numder left of you in the order divided by lenght of data set. 
	print cdf_value

	print np.where(5.0 < galactic)[0][0]/float(len(galactic))
	
	
	return None

if __name__ == "__main__":
	# main()
	SNID = 2635
	position = get_SN_HST_coord(SNID)
	
	# galaxy, SN = get_pixels(SNID, position) #a sup-part of `get_FPR`, for testing
	# print galaxy, SN

	fpr = get_FPR(SNID, position)#, position)
	print fpr





def get_SN_pixel(SN_num, position):

	for i, sn in enumerate(SN_num):
		#import HST image @todo(should be a function call)
		hdu = fits.open('data/HST - combined/SN{0}_combined.fits'.format(sn))
		data = hdu[1].data #the location of science data in HST multi extention FITS images
		data = data.byteswap(True).newbyteorder()
		#do I need to close this?

		#get pixel value of SN postion (as a skycoord)
		hdu = fits.open('data/HST - combined/SN{0}_combined.fits'.format(sn))
		w = WCS(hdu[1].header)

		#can't just do WCS('filename') because of HST has multihead fits files
		x, y = w.all_world2pix(position[i].ra.to(u.deg).value, position[i].dec.to(u.deg).value, 1) #this should be better, do it all at once?
		#all_world2pix gives non-whole numbers
		x, y = round(x), round(y) #can do np.round (returns lists) or round (returns value)
		print x, y

		#should it be data[x,y] or data[y,x]?
		#world2pix looks like (x,y) therefore data[y,x]
		print hdu[1].data[y,x]
		print hdu[1].data[x,y]

		#check if its outside the galaxy?


	return None