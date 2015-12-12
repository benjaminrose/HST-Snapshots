'''ancillary functions
'''

import numpy as np
from astropy.io import fits		#to read fits files
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob 					#for geting names from the file system.
import re

def import_fits(image, extention=1, hdu_return=False):
	'''
	Import a fits imagae. 

	# Parameters:
	image: str
		The relative path to the fits image.

	extention: int
		Allows for use of mutli-extention FITS images. HST puts data in extention=1.

	hdu_return: bool
		Defalts to false, and only returns the sciece data for the extention given. 
		But if set to true it returns the full file. 

	# Returns:
	data: np.array
		A 2D array of the science data.

	OR

	hdu: something
		The fits image as an object. Contains header and all extentions
			hdu.header
			hdu.data

	# Examples:

		sn = 1415
		hdu, data = import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn), hdu_return=True)
		data = import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn))
	'''

	hdu = fits.open(image)

	data = hdu[extention].data #the location of science data in HST multi extention FITS images
	data = data.byteswap(True).newbyteorder() 
	#I am unsure why this is needed, but it is!

	if hdu_return:
		return hdu, data

	return data

def get_sn_names(data_location = 'data/HST - combined/'):
	'''
	Get the SDSS transient ID number for all files stored in data_location.
	Default location contains names generated for raw data's fits-header-target name.
	Technically the action is to get characters from postition 2 through -14. 
	
	# Parameters
	data_location: string
		The string to the directory storeing the data. 

	# Returns
	names: np.array of strings formated as 'a6'
		The SDSS transient ID number for these SN. Can be padded with to be length of 6. 
	'''
	files = glob.glob(data_location + '*')
	names = np.zeros(len(files), dtype='a6') #'a6' allows for a 6 character sting. other wise we cant pad in place.

	for i, fil in enumerate(files):
		names[i] = fil[len(data_location)+2:-14]
		#asssume files are saves like 'SN1415_combined.fits' or 'SN15451_combined.fits'

	return names

#@todo(what is the SDSS error? can we avgerae over 9 HST pixels or do we not know the SN's location enough.)
def getSDSSPosition(SN):
	'''
	Imports SN position from SDSS data. Data is found *. Function returns an array of SkyCoords
	# Parameters
	SN: numpy array of strings of 6 characters ('a6')
		list of sdss 

	# Returns
	SN_position: np.array of SkyCoord
		The position of each SN.
	'''
	SN_position = np.zeros(len(SN), dtype=object)
	
	for i, sn in enumerate(SN):
		sn = str(sn).zfill(6) #pad with zeros to match SMP
		#zfill needs it to be a string, some whow np.array can mess with that.

		with open('data/SDSS - photometry/SMP_{0}.dat'.format(sn), 'r') as f:
			first_line = f.readline()
			#split on white space, convert to numpy array so np.where works
		split = np.array( re.split(r'\s+', first_line) )
		ra_val_where = np.where(np.array(split) == 'RA:')[0][0]+1
		dec_val_where = np.where(np.array(split) == 'DEC:')[0][0]+1
		# print ra_val_where, dec_val_where
		# print float(split[ra_val_where])*u.deg, type(float(split[ra_val_where]))
		SN_position[i] = SkyCoord(ra = float(split[ra_val_where])*u.deg, dec = float(split[dec_val_where])*u.deg)

	return SN_position

# I from fits files I, get sci data, read header, and get wcs. 
#Maybe these are simple enough with astropy I should not write my own functions?

if __name__ == "__main__":
	#for testing:
	getSDSSPosition(['8297'])