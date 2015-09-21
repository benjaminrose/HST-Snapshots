'''ancillary functions
'''

import numpy as np
from astropy.io import fits		#to read fits files
import glob 					#for geting names from the file system.

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