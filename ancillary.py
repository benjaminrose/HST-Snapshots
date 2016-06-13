"""ancillary.py -- a storage file for arbitrary resued fuctions.

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-06-13    #when comments were added
    Licesed under the MIT License
"""
import glob 					#for geting names from the file system.
import re

import numpy as np
from astropy.io import fits		#to read fits files
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

def import_fits(image, extention=0):
	'''
	Imports a fits imagae. Takes the path and returns both the hdu object and the science data, as a numpy array.

	# Parameters:
	image: str
		The relative path to the fits image file.

	extention: int
		This is the extention that you want `data` to be taken from. Default is `0` but a common exception is HST where it should be `1`. 

	# Returns:
	hdu: astropy.fits.hdu
		The fits image as an object, contains header and all extentions.
			hdu.header
			hdu.data

	data: np.array
		A 2D array of the science data.

	# Examples:
		sn = 1415
		hdu, data = import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn))
	'''

	hdu = fits.open(image)

	data = hdu[extention].data
	# I am unsure why this is needed, but it is!
	#Astropy talks about bid and little indian stuff
	data = data.byteswap(True).newbyteorder() 

	return hdu, data

def update_dataset():
	'''
	Not functioning yet.

	Updates `resources/dataset.csv`. It first finds if there are any new SN in 
	`data/HST - combined` then goes out and fills in `SDSS ID`, `RA`, `dec`. All other
	entries need to be entered someother way (`delta RA`, `delta Dec`, & `SDSS DR12 obj ID` 
	were all enetered by hand.)

	Assumes files names are like `SN1415_combined.fits` or `SN15451_combined.fits`
	'''
	#get names
	data_location = 'data/HST - combined/'
	try:
		files = glob.glob(data_location + '*')
	except Exception, e:
		import warnings
		warnings.warn('file import from `{0}` failed with warning: {1}'.format(data_location, e))

	names = np.zeros(len(files), dtype='str') #'a6' works instead of 'str', but this is more flexiable
	for i, fil in enumerate(files):
		names[i] = fil[len(data_location)+2:-14]

	#import csv
	try:
		data = Table.read('resources/dataset.csv', format='ascii.commented_header')
	except Exception, e:
		import warnings
		warnings.warn('importing `resources/dataset.csv` failed with warning: {0}'.format(e))
	#todo(unit support?)

	#test if new things
	print data['SDSS ID'][0], names

	for i in data['SDSS ID']:
		print i == any(names)
		for j in names:
			if str(i) == str(j):
				print 'match'
				break
			print 'no match'

		#get RA & Dec


	#save new `data` to file


	return None

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
	files = glob.glob(data_location + '*' + '_flux' + '*')
	names = np.zeros(len(files), dtype='a6') #'a6' allows for a 6 character sting. other wise we cant pad in place.

	for i, fil in enumerate(files):
		names[i] = fil[len(data_location)+2:-19]
		#asssume files are saves like 'SN1415_combined_flux.fits' or 'SN15451_combined_flux.fits'

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
	# getSDSSPosition(['8297'])
	update_dataset()