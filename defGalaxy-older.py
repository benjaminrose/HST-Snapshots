'''README
This file gets SN names and fits files. Runs `sep` to extract host galaxy shape
parameters.  Parameters are then saved at `resources/galaxies_{0}.csv` where `{0}`
is the sigma value used to define the edge of the galaxy.
'''
import numpy as np 				#cause numpy is better
import sep						#for source extraction
import ancillary 				#contains functions used in multiple scripts
from sys import exit

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits

def run_sep(data, sigma = 3, name = None, get_all = False, isSdssData = False):
	'''
	Takes the science data and a sigma and runs sep to search for the host.
	If the option `get_all` is set to `True`, that is all it does an it returns
	the full output of sep. If that option is set to `False` (default) then this
	function returns the paramteres for the most likely canidate for the SN's host.
	Host is largest with in a set radius. If none is found, the radius is
	incrementally increased. If more then one are the same size, the closest to the
	center wins.

	# Parameters
	data: 2D np.array
		The science data. Usually the output of `astropy.io.fits`

	sigma: int, float, ect.
		The amount about above background global rms you want to set a threshold
		for a detection pixel. This is part of determing `thresh` in `sep`, check
		its documentation.

	name: string
		The name or other identifer for the current SN being worked on. Used mosly when
		debugging when the fuction is calleded via `map()`.

	get_all: bool
		A flag to return all found objects or just the elliptical parameters of the most
		likely host (default).

	isSdssData: bool
		A flag that determins if this is run on sdss data. If so then the images are cropped 
		around the expected SN location sep is ran on an area of the sky similar to the HST images.
		This allows for the same aglorithum to be used to determine the correct host object is identified.

	# Returns
	host_galaxy: np.array
		An array of the elliptical parrameters of the likely host galaxy. Parameters are
		`x`, `y` (the center), `a`, `b` (the axes lenths, NOT radii!), `theta` (the rotation
		of the ellipse in radians).

	all_sources: np.array
		A structured array that is the standard output of `sep`. Nothing about it is changed.
	'''
	###############
	#convert SDSS Data to look like HST data
	if isSdssData:
		#get SN location
		#data is like data/SDSS - coadded/SN8297-i.fits or SN13354-g.fits
		try: 
			int(data[-12]) #if item -19 is an int, then ID is 5 numbers long.
			SDSS_id = data[-12:-7].zfill(6)
		except ValueError:
			SDSS_id = data[-11:-7].zfill(6)
		
		SN_location = ancillary.getSDSSPosition([SDSS_id])[0]
		#function needs it to be a iterable, and returns an array. Do it this way for only one locatoin.
		hdu = fits.open(data)

		extention = 0 #default for SDSS
		try:
			hdu[0].header.rename_keyword('RADECSYS', 'RADESYS')
		except ValueError:
			pass
		w = WCS(hdu[0].header) #can also use io.fits.getheader('filename')
		x_sn, y_sn = w.all_world2pix(SN_location.ra.to(u.deg).value, SN_location.dec.to(u.deg).value, 1)

		#crop data to look like HST
		#todo(check this out again.)
		img_width = 2000.0/4 #HST pixels * ratio of SDSS/HST pixel size
		img_hight = 2000.0/4

		sci_data = hdu[extention].data #the location of science data in HST multi extention FITS images
		sci_data = sci_data.byteswap(True).newbyteorder() 

		#but what if I hit the edge of the image?
		sci_data = sci_data[x_sn-img_width/2:x_sn+img_width/2, y_sn-img_hight/2:y_sn+img_hight/2]

		##save as new fits files?
	else:
		#data is like 'data/HST - combined/SN*_combined.fits' with * being 4 or 5 numbers
		# I write these lines all the time!!
		extention = 1 #default for HST
		hdu = fits.open(data)
		sci_data = hdu[extention].data #the location of science data in HST multi extention FITS images
		sci_data = sci_data.byteswap(True).newbyteorder() 



	###############
	#get backgournd and threshold
	bkg = sep.Background(sci_data)
	thresh = sigma*bkg.globalrms

	#get sources
	try:
		all_sources = sep.extract(sci_data, thresh, minarea=50)
		#default has minarea=5, SN1415 has 211 (I think with sigam=1)
		#we get too many sources in some images if we have sigma=1
	except Exception, e: # (`internal pixel buffer full` is the expected error)
		exit('error found in {0}: {1}'.format(name, e))

	if get_all:
		return all_sources

	#############
	#finding host:
	#largest obect of the sources s-far away from the middle?

	#find objects near middle
	center = (2090/2.0, 2108/2.0)
	search_r = 200

	x_obj = all_sources['x'] #the x value for all objects
	y_obj = all_sources['y']
	id_center = []
	for i, x in enumerate(x_obj):
		if ((x-center[0])**2 + (y_obj[i]-center[1])**2) < search_r**2:
			id_center.append(i)

	#if nothing exists, incrementally make the radius larger.
	while len(id_center) == 0:
		print "{0} did not find anything".format(name)
		search_r += 100
		for i, x in enumerate(x_obj):
			if ((x-center[0])**2 + (y_obj[i]-center[1])**2) < search_r**2:
				id_center.append(i)

	#get largest of middle objects.
	#where is it that 'npix' is the same as max of 'npix' of center objects
	idx = np.where( all_sources['npix']==max(all_sources['npix'][id_center]) )

	#if there are more then one, choose the one closset to the center.
	if len(idx[0]) > 1:
		new_id = ((x_obj[idx]-center[0])**2 + (y_obj[idx]-center[1])**2).argmin()
		idx = (np.array([ idx[0][new_id] ]), )
		#this is the stange from that idx was originally in
		print "{0} had more then one possible source".format(name)


	#collect only the elliptical parameters of the host galaxy!
	#convert it to a useable data type (via the `view` function)
	host_galaxy = all_sources[['x', 'y', 'a', 'b', 'theta']][idx].view(np.float64)
	#if I give it an array it returns an np.ndarray
	#if I give it an int it returns an np.void
	#do I just want a np.array or list?
	#What is all_sources? A structured array. I called for collumns by name, row by index,
		#and a sertain view that gives me a clean numpy at the end.

	return host_galaxy

def save_galaxies(SN, data, sigma):
	'''
	Files are saved in `resources/galaxies_{0}.csv` with `{0}` being the sigma. Note that simga
	may not be just an `int`. You might get a file named `resources/galaxies_1.5.csv`. The file is
	saved with a three line header, first column being `SN`, and the remaining columns being `data`.

	# Parameters
	SN: str
		The names or other identifers for the SN's. This gets prepended to each row of the `data`
		array before being saved.

	data: 2D np.array
		The science data. Usually the output of `astropy.io.fits`

	sigma: int, float, ect.
		The amount about above background global rms you want to set a threshold
		for a detection pixel. This is part of determing `thresh` in `sep`, check
		its documentation.

	# Returns
	None
	'''
	#add names
	#need to do this odd mess so that we can convert SN from string to int.
	#otherwise do something like https://gist.github.com/benjaminrose/3027d82d8c2bb021b562.
	#if `data[0]` is a np structured array, use `data[0].view(np.float64)` instead of just data[0]
	data_tosave = np.append(int(SN[0]), data[0])
	for i, j in zip(SN[1:], data[1:]):
		combined = np.append( int(i), j )
		data_tosave = np.vstack( (data_tosave, combined) )

	#save file
	filename = 'resources/galaxies_{0}.csv'.format(sigma)
	header = '''Definitions of the elliptical shape of galaxies found from using sep & sigma={0}
SN number, x, y, a, b, theta
int, pixel, pixel, pixel, pixel, radians'''.format(sigma)
	np.savetxt(filename, data_tosave, delimiter=',', header=header)#, fmt='%s')

	return None

def main(sigma = 1.5):
	print 'running defGalaxy.py'

	#get files names
	names = ancillary.get_sn_names()
	toimport = np.array([])
	for i in names:
		toimport = np.append(toimport, 'data/HST - combined/SN{}_combined.fits'.format(i))

	#import images into memeory
	data = map(ancillary.import_fits, toimport)
	#Opening all of these at once is ok. The scince data alone will not mess with RAM too much.

	# run sep
	sigma_iterate = np.ones(len(data))*sigma
	sep_results = np.array(map(run_sep, data, sigma_iterate, toimport))
	#need to convert to np.array to perform sliceing, not sure if it is still needed in final version
	#sep_results[0] -> the result for host-0

	#same info from sep
	save_galaxies(names, sep_results, sigma)

if __name__ == "__main__":
	main()


'''
Benjamin Rose
benjamin.rose@me.com
2015-05-20
Copyright (c) 2015 Benjamin Rose

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
