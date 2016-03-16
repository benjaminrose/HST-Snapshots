import numpy as np
import matplotlib.pyplot as plt
import sep
from sys import exit

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits     #to read fits files


def main():
	#get SN details
	SN = {
		'number' : '2635',
		'location' : SkyCoord(52.704247*u.deg, -1.238154*u.deg),
		#for contrast & etc of the visulaized data
		'vmin' : -0.03, #was -1 for fits file
		'vmax' : 0.7, #and 0.4; 0.9 in first python draft
		#for physical sky scale
		'x' : 980+np.arange(26), #need 26 to close the last conection.
		'y' : 1000*np.ones(26),
		'length' : '1"'
	}

	#import fits file
	data = 'data/HST - combined/SN{}_combined.fits'.format(SN['number'])
	hdu = fits.open(data)
	extention = 1 #default for HST
	sci_data = hdu[extention].data #the location of science data in HST multi extention FITS images
	sci_data = sci_data.byteswap(True).newbyteorder()

	#smooth image
	stdev = 3
	if stdev > 0:
		gauss = Gaussian2DKernel(stdev)
		smoothed_data = convolve(sci_data, gauss) 

	#bin (you need to combine the data so you get more light per pixel)
	'''
	This will take an input science image and binn it up so stdev becomes 1 new pixel.
	'''
	rebined = np.zeros(sci_data.shape) #keep it the same shape as

	#for each 10x10 block: np.sum(data[0:11,0:11])
	x_newPixel_edges = np.arange(0,sci_data.shape[0],10) #but I think I still need sci_data.shape[0]?
	x_newPixel_edges = np.append(x_newPixel_edges, sci_data.shape[0])

	y_newPixel_edges = np.arange(0,sci_data.shape[1],10) #but I think I still need sci_data.shape[1]?
	y_newPixel_edges = np.append(y_newPixel_edges, sci_data.shape[1])

	for i in np.arange(len(x_newPixel_edges)):
		if x_newPixel_edges[i] == sci_data.shape[0]:
			break #break out if last edge
		for j in np.arange(len(y_newPixel_edges)):
			if y_newPixel_edges[j] == sci_data.shape[1]:
				break #break out if last edge
			rebined[x_newPixel_edges[i]:x_newPixel_edges[i+1], y_newPixel_edges[j]:y_newPixel_edges[j+1]] =  np.sum(sci_data[ x_newPixel_edges[i]:x_newPixel_edges[i+1],y_newPixel_edges[j]:y_newPixel_edges[j+1] ])

	#Test image differnece
	print sci_data.shape, rebined.shape

	plt.figure(1)
	plt.imshow(sci_data, vmin=SN['vmin'], vmax=SN['vmax'], origin='lower', cmap='Greys')

	plt.figure(2)
	plt.imshow(rebined, vmin=SN['vmin'], vmax=10*SN['vmax'], origin='lower', cmap='Greys')

	#remove background
	# data = rebined #for testing between the two
	for data in [sci_data, smoothed_data, rebined]:
		sigma = 1.5
		bkg = sep.Background(data)
		thresh = sigma*bkg.globalrms

		#run sep
		# I have a try statement in defGalaxy_deprecated_20160108.py
		all_sources = sep.extract(data, thresh, minarea=50)
		
		#select center object as galaxy
		center = (2090/2.0, 2108/2.0)
		search_r = 200

		x_obj = all_sources['x'] #the x value for all objects
		y_obj = all_sources['y']
		id_center = []
		for i, x in enumerate(x_obj):
			if ((x-center[0])**2 + (y_obj[i]-center[1])**2) < search_r**2:
				id_center.append(i)

		# defGalaxy_deprecated_20160108.py has a much longer testing list
		idx = np.where( all_sources['npix']==max(all_sources['npix'][id_center]) )

		host_galaxy = all_sources[['x', 'y', 'a', 'b', 'theta']][idx].view(np.float64)
		print host_galaxy

	plt.show()

if __name__ == "__main__":
	main()