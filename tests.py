''' tests.py -- look at all of the tests!

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-03-03
    Licesed under the MIT License
'''

from __future__ import print_function, division
from datetime import datetime
from os import path, makedirs

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS

import ancillary
import fractionalRank

def correct_galaxy(n=25):
	"""
	This saves visual images to see if defGalaxy is working.

	# Parameters
	n : int
		The scaling of (a,b) because Source Extractor is stupid.
	"""
	currentTime = datetime.now()

	# Get galaxy information -- from csv file
	galaxies = Table.read('resources/hosts_25.csv', format='ascii.csv', header_start=2)
	#todo(make the table work better, it is getting the first column name as '# SN' rather then 'SN')
	galaxies['theta'].unit = u.radian

	# Get Positions
	positions = map(fractionalRank.get_SN_HST_coord, galaxies['# SN'].data)

	# for each SN
	for sn, position, x, y, a, b, theta in zip(galaxies['# SN'], 
		positions, galaxies['x'], galaxies['y'], galaxies['a'], 
		galaxies['b'], galaxies['theta'].quantity):

		# get scince data
		hdu, data = ancillary.import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn), extention=1)

		# get WCS data & convert SN position to pixel
		w = WCS(hdu[1].header)
		x_sn, y_sn = w.all_world2pix(position.ra.to(u.deg).value, position.dec.to(u.deg).value, 1) #this should be better, do it all at once?
		#all_world2pix gives non-whole numbers
		x_sn, y_sn = round(x_sn), round(y_sn) #can do np.round (returns lists) or round (returns value)

		###### Plot ######
		# figure out region to plot
		center = (x, y)
		r = np.ceil(max(n*a, n*b))
		r += 10     #add a bit of a padding
		x_plot = (center[0]-r, center[0]+r)    #x-range to plot over
		y_plot = (center[1]-r, center[1]+r)

		# set up figure
		fig = plt.figure(str(sn))
		ax = fig.add_subplot(111, projection=w)

		# plot data
		plt.imshow(data, origin='lower', cmap='cubehelix', clim=(0, 1))
		plt.colorbar()
		plt.plot(x_sn, y_sn, marker='*', markersize=15, markerfacecolor='r', markeredgecolor='w')

		# plot elipse
		e = Ellipse(center, n*a, n*b, theta.to(u.deg).value)
		e.set_facecolor([1,.733333333,0])
		e.set_alpha(0.25)
		plt.gca().add_patch(e)

		# clip figure and add labels
		plt.xlim(x_plot[0], x_plot[1])
		plt.ylim(y_plot[0], y_plot[1]) 
		plt.xlabel('RA')
		plt.ylabel('Dec')

		#save figure
		folder = 'test_results/correct_galaxy/'+currentTime.strftime("%Y-%m-%d %H:%M:%S")
		if not path.exists(folder): makedirs(folder)
		plt.savefig(folder+'/SN'+str(sn)+'_scaled'+str(n)+'.pdf')

		#close to same memory
		plt.close(str(sn))

if __name__ == "__main__":
	correct_galaxy()