''' look at all of the tests!
'''

import ancillary
import fractionalRank

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS

import time
from os import path, makedirs

def correct_galaxy(time, sigma=1.5):
	galaxies = Table.read('resources/galaxies_{0}.csv'.format(sigma), format='ascii.commented_header', header_start=1)
	#table is: x, y, a, b, theta	
	galaxies['theta'].unit = u.radian

	SN_num = ancillary.get_sn_names()
	positions = fractionalRank.get_SN_position(SN_num)

	for sn, position, x, y, a, b, theta in zip(SN_num, positions, galaxies['x'], galaxies['y'], galaxies['a'], galaxies['b'], galaxies['theta'].quantity):
		print sn, position, x, y, a, b, theta

		hdu, data = ancillary.import_fits('data/HST - combined/SN{0}_combined.fits'.format(sn), hdu_return=True)
		w = WCS(hdu[1].header)
		x_sn, y_sn = w.all_world2pix(position.ra.to(u.deg).value, position.dec.to(u.deg).value, 1) #this should be better, do it all at once?
		#all_world2pix gives non-whole numbers
		x_sn, y_sn = round(x_sn), round(y_sn) #can do np.round (returns lists) or round (returns value)

		center = (x, y)
		r = np.ceil(max(a, b))

		print x_sn, y_sn
		
		#elipse over data
		x_plot = np.arange(center[0]-r, center[0]+r)
		y_plot = np.arange(center[1]-r, center[1]+r)
		plt.figure(str(sn))
		plt.pcolormesh( x_plot, y_plot, data[y_plot[0]:y_plot[-1], x_plot[0]:x_plot[-1]])
		plt.plot(x_sn, y_sn, marker='*', markersize=15, markerfacecolor='r', markeredgecolor='w') #the location of SN1415
		plt.colorbar()
		e = Ellipse(center, a, b, theta.to(u.deg).value)
		e.set_facecolor([1,.733333333,0])
		e.set_alpha(0.25)
		plt.gca().add_patch(e)

		folder = 'test_results/correct_galaxy/'+time
		if not path.exists(folder): makedirs(folder)
		plt.savefig(folder+'/SN'+sn) #currently makes png's

		plt.close(str(sn))

if __name__ == "__main__":
	date = time.strftime("%Y-%m-%d-%H%M%S")
	correct_galaxy(date)