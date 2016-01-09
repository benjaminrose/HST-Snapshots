import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits		#to read fits files

SN_number = '14984'
SN_location = SkyCoord(313.833496*u.deg, -0.092829*u.deg)

# import three sdss images.
#make sdss color image from g-r-i or u-g-r, look at Sloan algorithm.
filters = ['i', 'r', 'g'] #order should be RGB?
data = ['data/SDSS - coadded/SN{}-{}.fits'.format(SN_number, filters[0]),
	'data/SDSS - coadded/SN{}-{}.fits'.format(SN_number, filters[1]),
	'data/SDSS - coadded/SN{}-{}.fits'.format(SN_number, filters[2])]

extention = 0 
hdu_0 = fits.open(data[0])
data_0 = hdu_0[extention].data
data_0 = data_0.byteswap(True).newbyteorder()

hdu_1 = fits.open(data[1])
data_1 = hdu_1[extention].data
data_1 = data_1.byteswap(True).newbyteorder()

hdu_2 = fits.open(data[2])
data_2 = hdu_2[extention].data
data_2 = data_2.byteswap(True).newbyteorder()

#are they centered the same?

#renominalize from [0.0,1.0]
def norm(a):
	return (a-a.min())/(a.max()-a.min())
data_0 = norm(data_0)
data_1 = norm(data_1)
data_2 = norm(data_2)

#combine science data 
d = np.stack([data_0,data_1,data_2], axis=2)
# d = np.stack([data_0,2*data_0,data_0], axis=2)

#plot
vmin = -13.0/150000
vmax = 140.0/150000
plt.imshow(d, origin='lower', vmin=vmin, vmax=vmax)

plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.show()