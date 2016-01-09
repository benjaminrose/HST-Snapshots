import numpy as np
import matplotlib.pyplot as plt

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits		#to read fits files

scale_up = 100.0

#define data 
#basic
supernova_number = '14984'
supernova_location = SkyCoord(313.833496*u.deg, -0.092829*u.deg)

#store in a class?
# class supernova(object):
# 	"""docstring for supernova"""
# 	def __init__(self, num):
# 		super(supernova, self).__init__()
# 		self.num = num

#store in a dictionary per SN?		
SN14984 = {
	'number' : '14984',
	'location' : SkyCoord(313.833496*u.deg, -0.092829*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.001, #was -0.36 for fits file
	'vmax' : 0.4, #and 0.68
	#for physical sky scale
	'x' : 990+np.arange(26), #was 980
	'y' : 1020*np.ones(26), #was 1000
	'length' : '1"'
}

SN2635 = {
	'number' : '2635',
	'location' : SkyCoord(52.704247*u.deg, -1.238154*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.03, #was -1 for fits file
	'vmax' : 0.7, #and 0.4; 0.9 in first python draft
	#for physical sky scale
	'x' : 980+np.arange(26), #need 26 to close the last conection.
	'y' : 1000*np.ones(26),
	'length' : '1"'
} #note frame needs to be 1/4 of 257 per side.
SN6057 = {
	'number' : '6057',
	'location' : SkyCoord(52.553547*u.deg, -0.974675*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.02, #was -0.82 for fits file
	'vmax' : 0.8 , #and 1.5 (-0.6 and 1.6)
	#for physical sky scale
	'x' : 975+np.arange(26),
	'y' : 1000*np.ones(26),
	'length' : '1"'
}
SN14279 = {
	'number' : '14279',
	'location' : SkyCoord(18.488691*u.deg, 0.371556*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.04, #was 0.005 for fits file -- but I would like to do it as a log not lin
	'vmax' : 1.0 , #and 15
	#for physical sky scale
	'x' : 950+np.arange(26),
	'y' : 1025*np.ones(26),
	'length' : '1"'
}
SN17886 = {
	'number' : '17886',
	'location' : SkyCoord(54.006268*u.deg, 1.103300*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.04, #was 0.05 for fits file -- but I would like to do it as a log not lin
	'vmax' : 0.8 , #and 1.871279
	#for physical sky scale
	'x' : 800+np.arange(51),
	'y' : 800*np.ones(51),
	'length' : '2"'
} #note frame needs to be 257x257 not half tha per side.
SN12874 = {
	'number' : '12874',
	'location' : SkyCoord(353.964478*u.deg, -0.177243*u.deg),
	#for contrast & etc of the visulaized data
	'vmin' : -0.01, #was -0.01 for fits file
	'vmax' : 0.3 , #and 0.3
	#for physical sky scale
	'x' : 875+np.arange(26),
	'y' : 975*np.ones(26),
	'length' : '1"'
}

SN = SN14984 #to be a for loop later. 

#open data
data = 'data/HST - combined/SN{}_combined.fits'.format(SN['number'])
hdu = fits.open(data)
extention = 1 #default for HST
sci_data = hdu[extention].data #the location of science data in HST multi extention FITS images
sci_data = sci_data.byteswap(True).newbyteorder()

#smooth data
stdev = 1
gauss = Gaussian2DKernel(stdev)
smoothed_data = convolve(sci_data, gauss) 
#should look into boundry issues, but not now 'cause my galaxy is in the center!

#get supernova pixel location
w = WCS(hdu[extention].header)
supernova_pixel = w.all_world2pix(SN['location'].ra.to(u.deg).value, SN['location'].dec.to(u.deg).value, 1)

'''
#save dataa
save_location = '/Users/brose/Desktop/SN{}_convolved.fits'.format(SN14984['number'])
tosave = hdu #copy over all WCS and Headers
tosave[extention].data = smoothed_data
tosave.writeto(save_location)
'''
#limit frame size
frame = 257/4 #sdss images are currently 515x515 from screen capture of skyserver
plt.xlim([supernova_pixel[0]-frame,supernova_pixel[0]+frame])
plt.ylim([supernova_pixel[1]-frame,supernova_pixel[1]+frame])

#display data
sqrt_data = np.sign(smoothed_data)*np.sqrt(np.absolute(smoothed_data))
# plt.imshow(sqrt_data, vmin=-0.4, vmax=0.68, origin='lower', cmap='Greys') #makes the sky noisy
plt.imshow(smoothed_data, vmin=SN['vmin'], vmax=SN['vmax'], origin='lower', cmap='Greys')
#or 0.35
# plt.colorbar()

#display SN
plt.scatter(supernova_pixel[0], supernova_pixel[1], s=60, c='r', marker='*')

#display length bar, 0.04 arcsec/pixel, 25 pixels = 1 arcsec

x = SN['x'] #need 26 to close the last conection. 
y = SN['y']
plt.plot(x,y,'k')
plt.text(x[0], y[0]+5, SN['length'])

# x = supernova_pixel[0]-frame+50+np.arange(151)
# y = (supernova_pixel[0]+frame-25)*np.ones(151)
# plt.plot(x,y,'k')
# plt.text(x[0], y[0]+5, '5"')

plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.show()