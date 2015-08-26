'''
# README
* Currently works with the HST images
* Orignal images can be found relative to this file 'data/HST - original/*.fits'
* took 54 mins with 122 objects on 2015-05-29.`re)

## Procedure
### rename(doc)
* doc -> file name of the image you want to rename
* takes fits files (*_drz.fits & *_flt.fits) and remnames them to <<SN>>_<<filter>>_<<img>>_<<count>>.fits
	* <<SN>> is the SDSS tranient number
	* <<filter>> is the HST filter (by number)
	* <<img>> is wether its drz or flt
	* <<count>> optional and only for flt, since there are two images per filter.
* saves renamed as 'data/HST - renamed/*.fits'

### merge(img1, img2)
* file name of the two fits images you want to merge
* comments in new header assume img1 == F475W and img2 == F625W
* combines the two (from each filter) drz images
	* <<filter>> becomes 'combined'
	* <<img>> and <<count>> are dropped
* saves combined as 'data/HST - combined/*.fits'

## Notes
* drz are final images, flt are pre-cosmic ray subtracted
'''

from astropy.io import fits
from astropy.wcs import WCS
from sys import exit
from numpy import shape
import glob

def rename(doc):
	#from file name, is it drz or flt
	img = doc[-8:-5] #assumes file name is *_drz.fits or *_flt.fits

	#open fits file
	hdu = fits.open(doc) 						# open file
	
	# get names
	sn = hdu[0].header['targname'] 				#might have trailing white space
	filter1 = hdu[0].header['filter1'] 			#F475W or F625W
	filter2 = hdu[0].header['filter2'] 			#this should be CLEAR2L
	
	# check that all is ok - could do more
	if filter2 != "CLEAR2L":							#if filter2 is not clear
		sys.exit('Not using clear filter in ' + doc)	# abort with error message

	# save renamed file
	if img == 'flt':				#If they are flt there should be 2 files, needs more work
		attempt = 1 					#Is this image 1 or 2?
		while True:										#make infinate loop
			name = sn+'_'+filter1+'_'+str(attempt)+doc[-9:]		#define name to save images 
			try:
				hdu.writeto('data/HST - renamed/'+name)			#can I save with this attemt?
				break											#if so break
			except IOError:					#if error
				attempt += 1 				#iterate up on attempt (then loop again)
				if attempt > 4: 			#if attempt is too high, exit with error message. NO INFINITE LOOPS!
					sys.exit('Over 3 images of the same SN/filter check '+ doc +' should be ' + name)
	else:							#for drz images
		name = sn+'_'+filter1+doc[-9:]		#define name to save images 
		hdu.writeto('data/HST - renamed/'+name)		#save fits file
	print 'saved '+name	

	# close file to clear memory
	hdu.close()
	return None

def merge(img1, img2):

	# open images	@todo(should I use try?)
	if img1[-5:] == '.fits':
		if img2[-5:] == '.fits':
			hdu_1 = fits.open(img1)
			hdu_2 = fits.open(img2)
		else:
		 exit('given none .fits files. Files given: ' + img1 + ' ' + img2)
	else:
		exit('given none .fits files. Files given: ' + img1 + ' ' + img2)

	# make final hdu
	#	copy from hdu_1 into final image to transfer all the information
	#	note that WCS is the same for both images. All is corrected (looked at SN1415 as an example)
	hdu_combined = hdu_1

	# Do they have the same wcs at (1000,1000)
	w1 = WCS(hdu_1[1].header)							#get the WCS of the images
	w2 = WCS(hdu_2[1].header)
	position1 = w1.all_pix2world(1000., 1000., 1)		#what is the sky possition of pixel 1000, 1000
	position2 = w2.all_pix2world(1000., 1000., 1)
	if (abs(position1[0] - position2[0] ) > 1e-5) or (abs(position1[1] - position2[1] ) > 1e-5):
		#if either RA or Dec's differeance is too large
		#1e-5 is good. 1e-6 does not show up in 1/1000th of a second?
		badWCS.append(hdu_1[0].header['targname']) #@todo(should we be using a list like this?)

	# Are they the same size
	if shape(hdu_1[1].data) != shape(hdu_2[1].data):
		#if not the same size
		#hard code the miss size issue in SN13038, SN15461 and SN3488. F474W_drz is (2076, 2126) and F625W_drz is (2075, 2126) 
		hdu_combined[1].data = hdu_1[1].data[:2075] + hdu_2[1].data
	else:
		#change the data in hdu_combined to be the addition of hdu_1 and hdu_2
		hdu_combined[1].data = hdu_1[1].data + hdu_2[1].data

	# Correct header
	hdu_combined[0].header['comment'] = 'This is a combined images from F475W and FF625W filters by brose3@nd.edu'
	hdu_combined[0].header['comment'] = 'Header is for F475W images'

	# Save images
	target = hdu_1[0].header['targname']
	name = target+'_combined.fits'
	hdu_combined.writeto('data/HST - combined/'+name)
	print 'saved '+name

	# Close files to clear memory
	hdu_1.close()
	hdu_2.close()
	hdu_combined.close()
	return None

def main():
	# Rename files
	drz = glob.glob('data/HST - original/*_drz.fits')
	flt = glob.glob('data/HST - original/*_flt.fits')
	data = drz + flt
	map(rename, data)

	# Combine files
	global badWCS 										#set up global varriable to store info if there is poor allignment of images
	badWCS = [] 
	F475W = glob.glob('data/HST - renamed/SN*F475W*drz.fits') 			
	F625W = glob.glob('data/HST - renamed/SN*F625W*drz.fits') 			
	map(merge, F475W, F625W) 							#merge images and save to combineddata folter
	print "bad WCS", badWCS 

	return None

if __name__ == '__main__':
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
