'''
# README
* I need to know what exact pixel is the SN in the HST image
* I have the SN position (from SDSS) and and HST WCS, do these match?
	* test by looking at guide stars

## Procedure
### findGuideStars(SN_location)
* takes an astropy SkyCoord and looks around it 
	* currently hardcoded.
* returns two lists for the HST and SDSS guide stars near there.

### ?
* Determins if two guide stars are the same. 



## Notes
* 
'''

from astropy.vo.client import conesearch
from astropy.coordinates import SkyCoord
#@todo(check all SkyCoords for the correct frame. What is the default equinox?)
from astropy import units as u
import numpy as np
import sys
import warnings
import re #for regular expressions!
import glob

#@todo(what are the errors on these guide stars, is the shift less than this?)
def findGuideStars(SN_location):
	'''look into Astropy.coordinates.[match_to_catalog_sky, match-coordinate_sky]. These
	two functions might be what I am doing. Search for astropy compare catalogs.
	''' 

	#double check these are the correct names with conesearch.list_catalogs()
	hst_cat = 'Guide Star Catalog v2 1' #thought it might be 'Guide Star Catalog 2.3 1'
	# hst_cat = 'The HST Guide Star Catalog, Version 1.1'
	#check to see if both cycles used the same guide stars. 
	'''
	*flt.fits say (in SCI header) that ORIGIN = 'HSTIO/CFITSIO March 2010'. This 
	looks like it might be the GSC 1.1 (http://tdc-www.harvard.edu/catalogs/hstgsc.html).
	Astopy should support this (in 1.0.1 docs), but I don't see it? It could be 
	GSC 2.3 thought, cause that looks like it was updated near March 2010.
	'''
	# sdss_cat = 'SDSS DR8 - Sloan Digital Sky Survey Data Release 8 2'
	sdss_cat = 'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 1' #Peter thinks an older DR might be better in Stripe 82. And it does! Down from ~20 fails to 3!

	#Search for something within 0.5 arcmin
	#some how np.append does not work, throws up an 'TypeError: invalid type promotion', so do it another way
	guides_hst  = np.zeros(len(SN_location), dtype=object) 
	guides_sdss = np.zeros(len(SN_location), dtype=object) 
	for i, SN in enumerate(SN_location):
		search_size = 162.0*u.arcsec/2.0 #[WFC FOV / 2](http://www.stsci.edu/hst/wfc3/design/at_a_glance/)
		#still seems way to big. think about this again?
		warnings.simplefilter('ignore') #making vo.client.conesearch quieter, can also mute in function (https://github.com/astropy/astropy-tutorials/blob/master/tutorials/vo/conesearch_tutorial.ipynb). 
		
		#find HST guide stars
		try:
			hst = conesearch.conesearch(SN, search_size, catalog_db=hst_cat)
		except:
			#I should do something, escpeciall if it is an error I expect for nothing found.
			sys.exit('hst')
		
		#find SGSS guide stars
		try:
			sdss = conesearch.conesearch(SN, search_size, catalog_db=sdss_cat)
		except:
			# raise
			# sys.exit('sdss')
			print 'adding to the search area for SDSS, for SN {0} at {1}'.format(i, SN)
			continue
			# search_size = search_size + 15*u.arcmin
			# sdss = conesearch.conesearch(SN, search_size, catalog_db=sdss_cat)
		
		warnings.simplefilter('default') #turn warning noise back on. 

		#add to_table() because conesearch gives a 'astropy.io.votable.tree.Table' not a 'astrop.table.Table'
		guides_hst[i] = hst.to_table()
		guides_sdss[i] = sdss.to_table()	
		


	#returns an array of astropy.tables, a table of guide stars near each SN, per catalog
	# guides_hst[0].show_in_browser(jsviewer = True)
	return guides_hst, guides_sdss

def findRADecShift(catalog_1, catalog_2, accuracy = 1*u.arcsec):
	'''
	Running through two list of astropy tables to find matching objects by position (maybe more later).
	Returns the RA & dec shift between these two objects. 
	Each iteration in the catalogs are done independataly. 

	# Parameters 
	catalog_1: an astropy table or a numpy array of astropy tables 
		containing one column of 'ra' and 'dec'. Needs to be the same size as catalog_2
	catalog_2: an astropy table or a numpy array of astropy tables 
		containing one column of 'ra' and 'dec'
	accuracy: astropy angle unit
		What accuracy should be considered the same object. 

	# Returns
	shift: np.array of tupoles of astropy.coordinates.angles.Angle
		each tuple contains the RA & dec shifts needed to convert from catolog_1 to catalog_2.
		therefore `catalog_1 + shift = catalog_2`
	
	sameobject: np.array of tuples of astropy.coordinates.angles.Angle
		stores both table rows of the two objects that are the same in a tuple for each SN
	'''
	#catalogs need to be the same length
	if len(catalog_1) != len(catalog_2):
		sys.exit('catalogs are not the same length')
	cat_size = len(catalog_1)

	#predefine intemeidate variables
	'''cat_1/2_postions: np.array
		is the positions found in the catalog_1/2, with the same sturcture as the input catalogs. '''
	cat_1_postions = np.zeros(cat_size, dtype=object)
	cat_2_postions = np.zeros(cat_size, dtype=object)

	sameobject = np.zeros(cat_size, dtype=object)
	'''shift: output'''
	shift = np.zeros(cat_size, dtype=object)

	#make an array of SkyCoord arrays of all the catalog_1 stars
	#@todo(can use http://astropy.readthedocs.org/en/latest/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.guess_from_table)
	for i, SN1 in enumerate(catalog_1):
		if type(SN1) != int:
			cat_1_postions[i] = SkyCoord( ra=np.array(SN1['ra'])*SN1['ra'].unit, dec=np.array(SN1['dec'])*SN1['dec'].unit, frame='icrs' )
		else:
			print 'SN interation {0} did not get information from catalog 1'.format(i)

	#make a SkyCoord of all the catalog_2 stars
	for i, SN2 in enumerate(catalog_2):
		if type(SN2) != int:
			cat_2_postions[i] = SkyCoord( ra=np.array(SN2['ra'])*SN2['ra'].unit, dec=np.array(SN2['dec'])*SN2['dec'].unit, frame='icrs' )
		else:
			print 'SN interation {0} did not get information from catalog 2'.format(i)

	
	#for each position in catalog_1, calculate ang_dist with all objects in catalog_2.
	#also check to see if separation is < accuracy, if so update sameobject array.

	#for each SN field, for each object (from catalog_1) in that field, check its seperation from
		#each object in catalog_2. Check if that seperation < accuracy, if so it is the same object.
		#save both table information for the object into `sameobject` and calculate seperation and
		#and save that to `shift`. 
		#@todo(test if we don't find match) @todo(what is there are multiple same objects)

	for i, objs_near_SN in enumerate(cat_1_postions):
		if type(objs_near_SN) == int:
			continue
		for j, position in enumerate(objs_near_SN):
			#calculate separation
			sep = position.separation(cat_2_postions[i])

			#check if seperation < accuracy
			if any(sep.to('deg').value < accuracy.to('deg').value):
				#store same object data from tables
				#somehow to get a table, need catalog[numer][array]. j is an int, but same is an array
				same = np.where(sep.to('deg').value<accuracy.to('deg').value)[0]
				sameobject[i] = (catalog_1[i][[j]], catalog_2[i][same])
				
				#save delta ra & dec shift
				#@todo(are these the same referance frame (icrs, J2000/fk5))
				d_ra = cat_2_postions[i][same].ra - cat_1_postions[i][j].ra
				d_dec = cat_2_postions[i][same].dec - cat_1_postions[i][j].dec
				shift[i] = (d_ra, d_dec)
				#if shift is .5 arcsec or 0.0001 degrees it 1/10,000th the size of the moon.

	return shift, sameobject

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

def getSNNames(data_location = 'data/HST - combined/'):
	'''
	getting the SDSS transient ID number for all files stored in data_location.
	Default location contains names generated for raw data's fits-header-target name.
	# Parameters
	data_location: string
		the string to the directory storeing the data. 

	# Returns
	names: np.array of strings formated as 'a6'
		The SDSS transient ID number for these SN.
	'''
	files = glob.glob(data_location + '*')
	names = np.zeros(len(files), dtype='a6') #'a6' allows for a 6 character sting. other wise we cant pad in place.

	for i, fil in enumerate(files):
		names[i] = fil[len(data_location)+2:-14]
		#asssume files are saves like 'SN1415_combined.fits' or 'SN15451_combined.fits'

	return names

def saveShift(name, shift):
	'''Saves the shifts (and SN name) to csv file for future referance and use.  
	Also allows for this code to get 90% and to fill in the rest by hand. 
	'''
	for i in shift:
		if type(i) == int:
			print i, 'int'
		else:
			# print i
			print i[0].to(u.arcsec), i[1].to(u.arcsec)

	return None

def main():
	names = getSNNames()
	hst, sdss = findGuideStars( getSDSSPosition(names) )
	print hst == sdss #why are the catalogs failing on the same objects in `findRADecShift()`, but not in `findGuideStars()`?
	s, ob = findRADecShift(hst, sdss) #returns shift and collection of the same object data
	#object 4 has two RA's & two dec's, how!!! and fix!!!
	for i in s:
		if type(i) == int:
			print i, 'int'
		else:
			# print i
			print i[0].to(u.arcsec), i[1].to(u.arcsec)

	# print ob

	'''
	print getSDSSPosition(np.array(['1415', '2102'], dtype='a6')) #'a6' allows for a 6 character sting. other wise we cant pad in place.
	
	#for testing
	SN = [SkyCoord(ra=6*u.degree, dec=0.56*u.degree, frame='icrs'), #SN1415
		SkyCoord(ra=8*u.degree, dec=-0.56*u.degree, frame='icrs')	#random & does not work with in 0.5*u.arcmin in sdss.
		]

	#findGuideStars
	hst, sdss = findGuideStars(SN)
	# for i, j in zip(hst, sdss):
	# 	print i.array['ra'].data, i.array['dec'].data
	# 	print j.array['ra'].data, j.array['dec'].data
	
	#findRADecShift
	print findRADecShift(hst, sdss)
	'''
	return None

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
furnished to do so, subject trdo the following conditions:

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
