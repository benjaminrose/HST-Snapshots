# Find HST and SDSS guide stars
#To detemine the offest of the HST WCS in data and position of SN from SDSS data.
#http://docs.astropy.org/en/stable/vo/conesearch/index.html is uesful
#might want to use `sgsc` or `imgsc` from WCSTools in iraf? I could not get it to work, but it might be better if it worked.

from astropy.vo.client import conesearch
from astropy import units as u
import numpy as np
import sys

def findGuides(SN_location):
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
	sdss_cat = 'SDSS DR8 - Sloan Digital Sky Survey Data Release 8 2'

	guides_hst  = np.array([])
	guides_sdss = np.array([])
	for i in SN_location:
		#need to have a check to see if there is notheing with in 0.5*u.arcmin
		try:
			guides_hst  = np.append(guides_hst,  conesearch.conesearch(i, 0.5*u.arcmin, catalog_db=hst_cat ) )
		except:
			#I should do something, escpeciall if it is an error I expect for nothing found.
			sys.exc_info()
		try:
			guides_sdss = np.append(guides_sdss, conesearch.conesearch(i, 0.5*u.arcmin, catalog_db=sdss_cat) )
		except:
			sys.exit()
	return guides_hst, guides_sdss

def findSameStar():
	return 0

if __name__ == "__main__":
	from astropy.coordinates import SkyCoord
	SN = [SkyCoord(ra=6*u.degree, dec=0.56*u.degree, frame='icrs'), #SN1415
		SkyCoord(ra=8*u.degree, dec=-0.56*u.degree, frame='icrs')	#random & does not work with in 0.5*u.arcmin in sdss.
		]

	hst, sdss = findGuides(SN)
	for i, j in zip(hst, sdss):
		print i.array['ra'].data, i.array['dec'].data
		print j.array['ra'].data, j.array['dec'].data
