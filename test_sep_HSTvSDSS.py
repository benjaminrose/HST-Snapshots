import numpy as np

import defGalaxy
import ancillary

#define fits files I want
Supernova = '1415'
sdssFilter = 'g' #currently not used in file names, but will be.
hstLocation = 'data/HST - combined/SN{0}_combined.fits'.format(Supernova)
sdssLocation = 'data/SDSS - coadded/SN{0}-{1}.fit'.format(Supernova, sdssFilter)

#from looking at fits images: the approxiamt middle of the galaxy desired
hstPixels = (1,1)
sdssPixels = (1,1)

#import sci data from hst 
hst = ancillary.import_fits(hstLocation)
print 'HST sources'
hst_sources = defGalaxy.run_sep(hst, sigma=1.5, get_all=True)
print hst_sources[['x', 'y', 'a', 'b', 'theta']]

#import sci data from sdss 
sdss = ancillary.import_fits(sdssLocation, extention=0)
print '\n\n SDSS Sources'
sdss_sources = defGalaxy.run_sep(sdss, sigma=1.5, get_all=True)
print sdss_sources[['x', 'y', 'a', 'b', 'theta']]

'''for SN1415
 HST - (1075.6214058144144, 1042.2975354911164, 30.01184844970703, 19.923686981201172, 1.3260775804519653)
 SDSS-g - (1700.2162240162525, 829.3687612510705, 7.608850955963135, 4.741743087768555, 0.7577474117279053)
 '''