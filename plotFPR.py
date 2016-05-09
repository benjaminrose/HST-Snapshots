""" plotFPR.py -- plot the fractional pixel rank

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 2
    2016-04-05
    Licesed under the MIT License
"""
from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# set plot parameters here, is it for presontation, journal, web?

# read in data
snID = 2635
saveLocation = 'resources/SN{0}/'.format(snID)

galaxyPixels = np.loadtxt(saveLocation+'SN{0}_host_pixel_values.csv'.format(snID), delimiter=',')
snPixelValue = np.loadtxt(saveLocation+'SN{0}_pixel_values.csv'.format(snID), delimiter=',')
fpr = np.loadtxt(saveLocation+'SN{0}_fpr.csv'.format(snID), delimiter=',')

# plot results
## histogram with a gaussian kernel density estimate overtop
# sns.distplot(galaxyPixels)
# plt.plot([snPixelValue, snPixelValue], [0.0, 6.0])

## CDF from a gaussian kernel density estimate overtop
sns.kdeplot(galaxyPixels, cumulative=True) #shade=True
plt.plot([snPixelValue, snPixelValue], [0.0, 1.0])
plt.xlabel('pixel count rate [counts/s]')
plt.ylabel('p(<pixel count rate)')     #the y-axis is the probabilty in a CDF

# save results
# plt.show()
plt.savefig('figures/2016-04-16-fpr-2635.pdf')