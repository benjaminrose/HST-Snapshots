""" plotFPR.py -- plot the fractional pixel rank

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-04-05
    Licesed under the MIT License
"""
from __future__ import print_function, division

from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#also need statsmodels
from astropy.table import Table

# set plot parameters here, is it for presontation, journal, web?

# read in data
def plot_individual(snID = 2635, show=True):
    """
    show : boolian
        if true shows image. if false saves image.
    """
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
    if show:
        plt.show()
    else:
        plt.savefig('figures/2016-04-16-fpr-2635.pdf')

def collect_data(flag):
    """
    flag : string
        do you wnat 'hst' or 'sdss'

    # Return
    fpr : np.array
        (x+1,2) array of strings. first column is SN names, second column is fpr. first row is a header statings previous sentance. Lenght is number of SN plus 1 (header row!)
    """
    location = glob('resources/SN*/SN*_'+flag+'_fpr.csv')

    # creating a header makes more sense then `if len(fprs)==0:`.
    fprs = np.array([['#sn', 'fpr']])
    for i in location:
        name = i[12:17]

        if name[-1] =='/':
            name = name[:-1]

        # for this to work, fprs needs something in it form the start
        fprs = np.append(fprs, [[name, np.genfromtxt(i)]], axis=0)

    return fprs

def main(flag, show=True):
    """
    flag : string
        do you wnat 'hst' or 'sdss'

    show : boolian
    """
    hst = collect_data('hst')
    hst = hst[1:]
    fpr_hst = np.array(hst[:,1], dtype=np.float)
    fpr_hst.sort()
    cdf_hst = np.array(range(fpr_hst.size))/(fpr_hst.size-1)

    sdss = collect_data('sdss')
    sdss = sdss[1:]
    fpr_sdss = np.array(sdss[:,1], dtype=np.float)
    fpr_sdss.sort()
    cdf_sdss = np.array(range(fpr_sdss.size))/(fpr_sdss.size-1)


    print('hst: ', fpr_hst)
    print('sdss: ', fpr_sdss)

    #plot images

    # #plot histogram 
    # plt.figure('histogram')
    # # sns.kdeplot(fpr_hst, cumulative=True) #shade=True
    # sns.distplot(fpr_hst, bins=10, rug=True)
    # plt.xlabel('Fractional Prixel Rank, of HST galaxies')
    # plt.ylabel('frequency-ish')    #but not really.
    # plt.xlim(-0.5, 1.5)

    # #plot sdss-cdf
    # plt.figure('sdss-cdf')
    # plt.plot(fpr_sdss, cdf_sdss)
    # plt.xlabel('Fractional Prixel Rank, of HST galaxies')
    # plt.ylabel('p(<pixel count rate)')     #the y-axis is the probabilty in a CDF

    # #plot hst-cdf
    # plt.figure('hst-cdf')
    # plt.plot(fpr_hst, cdf_hst)
    # plt.xlabel('Fractional Prixel Rank, of HST galaxies')
    # plt.ylabel('p(<pixel count rate)')     #the y-axis is the probabilty in a CDF

    #plot combined cdf
    plt.figure('combined-cdf')
    plt.plot(fpr_sdss, cdf_sdss, label='sdss')
    plt.plot(fpr_hst, cdf_hst, label='hst')
    plt.legend(loc=0)
    plt.xlabel('Fractional Prixel Rank')
    plt.ylabel('p(<pixel count rate)')     #the y-axis is the probabilty in a CDF

    # save results
    if show:
        plt.show()
    else:
        plt.savefig('figures/test.pdf')

if __name__ == '__main__':
    main('hst')
    # collect_data('hst')