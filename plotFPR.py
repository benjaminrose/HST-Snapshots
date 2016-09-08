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

def collect_data(flag, cutoff):
    """
    flag : string
        do you wnat 'hst' or 'sdss'

    cutoff : float
        This matches the cutoff value for the galaxy definition. Likely 26 for 
        `'hst'`, and either 1, 2, or 3 for `'sdss'`.

    # Return
    fpr : np.array
        (x+1,2) array of strings. first column is SN names, second column is fpr. first row is a header statings previous sentance. Lenght is number of SN plus 1 (header row!)
    """
    location = glob('resources/SN*/SN*_{}_{}_fpr.csv'.format(flag, cutoff))

    # creating a header makes more sense then `if len(fprs)==0:` in the for-loop.
    fprs = np.array([['#sn', 'fpr']])
    for i in location:
        name = i[12:17]

        if name[-1] =='/':
            name = name[:-1]

        # for this to work, fprs needs something in it form the start
        fprs = np.append(fprs, [[name, np.genfromtxt(i)]], axis=0)

    return fprs

def main(show=True, zeros=False):
    """
    show : boolian
        determins if you should show the figure or save it.

    zeros : boolian
        determins if you should plot zero-fpr's or not. default is to NOT show 
        zeros on plot.
    """
    #collect hst data to plot
    hst = collect_data('hst', 26)
    #read documentation of `collect_data()` to explain slicing.
    hst = hst[1:]
    fpr_hst = np.array(hst[:,1], dtype=np.float)
    fpr_hst.sort()
    if not zeros:
        fpr_hst = fpr_hst[fpr_hst.nonzero()]
    cdf_hst = np.array(range(fpr_hst.size))/(fpr_hst.size-1)

    #collect sdss 1-sigma data to plot
    sdss1 = collect_data('sdss', 1)
    #read documentation of `collect_data()` to explain slicing.
    sdss1 = sdss1[1:]
    fpr_sdss1 = np.array(sdss1[:,1], dtype=np.float)
    fpr_sdss1.sort()
    if not zeros:
        fpr_sdss1 = fpr_sdss1[fpr_sdss1.nonzero()]
    cdf_sdss1 = np.array(range(fpr_sdss1.size))/(fpr_sdss1.size-1)

    #collect sdss 2-sigma data to plot
    sdss2 = collect_data('sdss', 2)
    #read documentation of `collect_data()` to explain slicing.
    sdss2 = sdss2[1:]
    fpr_sdss2 = np.array(sdss2[:,1], dtype=np.float)
    fpr_sdss2.sort()
    if not zeros:
        fpr_sdss2 = fpr_sdss2[fpr_sdss2.nonzero()]
    cdf_sdss2 = np.array(range(fpr_sdss2.size))/(fpr_sdss2.size-1)

    #collect sdss 3-sigma data to plot
    sdss3 = collect_data('sdss', 3)
    #read documentation of `collect_data()` to explain slicing.
    sdss3 = sdss3[1:]
    fpr_sdss3 = np.array(sdss3[:,1], dtype=np.float)
    fpr_sdss3.sort()
    if not zeros:
        fpr_sdss3 = fpr_sdss3[fpr_sdss3.nonzero()]
    cdf_sdss3 = np.array(range(fpr_sdss3.size))/(fpr_sdss3.size-1)


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

    # make referance
    fpr_referacne = np.linspace(0, 1)
    cdf_referance = fpr_referacne**2

    #plot combined cdf
    plt.figure('combined-cdf')
    plt.plot(fpr_sdss1, cdf_sdss1, label=r'sdss, 1$\sigma$')
    plt.plot(fpr_sdss2, cdf_sdss2, label=r'sdss, 2$\sigma$')
    plt.plot(fpr_sdss3, cdf_sdss3, label=r'sdss, 3$\sigma$')
    plt.plot(fpr_hst, cdf_hst, label='hst')
    plt.plot(fpr_referacne, cdf_referance, label=r'$\propto$ luminosity')
    plt.legend(loc=0)
    plt.xlabel('Fractional Pixel Rank')
    plt.ylabel('Cumulative Distribution')     #the y-axis is the probabilty in a CDF

    # save results
    if show:
        plt.show()
    else:
        plt.savefig('figures/test.pdf')

if __name__ == '__main__':
    # main('hst', zeros=True)
    main()
    # main('hst', zeros=False, show=False)
    # collect_data('hst')