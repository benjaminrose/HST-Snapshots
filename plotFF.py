""" plotFF.py -- plot the fractional flux

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-09-04
    Licesed under the MIT License
"""
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#also need statsmodels
from astropy.table import Table

# set plot parameters here, is it for presontation, journal, web?

def collect_data(flag, cutoff, blockSize=None):
    """
    flag : string
        do you wnat 'hst' or 'sdss'

    cutoff : float
        This matches the cutoff value for the galaxy definition. Likely 26 for 
        `'hst'`, and either 1, 2, or 3 for `'sdss'`.

    blockSize : int
        The integer block size. Used to downsample a data array by applying a 
        function to local blocks.

    # Return
    ff : np.array
        (x+1,2) array of strings. first column is SN names, second column is fractional flux. First row is a header statings previous sentance. Lenght is number of SN plus 1 (header row!)
    """
    if blockSize is None:
        location = glob('resources/SN*/SN*_{}_{}_fractionalFlux.csv'.format(flag, cutoff))
    else:
        location = glob('resources/SN*/SN*_{}_{}_{}_fractionalFlux.csv'.format(flag, cutoff, blockSize))

    # creating a header makes more sense then `if len(fprs)==0:` in the for-loop.
    ffs = np.array([['#sn', 'fractional flux']])
    for i in location:
        name = i[12:17]

        if name[-1] =='/':
            name = name[:-1]

        # for this to work, fprs needs something in it form the start
        ffs = np.append(ffs, [[name, np.genfromtxt(i)]], axis=0)

    return ffs

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
    ff_hst = np.array(hst[:,1], dtype=np.float)
    ff_hst.sort()
    if not zeros:
        ff_hst = ff_hst[ff_hst.nonzero()]
    cdf_hst = np.array(range(ff_hst.size))/(ff_hst.size-1)

    hstReduced = collect_data('hst', 26, 8)
    #read documentation of `collect_data()` to explain slicing.
    hstReduced = hstReduced[1:]
    ff_hstReduced = np.array(hstReduced[:,1], dtype=np.float)
    ff_hstReduced.sort()
    if not zeros:
        ff_hstReduced = ff_hstReduced[ff_hstReduced.nonzero()]
    cdf_hstReduced = np.array(range(ff_hstReduced.size))/(ff_hstReduced.size-1)

    #collect sdss 1-sigma data to plot
    sdss1 = collect_data('sdss', 1)
    #read documentation of `collect_data()` to explain slicing.
    sdss1 = sdss1[1:]
    ff_sdss1 = np.array(sdss1[:,1], dtype=np.float)
    ff_sdss1.sort()
    if not zeros:
        ff_sdss1 = ff_sdss1[ff_sdss1.nonzero()]
    cdf_sdss1 = np.array(range(ff_sdss1.size))/(ff_sdss1.size-1)

    #collect sdss 2-sigma data to plot
    sdss2 = collect_data('sdss', 2)
    #read documentation of `collect_data()` to explain slicing.
    sdss2 = sdss2[1:]
    ff_sdss2 = np.array(sdss2[:,1], dtype=np.float)
    ff_sdss2.sort()
    if not zeros:
        ff_sdss2 = ff_sdss2[ff_sdss2.nonzero()]
    cdf_sdss2 = np.array(range(ff_sdss2.size))/(ff_sdss2.size-1)

    #collect sdss 3-sigma data to plot
    sdss3 = collect_data('sdss', 3)
    #read documentation of `collect_data()` to explain slicing.
    sdss3 = sdss3[1:]
    ff_sdss3 = np.array(sdss3[:,1], dtype=np.float)
    ff_sdss3.sort()
    if not zeros:
        ff_sdss3 = ff_sdss3[ff_sdss3.nonzero()]
    cdf_sdss3 = np.array(range(ff_sdss3.size))/(ff_sdss3.size-1)

    # make referance
    ff_referance = np.linspace(0, 1)

    #plot images

    #plot combined cdf
    plt.figure('combined-cdf')
    plt.plot(ff_sdss1, cdf_sdss1, label=r'sdss, 1$\sigma$')
    plt.plot(ff_sdss2, cdf_sdss2, label=r'sdss, 2$\sigma$')
    plt.plot(ff_sdss3, cdf_sdss3, label=r'sdss, 3$\sigma$')
    plt.plot(ff_hst, cdf_hst, label='hst')
    plt.plot(ff_hstReduced, cdf_hstReduced, label='hst reduced')
    plt.plot(ff_referance, ff_referance, label=r'$\propto$ luminosity')
    plt.legend(loc=0)
    plt.xlabel('Fractional Flux')
    plt.ylabel('Cumulative Distribution')      #the y-axis is the probabilty in a CDF

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