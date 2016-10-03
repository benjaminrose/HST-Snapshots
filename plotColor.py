""" plotColor.py -- creats plots

    Benjamin Rose
    benjamin.rose@me.com
    Universtiy of Notre Dame
    Python 3
    2016-08-12
    Licesed under the MIT License
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn

from scipy.optimize import curve_fit


def plotCDF(save=False):
    """
    """
    data = pd.read_csv('resources/color_g-r.csv', header=1, skiprows=[2], index_col=0, skipinitialspace=True)

    color_hst = data['hst color'].dropna()
    color_hst.sort_values(inplace=True)
    cdf_hst = np.array(range(color_hst.size))/(color_hst.size-1)
    color_sdss = data['sdss color'].dropna()
    color_sdss.sort_values(inplace=True)
    cdf_sdss = np.array(range(color_sdss.size))/(color_sdss.size-1)
    color_best = data['color'].dropna()
    color_best.sort_values(inplace=True)
    cdf_best = np.array(range(color_best.size))/(color_best.size-1)

    # print(len(cdf_hst), len(color_hst))
    print(color_hst)

    fig = plt.figure('color cdf')
    ax = fig.add_subplot(111)
    # ax.axis('equal')
    ax.plot(color_hst, cdf_hst, label='hst')
    ax.plot(color_sdss, cdf_sdss, label='sdss')
    ax.plot(color_best, cdf_best, label='best')
    plt.ylabel('Cumulative Probability')
    plt.xlabel('color [g-r, or equivelent]')
    plt.legend(loc=0)
    # plt.xlim(0.1, 1.6)
    # plt.ylim(0.1, 1.6)
    # plt.title('Same size "pixels", SNR > 20')

    if save:
        # plt.savefig('figures/2016-08-31-color.pdf')
        plt.savefig('figures/temp.pdf')
    else:
        plt.show()


def plotComparison(save=False):
    """
    """
    data = pd.read_csv('resources/color_color_cord.csv', header=1, skiprows=[2], index_col=0, skipinitialspace=True)
    data.dropna()
    print(type(data['sdss color'].values))

    def chi_sqr_linear(params, x, y, yUncert):
        parvals = params.valuesdict()
        m = parvals['slope']
        b = parvals['intercept']

        model = m*x + b
        return (model-y)/yUncert

    fig = plt.figure('color cdf')
    ax = fig.add_subplot(111)
    # ax.axis('equal')
    # ax.scatter(data['sdss color'], data['hst color'])
    # ax.errorbar(data['sdss color'], data['hst color'],
                # xerr=data['sdss color uncert'], yerr=data['hst color uncert'])
    ax.errorbar(data['hst color'], data['sdss color'],
                xerr=data['hst color uncert'].values,
                yerr=data['sdss color uncert'].values, marker='*')
    plt.xlabel('SDSS g-r [mag]')
    plt.ylabel('HST F475W-F625W [mag]')
    # plt.legend(loc=0)
    # plt.xlim(0.1, 1.6)
    # plt.ylim(0.1, 1.6)
    plt.title('Same size "pixels", SNR > 20')

    if save:
        # plt.savefig('figures/2016-08-31-color.pdf')
        plt.savefig('figures/temp.pdf')
    else:
        plt.show()



if __name__ == '__main__':
    # plotComparison()
    plotCDF()


'''
data20 = np.genfromtxt('resources/color_snr20.csv', delimiter=',', names=True, skip_header=1)
data = data20

# fit color data
def func(x, m, b):
    return m*x + b
color = np.linspace(-0.2, 2.0)      #to plot resulting function
badData1 = np.isnan(data['color'])
badData2 = np.isnan(data['sdss_color'])
badData = np.logical_or(badData1, badData2)
hst = data['color'][~badData]
sdss = data['sdss_color'][~badData]
popt, pcov = curve_fit(func, sdss, hst)

print(popt)
print(pcov)
fig = plt.figure('color-color correleation')
ax = fig.add_subplot(111)
# ax.axis('equal')
ax.scatter(data['sdss_color'],data['color'])
ax.plot(color, func(color, *popt))
plt.ylabel('hst [475W-625W]')
plt.xlabel('sdss [g-r]')
plt.xlim(0.1, 1.6)
plt.ylim(0.1, 1.6)
ax.text(0.3, 1.5, r"y = m*x + b")
ax.text(0.3, 1.4, r"m = {0:.3f} +/- {1:.3f}".format(popt[0], pcov[0,0]**0.5))
ax.text(0.3, 1.3, r"b = {0:.3f} +/- {1:.3f}".format(popt[1], pcov[1,1]**0.5))
plt.title('Same size "pixels", SNR > 20')
'''
# plt.figure('F475W-g correleation')
# plt.scatter(data['blue_mag'], data['g_band'])
# plt.xlabel('F475W [mag]')
# plt.ylabel('g [mag]')

# plt.figure('F625W-r correleation')
# plt.scatter(data['red_mag'], data['r_band'])
# plt.xlabel('F625W [mag]')
# plt.ylabel('r [mag]')

# http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html
