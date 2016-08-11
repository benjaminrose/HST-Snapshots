import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn

from scipy.optimize import curve_fit

# data = pd.read_csv('hst_sdss_color_snr5.csv')
# data = pd.read_csv('hst_sdss_color_snr2.csv')
# data = pd.read_csv('resources/hst_color.csv',header=1, skiprows=[2])
data = np.genfromtxt('resources/hst_color.csv', delimiter=',', names=True, skip_header=1)

# fit color data
def func(x, m, b):
    return m*x + b
color = np.linspace(-0.2, 2.0)      #to plot resulting function
badData1 = np.isnan(data['color'])
badData2 = np.isnan(data['sdss_color'])
badData = np.logical_or(badData1, badData2)
hst = data['color'][~badData]
sdss = data['sdss_color'][~badData]
popt, pcov = curve_fit(func, hst, sdss)

print(popt)
print(pcov)
fig = plt.figure('color-color correleation')
ax = fig.add_subplot(111)
ax.scatter(data['color'], data['sdss_color'])
ax.plot(color, func(color, *popt))
plt.xlabel('hst [475W-625W]')
plt.ylabel('sdss [g-r]')
plt.xlim(0.2, 1.5)
ax.text(0.3, 1.5, r"y = m*x + b")
ax.text(0.3, 1.4, r"m = {0:.3f} +/- {1:.3f}".format(popt[0], pcov[0,0]**0.5))
ax.text(0.3, 1.3, r"b = {0:.3f} +/- {1:.3f}".format(popt[1], pcov[1,1]**0.5))
plt.title('Same size "pixels", SNR > 20')

# plt.figure('F475W-g correleation')
# plt.scatter(data['blue_mag'], data['g_band'])
# plt.xlabel('F475W [mag]')
# plt.ylabel('g [mag]')

# plt.figure('F625W-r correleation')
# plt.scatter(data['red_mag'], data['r_band'])
# plt.xlabel('F625W [mag]')
# plt.ylabel('r [mag]')

plt.show()

# http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html