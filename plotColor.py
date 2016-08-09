import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# data = pd.read_csv('hst_sdss_color_snr5.csv')
# data = pd.read_csv('hst_sdss_color_snr2.csv')
# data = pd.read_csv('resources/hst_color.csv',header=1, skiprows=[2])
data = np.genfromtxt('resources/hst_color.csv', delimiter=',', names=True, skip_header=1)

plt.figure('color-color correleation')
plt.scatter(data['color'], data['sdss_color'])
plt.xlabel('hst [475W-625W]')
plt.ylabel('sdss [g-r]')

plt.figure('F475W-g correleation')
plt.scatter(data['blue_mag'], data['g_band'])
plt.xlabel('F475W [mag]')
plt.ylabel('g [mag]')

plt.figure('F625W-r correleation')
plt.scatter(data['red_mag'], data['r_band'])
plt.xlabel('F625W [mag]')
plt.ylabel('r [mag]')

plt.show()

# http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html