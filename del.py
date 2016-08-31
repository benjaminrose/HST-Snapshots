from __future__ import print_function, division

from glob import glob
import numpy as np
#also need statsmodels
from astropy.table import Table

''' # getting fpr's from their file sturcture
def collect_data(flag):
    """
    flag : string
        do you wnat 'hst' or 'sdss'

    # Return
    fpr : np.array
        (x+1,2) array of strings. first column is SN names, second column is fpr. first row is a header statings previous sentance. Lenght is number of SN plus 1 (header row!)
    """
    location = glob('resources/SN*/SN*_'+flag+'_fpr.csv')

    # creating a header makes more sense then `if len(fprs)==0:` in the for-loop.
    fprs = np.array([['#sn', 'fpr']])
    for i in location:
        name = i[12:17]

        if name[-1] =='/':
            name = name[:-1]

        # for this to work, fprs needs something in it form the start
        fprs = np.append(fprs, [[name, np.genfromtxt(i)]], axis=0)

    return fprs

hst = collect_data('hst')
hst = hst[1:]
fpr_hst = np.array(hst[:,1], dtype=np.float)

sdss = collect_data('sdss')
sdss = sdss[1:]
fpr_sdss = np.array(sdss[:,1], dtype=np.float)


print(fpr_hst)
print('\n\n')
print(fpr_sdss)

'''




#getting HST ra & dec
import fractionalRank
import ancillary

SN = ancillary.get_sn_names()
SN = np.array(SN, dtype=np.int)

for sn in SN:
    pos = fractionalRank.get_SN_HST_coord(sn)
    # print(pos.to_string('hmsdms'))