''' minimizingFspsResultsTesting.py - Figureing out how to go from model to age +/- uncertantiy

Benjamin Rose
brose3@nd.edu/benjamin.rose@me.com
python3
2016-12-19
'''
import re
import warnings
from glob import glob
from copy import deepcopy
from datetime import datetime

import numpy as np
import pandas as pd
from astropy.io import fits

import matplotlib.pyplot as plt
import seaborn as sns

#for model creating
import fsps

#for "Comparing Varriables"
from scipy import integrate
# from astropy.cosmology import WMAP9 as cosmo   # or make my own
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27) ## Kessler 2009a but flat and 2-sig figs (Gupta 2011, §2.3)


def getData():
    #INPUT 1
    #note that data is in second header, `cause Yeah fits!
    fCampbell = fits.open('data/SDSS_Photometric_SNe_Ia.fits')

    #fix `ValueError: Big-endian buffer not supported on little-endian compiler`
    fCampbell[1].data = fCampbell[1].data.byteswap(True).newbyteorder() 
    data = deepcopy(fCampbell[1].data)  

    #make data into a `DataFrame`
    SNIaCam = pd.DataFrame(data)
    #might not want to `set_index` since this is a column we want.
    # SNIaCam.set_index('CID', inplace=True)

    #INPUT 2
    # add columns for [u, u uncert, g, g uncert, r, r uncert, i, i uncert, z, z uncert] 
    SNIaCam['u'] = np.nan
    SNIaCam['u uncert'] = np.nan
    SNIaCam['g'] = np.nan
    SNIaCam['g uncert'] = np.nan
    SNIaCam['r'] = np.nan
    SNIaCam['r uncert'] = np.nan
    SNIaCam['i'] = np.nan
    SNIaCam['i uncert'] = np.nan
    SNIaCam['z'] = np.nan
    SNIaCam['z uncert'] = np.nan

    #INPUT 3
    dataFile = 'data/SDSS - photometry/SMP_{}.dat'    #need to add a s6, error padded integer

    for i, CID in enumerate(SNIaCam['CID']):
        CID = str(int(CID)).zfill(6)
        magitudes = readMagnitudes(dataFile.format(CID))
        SNIaCam.loc[i, ['u', 'u uncert', 'g', 'g uncert', 'r', 'r uncert', 'i', 'i uncert', 'z', 'z uncert']] = magitudes

    return SNIaCam

def readMagnitudes(dataFile):
    """
    This reads and returns the mangitudes from an SDSS SMP data file
    
    # Parameters
    dataFile : str
        The name of the SMP data file. Include any needed file paths.
        
    # Returns
    magnitudes : np.array
        A structured array made some stupid repeating BS. So it got tossed. 
        The aruments come out in [u, u uncert, g, g uncert, r, r uncert, 
        i, i uncert, z, z uncert] order.
    """
    with open(dataFile, 'r') as f:
        f.readline()
        secondLine = f.readline()
        thirdLine = f.readline()
    
    asinhmag = False
    if asinhmag:
        dataLine = secondLine
    else:
        dataLIne = thirdLine
        
    #clean up data line and extract desired info
    split = np.array( re.split(r'\s+', dataLIne) )     #split on whitespace
    #the index of the g (r) values are 3 (4) spaces over from the first '='
    index = np.arange(5)+ 5       #(np.where(split == '=')[0]+3)
    SB = split[index].astype(float)
    
    if any(SB < 0):
        #send warning if we got a negative flux! But still out put the data as is.
        warnings.warn(r'{} has a negative flux: u,g,r,i,z = {} μJy/sqr-arcsec'.format(dataFile, SB))
        #todo(fix warning so that it outputs the correct units depending on any `asinhmag`)
    #         magnitudes = np.nan*np.ones(10)

    uncertIndex = index + 8
    uncert = split[uncertIndex].astype(float)

    if not asinhmag:
        SB = -2.5*np.log10(SB*1e-6/3631)
        # or could be 1.0857*uncertG/gSB
        uncert = uncert*np.abs(2.5/(SB*np.log(10)))

    #combine as [value, uncert, ...], and save as a structured array
    magnitudes = np.array(np.stack((SB, uncert), axis=1).flatten(),) 
    
    return magnitudes

def readModelSpace(fileName='resources/fspsModelSpace.csv'):
    """
    reads saved, via `saveModelSapce()`, model space and reshapes it.
    
    # Parameters
    
    fileName : str
        
    # Returns
    
    modelSpace : multidimentional numpy array
        The *ugriz* outputs of the poly dimentional model space.
    """
    #get data
    data = np.loadtxt(fileName, delimiter=',')
    
    #get final shape
    with open(fileName, 'r') as f:
        firstLine = f.readline()
    
    shapeString = re.search(r'\d,.*\d', firstLine).group(0)
    shape = np.fromstring(shapeString, dtype=int, sep=',')
    
    #create modelSpace object
    modelSpace = data.reshape(shape)
    return modelSpace

if __name__ == '__main__':
    #full data
    SNIaCam = getData()
    #for clean data: `SNIaCam.dropna()`

    fileName = 'resources/fspsModelSpace/fspsModelSpace_0.06539160013198853.csv'
    modelSpace = readModelSpace(fileName)
    print(modelSpace)