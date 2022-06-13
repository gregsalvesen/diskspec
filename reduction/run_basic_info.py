import numpy as np
from basic_info import *

'''
PURPOSE:
--------

'''

datadir = '/Users/salvesen/research/fcolabs/data/pileup/'
obsIDs  = ['00030009021', '00030009022', '00030009023', '00030009025', '00030009026']
phacuts = [50, 1000]
Nobs    = len(obsIDs)

# Piled-up data
for i in np.arange(Nobs):
    obsID = obsIDs[i]
    fsrc  = datadir + obsID + '/source_Rcut0_grd0-2_' + obsID + '.pha'
    fbkg  = datadir + obsID + '/background_Rcut0_grd0-2_' + obsID + '.pha'
    getinfo(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts)

# Cleaned data
