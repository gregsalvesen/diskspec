import numpy as np
import glob
from reduce_data import *

'''
PURPOSE:
--------
Clean all of the data for the pile-up analysis.
'''

datdir  = '/Users/salvesen/research/fcolabs/data/'
wrtdir  = '/Users/salvesen/research/fcolabs/data/pileup/'
ra      = '"16 54 00.137"'
dec     = '"-39 50 44.90"'
regsrc  = glob.glob(datdir + 'regions/' + 'source_annulus_Rcut*.reg')    # <-- Rout = 30 pixels
dRsrc   = glob.glob(datdir + 'regions/' + 'source_annulus_R[0-9]*.reg')  # <-- dR = 1 pixel
regbkg  = datdir + 'regions/' + 'background_annulus.reg'
grades  = ['0', '0-2']
phacuts = None
grpmin  = 20
syserr  = 0.03
obsIDs  = ['00030009021', '00030009022', '00030009023', '00030009025', '00030009026']
Nobs    = len(obsIDs)
Nreg    = len(regsrc)
Ngrd    = len(grades)

# Loop through each observation
for i in range(Nobs):

    # Loop through each region
    for j in range(Nreg):

        # Loop through different event grade choices (0 vs. 0-2)
        for k in range(Ngrd):

            # Rcut <-- Rout = 30 pixels
            Rcut = regsrc[j].split('/')[-1].split('.')[0].split('Rcut')[-1]
            tag  = 'Rcut' + Rcut + '_grd' + grades[k] + '_'
            print "Reducing data with tag ", tag
            reduce_data(obsID=obsIDs[i], regsrc=regsrc[j], regbkg=regbkg, ra=ra, dec=dec, \
                tag=tag, grade=grades[k], grpmin=grpmin, syserr=syserr, phacuts=phacuts, \
                datdir=datdir, wrtdir=wrtdir, runxrt=False, runimg=False)

            # R <-- dR = 1 pixel
            Rin  = dRsrc[j].split('/')[-1].split('.')[0].split('R')[-1].split('-')[0]
            Rout = dRsrc[j].split('/')[-1].split('.')[0].split('R')[-1].split('-')[1]
            tag  = 'R' + Rin + '-' + Rout + '_grd' + grades[k] + '_'
            print "Reducing data with tag ", tag
            reduce_data(obsID=obsIDs[i], regsrc=dRsrc[j], regbkg=regbkg, ra=ra, dec=dec, \
                tag=tag, grade=grades[k], grpmin=grpmin, syserr=syserr, phacuts=phacuts, \
                datdir=datdir, wrtdir=wrtdir, runxrt=False, runimg=False)

