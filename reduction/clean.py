import numpy as np
from reduce_data import *

#=============================================================================#
'''
PURPOSE:
--------
Generate cleaned Swift/XRT WT mode GRO J1655-40 data ready for spectral fitting in XSPEC.
'''

#=============================================================================#
# Unlike in pileup_clean(), we do not filter the pha cannels here.
# We will mark 0-29 as bad in grppha later (per usual) and then ignore below 0.5 keV in XSPEC.

# Raw data directory and clean data directory to write to
datdir = '/Users/salvesen/research/fcolabs/data/'
wrtdir = '/Users/salvesen/research/fcolabs/data/clean/'

# Source (RA, DEC)
ra  = '"16 54 00.137"'
dec = '"-39 50 44.90"'

# Region files for the source and background
regsrc = datdir + 'regions/' + 'source_annulus.reg'
regbkg = datdir + 'regions/' + 'background_annulus.reg'

# Filters: Event grades, grouping, systematic error
grade  = '0-2'
grpmin = 20
syserr = 0.03

# Loop through each observation and reduce the data
obsIDs = ['00030009021', '00030009022', '00030009023', '00030009025', '00030009026']
Nobs   = len(obsIDs)
for i in range(Nobs):
    reduce_data(obsID  = obsIDs[i],
                regsrc = regsrc,
                regbkg = regbkg,
                grade  = grade,
                grpmin = grpmin,
                syserr = syserr,
                ra     = ra,
                dec    = dec,
                datdir = datdir,
                wrtdir = wrtdir,
                runxrt = False)

#=============================================================================#
