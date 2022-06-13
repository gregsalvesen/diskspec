import numpy as np
from reduce_data import *

'''
NEEDS WORK!!!

PURPOSE:
--------
Generate cleaned Swift/XRT WT mode V4641 Sgr data ready for spectral fitting in XSPEC.
'''

# Unlike in pileup_clean(), we do not filter the pha cannels here. We will mark 0-29 as bad in grppha later (per usual) and then ignore below 0.5 keV in XSPEC.

# Raw data directory and clean data directory to write to
datdir = '/Users/salvesen/research/fcolabs/data/'
wrtdir = '/Users/salvesen/research/fcolabs/data/clean/'

# Source (RA, DEC)
ra     = '"18 19 21.634"'
dec    = '"-25 24 25.849"'

# Region files for the source and background
regsrc = datdir + 'regions/' + 'V4641_source_annulus.reg'
regbkg = datdir + 'regions/' + 'V4641_background_annulus.reg'

# Filters: Event grades, grouping, systematic error
grade  = '0-2'
grpmin = 20
syserr = 0.03

# Loop through each observation and reduce the data
obsIDs = ['00030111050', '00030111051', '00030111052', '00030111053', '00030111054', '00030111055', '00030111056', '00030111057']
Nobs   = len(obsIDs)
for i in range(Nobs):
    reduce_data(obsID=obsIDs[i], regsrc=regsrc, regbkg=regbkg, grade=grade, grpmin=grpmin, \
        syserr=syserr, ra=ra, dec=dec, datdir=datdir, wrtdir=wrtdir, runxrt=False)
