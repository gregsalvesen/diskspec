import numpy as np
from reduce_data import *

'''
PURPOSE:
--------
Produce clean data and an image with no filters

NOTES:
------
- Must initialize HEAsoft
'''

# V4641 Sgr (For now, skip obsID 00030111048, which is in PC mode)
datdir = '/Users/salvesen/research/fcolabs/data/'
wrtdir = '/Users/salvesen/research/fcolabs/data/clean/'
#wrtdir = '/Users/salvesen/research/fcolabs/data/pileup/'
ra     = '"18 19 21.634"'
dec    = '"-25 24 25.849"'
obsIDs = ['00030111050', '00030111051', '00030111052', '00030111053', '00030111054', '00030111055', '00030111056', '00030111057']
Nobs   = len(obsIDs)
for i in range(Nobs):
    reduce_data(obsID=obsIDs[i], regsrc=None, regbkg=None, ra=ra, dec=dec, datdir=datdir, wrtdir=wrtdir, \
        runxrt=True, runimg=True, runsrc=False, runbkg=False, runarf=False, runrmf=False, rungrp=False)
