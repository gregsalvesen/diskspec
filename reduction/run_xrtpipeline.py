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
#=============================================================================#
# GRO J1655-40
datdir = '/Users/salvesen/research/fcolabs/data/'
#wrtdir = '/Users/salvesen/research/fcolabs/data/clean/'
wrtdir = '/Users/salvesen/research/fcolabs/data/pileup/'
ra     = '"16 54 00.137"'
dec    = '"-39 50 44.90"'
obsIDs = ['00030009021',
          '00030009022',
          '00030009023',
          '00030009025',
          '00030009026']

Nobs = len(obsIDs)
for i in range(Nobs):
    reduce_data(obsID  = obsIDs[i],
                regsrc = None,
                regbkg = None,
                ra     = ra,
                dec    = dec,
                datdir = datdir,
                wrtdir = wrtdir,
                runxrt = True,
                runimg = True,
                runsrc = False,
                runbkg = False,
                runarf = False,
                runrmf = False,
                rungrp = False)

#=============================================================================#
