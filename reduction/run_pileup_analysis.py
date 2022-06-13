import numpy as np
from pileup_grades_ratio import *
from pileup_radial_intensity import *
from pileup_spectral_distortion import *

#====================================================================================================

# Set some parameters that are constant across all observations
outdir  = '/Users/salvesen/research/fcolabs/results/'
pupdir  = '/Users/salvesen/research/fcolabs/data/pileup/'
phacuts = [50, 1000]  # <-- 0.5-10 keV (1 XRT channel = 10 eV)
Nreg    = 26          # <-- Rcut = [0,25] pixels; Rin = [0,25] and Rout = [1,26] pixels

# Annular regions with excised inner circle of radius Rcut and fixed outer radius Rout = 30 pixels
# (Used for the grades ratio and spectral distortion analyses)
Rcut = np.arange(Nreg)  # [pixels]

# Annular regions of 1-pixel width, with inner radius Rin and outer radius Rout
# (Used for the radial intensity analysis)
Rin  = np.arange(Nreg)      # [pixels]
Rout = np.arange(Nreg) + 1  # [pixels]

# Observation IDs
obsID_list = ['00030009021', '00030009022', '00030009023', '00030009025', '00030009026']
Nobs       = len(obsID_list)

print ""
print "PILE-UP ANALYSIS"

#====================================================================================================
'''
# GRADES RATIO

# Loop through each observation ID
for i in np.arange(Nobs):
    
    # Current observation ID
    obsID = obsID_list[i]
    print "...Calculating Grades Ratio for ObsID: ", obsID
    
    # Initialize the source/background file lists
    fsrc0_list  = []
    fbkg0_list  = []
    fsrc02_list = []
    fbkg02_list = []
    # Populate the source/background file lists
    for j in np.arange(Nreg):
        Rstr = 'Rcut' + str(Rcut[j])
        fsrc0_list.append(pupdir  + obsID + '/' + 'source_'     + Rstr + '_grd0_'   + obsID + '.pha')
        fbkg0_list.append(pupdir  + obsID + '/' + 'background_' + Rstr + '_grd0_'   + obsID + '.pha')
        fsrc02_list.append(pupdir + obsID + '/' + 'source_'     + Rstr + '_grd0-2_' + obsID + '.pha')
        fbkg02_list.append(pupdir + obsID + '/' + 'background_' + Rstr + '_grd0-2_' + obsID + '.pha')

    # Calculate the grades ratio for the current observation ID
    fout = outdir + obsID + '/' + 'pileup_grades_ratio_' + obsID + '.h5'
    write_grades_ratio(fout=fout, Rcut=Rcut, phacuts=phacuts, \
        fsrc0_list=fsrc0_list, fbkg0_list=fbkg0_list, \
        fsrc02_list=fsrc02_list, fbkg02_list=fbkg02_list)
'''
#====================================================================================================
'''
# RADIAL INTENSITY PROFILE
# --> PROBLEM: THE CALIBRATED KING RPOFILES ARE NOT FITTING THE DATA!!!
#phacuts = [50, 100]
PSFmin  = 10.0  # [counts / sec]
dE      = None  # [keV]

# Loop through each observation ID
for i in np.arange(Nobs):
    
    # Current observation ID
    obsID = obsID_list[i]
    print "...Calculating Radial Intensity Profile for ObsID: ", obsID

    # Initialize the source/background file lists
    fsrc_list  = []
    fbkg_list  = []
    # Populate the source/background file lists
    for j in np.arange(Nreg):
        Rstr = 'R' + str(Rin[j]) + '-' + str(Rout[j])
        fsrc_list.append(pupdir + obsID + '/' + 'source_'     + Rstr + '_grd0-2_' + obsID + '.pha')
        fbkg_list.append(pupdir + obsID + '/' + 'background_' + Rstr + '_grd0-2_' + obsID + '.pha')

    # Calculate the radial intensity profile for the current observation ID
    fout = outdir + obsID + '/' + 'pileup_radial_intensity_' + obsID + '.h5'
    write_radial_intensity(fout=fout, Rin=Rin, Rout=Rout, phacuts=phacuts, PSFmin=PSFmin, dE=dE, \
        fsrc_list=fsrc_list, fbkg_list=fbkg_list)
'''
#====================================================================================================
# SPECTRAL DISTORTIONS
cntrate0=100.0

# Loop through each observation ID
for i in np.arange(Nobs):
    
    # Current observation ID
    obsID = obsID_list[i]
    print "...Calculating Spectral Distortions for ObsID: ", obsID

    # Initialize the source/background file lists
    fgrp_list = []
    fsrc_list = []
    fbkg_list = []
    # Populate the source/background file lists
    for j in np.arange(Nreg):
        Rstr = 'Rcut' + str(Rcut[j])
        fgrp_list.append(pupdir + obsID + '/' + 'source_grp20_serr3_' + Rstr + '_grd0-2_' + obsID + '.pha')
        fsrc_list.append(pupdir + obsID + '/' + 'source_'     + Rstr + '_grd0-2_' + obsID + '.pha')
        fbkg_list.append(pupdir + obsID + '/' + 'background_' + Rstr + '_grd0-2_' + obsID + '.pha')

    # Calculate the spectral distortion for the current observation ID
    fout = outdir + obsID + '/' + 'pileup_spectral_distortion_' + obsID + '.h5'
    write_spectral_distortion(fout=fout, Rcut=Rcut, cntrate0=cntrate0, phacuts=phacuts, \
        fgrp_list=fgrp_list, fsrc_list=fsrc_list, fbkg_list=fbkg_list)

