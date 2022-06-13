import numpy as np
from Swift_XRT_WTmode import *
from xspec_models import *

'''
PURPOSE:
--------
For a Swift/XRT observation in WT mode, fit the spectrum in XSPEC with the model "TBabs * diskbb" as a function of the size of the excised inner hole when extracting the source with an annular region. Departures from a roughly constant best-fit "Tin" when the inner hole extends close to the source center will inform us where pile-up is an issue. (See Romano et al 2006).

NOTES:
------
- Must initialize HEAsoft
- These routines are run in the script: run_pileup_analysis.py
'''

#----------------------------------------------------------------------------------------------------
# Fit the spectrum in XSPEC with the model "TBabs * diskbb"
def spectral_distortions():
    '''
    '''

    return

#----------------------------------------------------------------------------------------------------
# Thaw a parameter by setting its delta > 0
def thaw(x, delta=0.0001):
    dx = np.abs(x * delta)  # Absolute value to avoid freezing parameters that can have negative values
    return dx

#----------------------------------------------------------------------------------------------------
# Write the spectral distortion best-fit results for all Rcut values to an HDF5 file
def write_spectral_distortion(fout, fgrp_list, fsrc_list, fbkg_list, Rcut, cntrate0=100.0, \
        E_min=0.5, E_max=10.0, NE=1000, nHGal=0.74, phacuts=None):
    '''
    fout       - HDF5 output filename containing the results of the spectral distortion analysis
    froot      -
    fgrp_list  - Grouped spectrum file list <-- XSPEC ready
    fsrc_list  - Source spectrum file list
    fbkg_list  - Background spectrum file list
    Rcut       - Array of inner excised radii values corresponding to the source/background lists [pixels]
    phacuts    - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    '''
    # Output directory
    dout = ('/').join(fout.split('/')[0:-1])
    
    # Loop through each region (i.e., Rcut value)...
    #...and compute the background-subtracted count rate (w/ stdev)
    Nreg    = len(Rcut)
    cntrate = np.zeros(Nreg)
    stdrate = np.zeros(Nreg)
    for i in np.arange(Nreg):
        fsrc = fsrc_list[i]
        fbkg = fbkg_list[i]
        cntrate[i], stdrate[i] = bkgsub_count_rate(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts, stdev=True)

    # Find the region with the smallest Rcut that does not exceed the max count rate for pile-up
    iRcut0 = np.where(cntrate <= cntrate0)[0][0]
    Rcut0  = Rcut[iRcut0]

    # Fit the spectrum corresponding to the above region with nH free to determine its value
    fgrp0 = fgrp_list[iRcut0]
    fgrp0_models = xspec_models(fdata=fgrp0, E_min=E_min, E_max=E_max, NE=NE)
    froot        = dout + '/' + fgrp0.split('/')[-1].split('.pha')[0]
    vals_nH      = [nHGal, thaw(1.0)]
    vals_Tin     = [1.0,   thaw(1.0)]
    vals_normbb  = [1.0,   thaw(1.0)]
    fgrp0_models.TBabs_diskbb(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb)

    # Collect nH from the output best-fit parameter file
    quit()

    # Loop through each region (i.e., Rcut value)...
    #...fit the spectrum (nH fixed) and retrieve the best-fit values...
    nH        = np.zeros(Nreg)
    Tin       = np.zeros(Nreg)
    TinErr    = np.zeros(Nreg)
    normbb    = np.zeros(Nreg)
    normbbErr = np.zeros(Nreg)
    chisq     = np.zeros(Nreg)
    dof       = np.zeros(Nreg)
    for i in np.arange(Nreg):
        # Fit spectrum
        nH[i], Tin[i], TinErr[i], normbb[i], normbbErr[i], chisq[i], dof[i] = \
            fit_spectrum(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts)


    # FIND THE SPECTRUM WITH A COUNT RATE BELOW THRESHOLD

    # Write out the results to an HDF5 file
    f = h5py.File(fout, 'w')
    f.create_dataset('Rcut',      data=Rcut)       # [pixels]
    f.create_dataset('nH',        data=nH)         # [atoms cm^-2]
    f.create_dataset('Tin',       data=Tin)        # [keV]
    f.create_dataset('TinErr',    data=TinErr)     # [keV]
    f.create_dataset('normbb',    data=normbb)     # [...]
    f.create_dataset('normbbErr', data=normbbErr)  # [...]
    f.create_dataset('chisq',     data=chisq)      # [-]
    f.create_dataset('dof',       data=dof)        # [-]
    f.create_dataset('cntrate',   data=cntrate)    # [counts/second]
    f.create_dataset('stdrate',   data=stdrate)    # [counts/second]
    f.close()

#----------------------------------------------------------------------------------------------------

