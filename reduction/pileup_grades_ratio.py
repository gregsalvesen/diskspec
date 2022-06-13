import numpy as np
import h5py
from Swift_XRT_WTmode import *

'''
PURPOSE:
--------
For a Swift/XRT observation in WT mode, calculate the counts ratio of event grade 0 / grades 0-2 as a function of the size of the excised inner hole when extracting the source with an annular region. Departures from a roughly constant grades ratio when the inner hole extends close to the source center will inform us where pile-up is an issue. (See Romano et al 2006).

NOTES:
------
- Must initialize HEAsoft
- These routines are run in the script: run_pileup_analysis.py
'''

#----------------------------------------------------------------------------------------------------
# Calculate the background-subtracted count fraction (grade 0) / (grade 0-2) w/ standard deviation as a function of excised inner annulus radius
def grades_ratio(fsrc0, fbkg0, fsrc02, fbkg02, phacuts=None):
    '''
    fsrc0   - Source spectrum file for grade 0
    fbkg0   - Background spectrum file for grade 0
    fsrc02  - Source spectrum file for grade 0-2
    fbkg02  - Background spectrum file for grade 0-2
    phacuts - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    '''
    # Collect the number of background-subtracted counts for grade 0 and 0-2
    Ncounts0,  stdcnts0  = bkgsub_counts(fsrc=fsrc0,  fbkg=fbkg0,  phacuts=phacuts)
    Ncounts02, stdcnts02 = bkgsub_counts(fsrc=fsrc02, fbkg=fbkg02, phacuts=phacuts)

    # Ratio of (grade 0) / (grade 0-2) background-subtracted counts (and standard deviation)
    frac = Ncounts0 / Ncounts02
    fstd = frac * np.sqrt((stdcnts0 / Ncounts0)**2 + (stdcnts02 / Ncounts02)**2)
    return frac, fstd

#----------------------------------------------------------------------------------------------------
# Write the grades ratio for each pair of source/background files (all Rcut values) to an HDF5 file
def write_grades_ratio(fout, fsrc0_list, fbkg0_list, fsrc02_list, fbkg02_list, Rcut, phacuts=None):
    '''
    fout        - HDF5 output filename containing the results of the event grades ratio analysis
    fsrc0_list  - Source spectrum file list for grade 0
    fbkg0_list  - Background spectrum file list for grade 0
    fsrc02_list - Source spectrum file list for grade 0-2
    fbkg02_list - Background spectrum file list for grade 0-2
    Rcut        - Array of inner excised radii values corresponding to the source/background lists [pixels]
    phacuts     - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    '''
    # Loop through each region (i.e., Rcut value)...
    #...compute the grades ratio (w/ stdev) and background-subtracted count rate (w/ stdev)
    Nreg    = len(Rcut)
    frac    = np.zeros(Nreg)
    fstd    = np.zeros(Nreg)
    cntrate = np.zeros(Nreg)
    stdrate = np.zeros(Nreg)
    for i in np.arange(Nreg):
        fsrc0   = fsrc0_list[i]
        fbkg0   = fbkg0_list[i]
        fsrc02  = fsrc02_list[i]
        fbkg02  = fbkg02_list[i]
        # Grades ratio
        frac[i], fstd[i] = grades_ratio(fsrc0=fsrc0, fbkg0=fbkg0, fsrc02=fsrc02, fbkg02=fbkg02, phacuts=phacuts)
        # Count rate [counts/second]
        cntrate[i], stdrate[i] = bkgsub_count_rate(fsrc=fsrc02, fbkg=fbkg02, phacuts=phacuts, stdev=True)

    # Write out the results to an HDF5 file
    f = h5py.File(fout, 'w')
    f.create_dataset('Rcut',                    data=Rcut)     # [pixels]
    f.create_dataset('grades_ratio',            data=frac)     # [-]
    f.create_dataset('grades_ratio_stdev',      data=fstd)     # [-]
    f.create_dataset('bkgsub_count_rate',       data=cntrate)  # [counts / second]
    f.create_dataset('bkgsub_count_rate_stdev', data=stdrate)  # [counts / second]
    f.close()

#----------------------------------------------------------------------------------------------------
