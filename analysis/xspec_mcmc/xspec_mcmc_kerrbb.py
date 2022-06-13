import numpy as np
import subprocess
import os
from xspec import *

'''
NAME:
-----
xspec_mcmc_kerrbb.py

PURPOSE:
--------
Use emcee to fit the GRO J1655-40 spectrum with the model "TBabs*kerrbb" w/ nH frozen.

OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
--------
xspec_mcmc_kerrbb.xcm        - parameter file used to initialize the walkers
xspec_chain_burn_kerrbb.hdf5 - burn-in chain file
xspec_chain_mcmc_kerrbb.hdf5 - emcee chain file

NOTES:
------
Must initialize HEAsoft before importing Xspec...but the command below does not do the job
subprocess.call(". $HEADAS/headas-init.sh", shell=True)
'''

#====================================================================================================
# DEFINE SOME VARIABLES

# MCMC inputs
nsteps = 1000
nwalk  = 100                  # Number of walkers
nburn  = int(nwalk * nsteps)  # Number of iterations, or steps, for each walker for the burn-in
nmcmc  = int(nwalk * nsteps)  # Number of iterations, or steps, for each walker for the production run

# Input/Output files
pre   = "Nwalk_" + str(nwalk) + "_"
dout  = "/Users/salvesen/research/fcolabs/results/"
fchi  = dout + "xspec_kerrbb.xcm"                           # Best-fit parameters from chi^2 method
fxcm  = dout + pre + "xspec_mcmc_kerrbb.xcm"         # Parameter file used to initialize the walkers
fburn = dout + pre + "xspec_chain_burn_kerrbb.fits"  # Burn-in chain file
fmcmc = dout + pre + "xspec_chain_mcmc_kerrbb.fits"  # Burn-in chain file

# GRO J1655-40 known parametes and their 1-sigma uncertainties for MCMC priors
mean_a              = 0.92         # Adopt a uniform prior
mean_i,   sigma_i   = 68.65, 0.91  # Adopt a Gaussian prior
mean_Mbh, sigma_Mbh = 5.40,  0.18  # Adopt a Gaussian prior
mean_Dbh, sigma_Dbh = 3.20,  0.2   # Adopt a Gaussian prior

# GRO J1655-40 parameter ranges
a_min,   a_max   = 0.9, 0.9999
i_min,   i_max   = 0.0, 85.0
Mbh_min, Mbh_max = 0.0, 50.0
Mdd_min, Mdd_max = 0.0, 10.0
Dbh_min, Dbh_max = 0.0, 50.0
hd_min,  hd_max  = 1.0, 10.0

#====================================================================================================
# LOAD IN THE DATASET

# Load the best fit parameter file returned by the script: xspec_kerrbb.sl
Xset.restore(fchi)

# Model components
m1    = AllModels(1)
comp1 = m1.TBabs
comp2 = m1.kerrbb

# Parameters: [val, delta, min, bot, top, max]

# Thaw parameters and restrict their range
comp2.a.frozen   = False
comp2.i.frozen   = False
comp2.Mbh.frozen = False
comp2.Mdd.frozen = False
comp2.Dbh.frozen = False
comp2.hd.frozen  = False
comp2.a.values   = ', , ' + str(a_min)   + ',' + str(a_min)   + ',' + str(a_max)   + ',' + str(a_max)
comp2.i.values   = ', , ' + str(i_min)   + ',' + str(i_min)   + ',' + str(i_max)   + ',' + str(i_max)
comp2.Mbh.values = ', , ' + str(Mbh_min) + ',' + str(Mbh_min) + ',' + str(Mbh_max) + ',' + str(Mbh_max)
comp2.Mdd.values = ', , ' + str(Mdd_min) + ',' + str(Mdd_min) + ',' + str(Mdd_max) + ',' + str(Mdd_max)
comp2.Dbh.values = ', , ' + str(Dbh_min) + ',' + str(Dbh_min) + ',' + str(Dbh_max) + ',' + str(Dbh_max)
comp2.hd.values  = ', , ' + str(hd_min)  + ',' + str(hd_min)  + ',' + str(hd_max)  + ',' + str(hd_max)

# Specify Gaussian priors for i, Mbh, Dbh
comp2.i.prior   = 'gauss, ' + str(mean_i)   + ', ' + str(sigma_i)
comp2.Mbh.prior = 'gauss, ' + str(mean_Mbh) + ', ' + str(sigma_Mbh)
comp2.Dbh.prior = 'gauss, ' + str(mean_Dbh) + ', ' + str(sigma_Dbh)

'''
# Initialize the walkers in a "tiny Gaussian ball" around the best fit result
a_walk0   = [comp2.a.values[0]   + 1e-4*np.random.randn() for i in range(nwalk)]
i_walk0   = [comp2.i.values[0]   + 1e-4*np.random.randn() for i in range(nwalk)]
Mbh_walk0 = [comp2.Mbh.values[0] + 1e-4*np.random.randn() for i in range(nwalk)]
Mdd_walk0 = [comp2.Mdd.values[0] + 1e-4*np.random.randn() for i in range(nwalk)]
Dbh_walk0 = [comp2.Dbh.values[0] + 1e-4*np.random.randn() for i in range(nwalk)]
hd_walk0  = [comp2.hd.values[0]  + 1e-4*np.random.randn() for i in range(nwalk)]
'''

# Remove any previous chain file of the same name
try: os.remove(fburn)
except OSError: pass

# Create a chain object and run!
c1 = Chain(fileName=fburn, burn=0, runLength=nburn, walkers=nwalk, rand=False, algorithm='gw')

#
#bfpars = AllChains.best()

#
#AllChains.dic()

# Need to remember that we have more free params now when calculating the reduced chi2

# Get the 1D distributions
#AllChains.margin('3, 0.9, 0.9999, 10')
#a_arr  = AllChains.marginResults(3)
#a_prob = AllChains.marginResults('probability')
#a_frac = AllChains.marginResults('fraction')

#The chain stat command writes out the repeat fraction.
#a_stat = AllChains.stat(comp2.a.index)  # <-- use this index thing to reference pars

# Save the parameter file used to initialize the walkers
try: os.remove(fxcm)  # Remove a pre-existing XCM file
except OSError: pass
Xset.save(fxcm)

# Burn-in xspec_emcee()
#xspec_emcee.emcee(xcm=fxcm, nwalkers=nwalk, nburn=0, niters=nburn, outhdf5=fburn)
