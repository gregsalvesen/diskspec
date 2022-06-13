import numpy as np
import subprocess
import os
import h5py
from xspec import *
import xspec_emcee

'''
NAME:
-----
emcee_kerrbb.py

PURPOSE:
--------
Use EMCEE to fit the GRO J1655-40 spectrum with the model "TBabs*kerrbb" w/ nH free.

OUTPUTS:
--------
*** To the directory /Users/salvesen/research/fcolabs/analysis/xspec/ ***
xspec_emcee_kerrbb.xcm         - XCM file used to initialize emcee
xspec_emcee_kerrbb.pars        - Parameter file used to initialize the walkers
xspec_walkers_burn_kerrbb.txt  - File initializing walkers for the burn-in
xspec_walkers_emcee_kerrbb.txt - File initializing walkers for the production run
xspec_chain_burn_kerrbb.chain  - Chain file for the burn-in, used for restarts
xspec_chain_emcee_kerrbb.chain - Chain file for the production run, used for restarts

*** To the directory /Users/salvesen/research/fcolabs/results/ ***
xspec_chain_burn_kerrbb.hdf5   - Chain file for the burn-in, used for analysis
xspec_chain_emcee_kerrbb.hdf5  - Chain file for the production run, used for analysis

NOTES:
------
You must initialize HEAsoft before importing Xspec (command below does work)
subprocess.call(". $HEADAS/headas-init.sh", shell=True)
'''

#====================================================================================================
# DEFINE SOME VARIABLES

# EMCEE inputs
nwalk  = 100#100   # Number of walkers
nburn  = 1000#1000  # Number of iterations, or steps, for each walker for the burn-in
nmcmc  = 1000#1000  # Number of iterations, or steps, for each walker for the production run

# Input/Output files
ftag   = "iBin_"
obsID  = "00030009021"
dout   = "/Users/salvesen/research/fcolabs/results/" + obsID + "/"  # Output directory
fchi   = dout + "TBabs_kerrbb_" + ftag + obsID + ".xcm"      # Input best-fit parameters from chi^2 method
#fchi   = dout + "TBabs_kerrbb_iBin_" + obsID + ".xcm"      # Input best-fit parameters from chi^2 method
fxcm   = dout + ftag + "emcee_TBabs_kerrbb" + "_" + obsID + ".xcm"  # XCM file used to initialize emcee

fprior = dout + ftag + "priors"       + "_" + obsID + ".txt"    # Text file specifying the parameter priors
fwalk0 = dout + ftag + "walkers_burn" + "_" + obsID + ".txt"    # Parameter file used to initialize the walkers
fwalk  = dout + ftag + "walkers_mcmc" + "_" + obsID + ".txt"    # Parameter file used to re-initialize the walkers
fburn  = dout + ftag + "chain_burn"   + "_" + obsID + ".hdf5"   # Burn-in chain file (for analysis)
fburnR = dout + ftag + "chain_burn"   + "_" + obsID + ".chain"  # Burn-in chain file (for restarts)
fmcmc  = dout + ftag + "chain_mcmc"   + "_" + obsID + ".hdf5"   # Production run chain file (for analysis)
fmcmcR = dout + ftag + "chain_mcmc"   + "_" + obsID + ".chain"  # Production run chain file (for restarts)

# GRO J1655-40 known parametes and their 1-sigma uncertainties for EMCEE priors
mean_i,   sigma_i   = 68.65, 0.91#85.0,  2.0   # Adopt a Gaussian prior (Alternative: 68.65, 0.91)
mean_Mbh, sigma_Mbh = 5.40,  0.18  # Adopt a Gaussian prior
mean_Dbh, sigma_Dbh = 3.20,  0.2   # Adopt a Gaussian prior

# GRO J1655-40 parameter ranges
nH_min,   nH_max   = 0.0, 10.0
a_min,    a_max    = -1.0, 0.9999
i_min,    i_max    = 0.0, 85.0
Mbh_min,  Mbh_max  = 0.0, 10.0
Mdd_min,  Mdd_max  = 0.0, 50.0
Dbh_min,  Dbh_max  = 0.0, 10.0
hd_min,   hd_max   = 1.0, 10.0
norm_min, norm_max = 0.0, 10.0

#====================================================================================================
# LOAD IN THE DATA AND MODEL, FREE PARAMETERS AND SET THEIR RANGES

# Load the best fit parameter file returned by the script: xspec_kerrbb.sl
Xset.restore(fchi)

# Model components
m1    = AllModels(1)
comp1 = m1.TBabs
comp2 = m1.kerrbb

# Thaw parameters
comp1.nH.frozen   = False
comp2.a.frozen    = False
comp2.i.frozen    = False
comp2.Mbh.frozen  = False
comp2.Mdd.frozen  = False
comp2.Dbh.frozen  = False
comp2.hd.frozen   = False
comp2.norm.frozen = False

# Restrict parameter ranges: [val, delta, min, bot, top, max]
comp1.nH.values   = ', , ' + str(nH_min)   + ',' + str(nH_min)   + ',' + str(nH_max)   + ',' + str(nH_max)
comp2.a.values    = ', , ' + str(a_min)    + ',' + str(a_min)    + ',' + str(a_max)    + ',' + str(a_max)
comp2.i.values    = ', , ' + str(i_min)    + ',' + str(i_min)    + ',' + str(i_max)    + ',' + str(i_max)
comp2.Mbh.values  = ', , ' + str(Mbh_min)  + ',' + str(Mbh_min)  + ',' + str(Mbh_max)  + ',' + str(Mbh_max)
comp2.Mdd.values  = ', , ' + str(Mdd_min)  + ',' + str(Mdd_min)  + ',' + str(Mdd_max)  + ',' + str(Mdd_max)
comp2.Dbh.values  = ', , ' + str(Dbh_min)  + ',' + str(Dbh_min)  + ',' + str(Dbh_max)  + ',' + str(Dbh_max)
comp2.hd.values   = ', , ' + str(hd_min)   + ',' + str(hd_min)   + ',' + str(hd_max)   + ',' + str(hd_max)
comp2.norm.values = ', , ' + str(norm_min) + ',' + str(norm_min) + ',' + str(norm_max) + ',' + str(norm_max)

'''
# Specify Gaussian priors for Mbh, Dbh
try: os.remove(fprior)  # Remove a pre-existing priors file
except OSError: pass
paridx   = [comp2.Mbh.index, comp2.Dbh.index]
parmean  = [mean_Mbh,  mean_Dbh]
parsigma = [sigma_Mbh, sigma_Dbh]
nprior   = len(paridx)
priors   = [[paridx[i], parmean[i], parsigma[i]] for i in range(nprior)]
# Write priors out to a text file that will be read in by xspec_emcee
with open(fprior, 'w') as f:
    f.write('\n'.join('%s %s %s' % tuple(x) for x in priors))
f.close()
'''

# Specify Gaussian priors for i, Mbh, Dbh
try: os.remove(fprior)  # Remove a pre-existing priors file
except OSError: pass
paridx   = [comp2.i.index, comp2.Mbh.index, comp2.Dbh.index]
parmean  = [mean_i,  mean_Mbh,  mean_Dbh]
parsigma = [sigma_i, sigma_Mbh, sigma_Dbh]
nprior   = len(paridx)
priors   = [[paridx[i], parmean[i], parsigma[i]] for i in range(nprior)]
# Write priors out to a text file that will be read in by xspec_emcee
with open(fprior, 'w') as f:
    f.write('\n'.join('%s %s %s' % tuple(x) for x in priors))
f.close()

# Choose Wilms et al. (2000) abundances
Xset.abund = "wilm"  # <-- Careful! 'wilms' is not acceptable

# Save the parameter file used to initialize the walkers
try: os.remove(fxcm)  # Remove a pre-existing XCM file
except OSError: pass
Xset.save(fxcm, info='a')

#====================================================================================================
# BURN-IN

# Initialize the walkers in a "tiny Gaussian ball" around the best-fit parameters: [nH, a, i, Mbh, Mdd, Dbh, hd, norm]
pars0    = [comp1.nH.values[0], comp2.a.values[0], comp2.i.values[0], comp2.Mbh.values[0], comp2.Mdd.values[0], comp2.Dbh.values[0], comp2.hd.values[0], comp2.norm.values[0]]
npars    = len(pars0)
walkers0 = [pars0 + 1e-4*np.random.randn(npars) for i in range(nwalk)]
with open(fwalk0, 'w') as f:
    f.write('\n'.join('%s %s %s %s %s %s %s %s' % tuple(x) for x in walkers0))
f.close()

# Remove any previous burn-in chain files of the same name
try: os.remove(fburn)
except OSError: pass
try: os.remove(fburnR)
except OSError: pass

# Burn-in (doing this my way so I have control, rather than default way)
cmd_burn = "python xspec_emcee.py " + fxcm + \
    " --initial-parameters " + fwalk0 + \
    " --nwalkers " + str(nwalk) + \
    " --nburn 0" + \
    " --niters " + str(nburn) + \
    " --output-hdf5 " + fburn + \
    " --output-chain " + fburnR + \
    " --parameter-priors " + fprior
subprocess.call(cmd_burn, shell=True)

#====================================================================================================
# PRODUCTION RUN

# Considering the last 10% of the burn-in chain...
# ...re-initialize the walkers about the maximal likelihood parameter set...
# ...with a spread in each parameter equal to its standard deviation.
# (This neglects any covariances between parameters that EMCEE has uncovered so far, but w/e)
# (We could return the covariance matrix and use this!)

# Starting index for the last 10% of the burn-in chain
i0 = int(nburn * 0.9)

# Collect the last 10% of the burn-in chains and lnprobs
f = h5py.File(fburn, 'r')
chains  = f['chain'][:,i0:,:]   # [walker, iteration, parameter]
lnprobs = f['lnprob'][:,i0:]    # [walker, iteration]
f.close()

# Find the maximum likelihood parameter set
ijLmax   = np.unravel_index(lnprobs.argmax(), lnprobs.shape)  # Maximum likelihood index (i,j)
parsLmax = chains[ijLmax]                                     # Maximum likelihood parameter set

# Calculate parameter standard deviations
std_nH   = np.std(chains[:,:,0].flatten())
std_a    = np.std(chains[:,:,1].flatten())
std_i    = np.std(chains[:,:,2].flatten())
std_Mbh  = np.std(chains[:,:,3].flatten())
std_Mdd  = np.std(chains[:,:,4].flatten())
std_Dbh  = np.std(chains[:,:,5].flatten())
std_hd   = np.std(chains[:,:,6].flatten())
std_norm = np.std(chains[:,:,7].flatten())
std_pars = [std_nH, std_a, std_i, std_Mbh, std_Mdd, std_Dbh, std_hd, std_norm]

# Re-initialize the walkers
walkers = [np.random.normal(loc=parsLmax, scale=std_pars) for i in range(nwalk)]
with open(fwalk, 'w') as f:
    f.write('\n'.join('%s %s %s %s %s %s %s %s' % tuple(x) for x in walkers))
f.close()

# Prodution run
cmd_mcmc = "python xspec_emcee.py " + fxcm + \
    " --initial-parameters " + fwalk + \
    " --nwalkers " + str(nwalk) + \
    " --nburn 0" + \
    " --niters " + str(nmcmc) + \
    " --output-hdf5 " + fmcmc + \
    " --output-chain " + fmcmcR + \
    " --parameter-priors " + fprior
subprocess.call(cmd_mcmc, shell=True)

