import numpy as np
import subprocess
from xspec import *

'''
NAME:
-----
kerrbb.py

PURPOSE:
--------
Fit the GRO J1655-40 spectrum with the model "TBabs*kerrbb" w/ nH frozen.
Only two free parameters: "Mdd" and "hd"

OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
--------
kerrbb.xcm - Model fit file

NOTES:
------
Must initialize HEAsoft before importing Xspec
'''

#====================================================================================================
# DEFINE SOME VARIABLES

# Output filenames
dout  = "/Users/salvesen/research/fcolabs/results/"
fpars = dout + "kerrbb.par"
fbfit = dout + "kerrbb.fit"
fxcm  = dout + "kerrbb.xcm"

# Data file
fdata = "source_annulus_grp10.pha"

# Energy range for fitting
Emin = 0.5   # [keV]
Emax = 10.0  # [keV]
NE   = 100

# Galactic value for nH (Dickey & Lockman 1990)
nHGal = 0.74  # [10^22 atoms cm^-2]

# GROJ1655-40 kerrbb parameters
norm  = 1.0    # [-]
eta   = 0.0    # [-]
a     = 0.92   # [-]
i     = 68.65  # [degrees]
Mbh   = 5.4    # [M_sun]
Mdd   = 1.0    # <-- Unknown (free) [10^18 g s^-1]
Dbh   = 3.2    # [kpc]
hd    = 1.7    # <-- Unknown (free) [-]
rflag = 1      # Self-irradiation flag
lflag = 1      # Limb-darkening flag

#====================================================================================================
# LOAD IN THE DATASET

# Set the Xspec chatter level to zero
XspecSettings.chatter = 0

# Load in the data
Spectrum(fdata)

# Ignore bad data
AllData.ignore("bad")

#====================================================================================================
# FIT THE MODEL "TBabs*kerrbb"

# Choose Wilms et al. (2000) abundances
Xset.abund = "wilm"  # <-- Careful! 'wilms' is not acceptable

# Define the model as TBabs*kerrbb
Model("TBabs*kerrbb")
m1 = AllModels(1)

# Adjust the energy range for the model according to the input Emin, Emax, NE
setEnergies = str(Emin) + " " + str(Emax) + " " + str(NE) + " " + "log"
AllModels.setEnergies(setEnergies)

# Model components
comp1 = m1.TBabs
comp2 = m1.kerrbb

# Set the prefit KERRBB parameters
comp1.nH.values    = nHGal
comp2.eta.values   = eta
comp2.a.values     = a
comp2.i.values     = i
comp2.Mbh.values   = Mbh
comp2.Mdd.values   = Mdd
comp2.Dbh.values   = Dbh
comp2.hd.values    = hd
comp2.rflag.values = rflag
comp2.lflag.values = lflag
comp2.norm.values  = norm

# Freeze the static system parameters
comp1.nH.frozen    = False
comp2.eta.frozen   = True
comp2.a.frozen     = True
comp2.i.frozen     = True
comp2.Mbh.frozen   = True
comp2.Mdd.frozen   = False
comp2.Dbh.frozen   = True
comp2.hd.frozen    = False
comp2.rflag.frozen = True
comp2.lflag.frozen = True
comp2.norm.frozen  = True


# Fit the spectrum
Fit.statMethod  = "chi"
Fit.nIterations = 1000
Fit.perform()

# Save the data + model + fit as an XCM file
subprocess.call("rm " + fxcm, shell=True)  # Remove a pre-existing XCM file
Xset.save(fxcm)
