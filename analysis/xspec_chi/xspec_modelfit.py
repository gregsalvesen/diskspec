import numpy as np
import subprocess
import os
import h5py
from xspec import *
from pvalue2sigma import *

'''
NAME:
-----
xspec_modelfit

PURPOSE:
--------
Fit a spectrum with an Xspec model specified by the user.

INPUTS:
-------
theModel - Xspec model to use for spectral fitting (string; e.g., "TBabs*kerrbb")
setup    - (see below)
froot    - Root filename for the outputs
fdata    - Data filename (.pha)
E_min    - Minimum photon energy (keV)
E_max    - Maximum photon energy (keV)
NE       - Number of energy bins to be logarithmically spaced
dfstat   - Delta fit statistic to use for the Xspec "error" command
link     - (see below)

"setup" must be a list of tuples...
...where each tuple specifies a model component (1st element)...
...and a dictionary of initial parameter names/values (2nd element)
"setup" Format: [('componentName', {'parameter1':values1, 'parameter2':values2, ...}), ...]
Example:
>> theModel = "TBabs*diskbb"
>> setup    = [("TBabs", {'nH':0.74}), ("diskbb", {'Tin':1.0, 'norm':1.0})]

"link" must be a list of tuples, where each tuple has 4-elements as follows...
... 1st, 2nd element --> component name, parameter name being linked FROM
... 3rd, 4th element --> component name, parameter name being linked TO
"link" Format: [('componentLinkFrom', 'parameterLinkFrom', 'componentLinkTo', 'parameterLinkTo'), ...]
Example:
>> theModel = "TBabs*compTT"
>> link     = [('compTT', 'T0', 'diskbb', 'Tin')]  # <-- This links compTT.T0 to diskbb.Tin


OUTPUTS:
--------
froot.init - Initial parameter file
froot.fit  - Best-fit parameter file
froot.xcm  - Data and model configuration --> output of Xset.save(info='a')
froot.hdf5 - HDF5 file containing the best-fit model spectrum and data spectrum

NOTES:
------
- Must initialize HEAsoft before importing Xspec
- Wilms et al. (2000) abundances are set if 'TBabs' is a model component
- The choice for abundances turns out to make a HUGE difference in the fitted nH results.
- We do not extend the energy array <-- DO THIS IF YOU WANT TO CALCULATE MODEL FLUXES
'''

#====================================================================================================
def xspec_modelfit(theModel, setup, froot, fdata, E_min=0.5, E_max=10.0, NE=1000, dfstat=1.0, link=None):

    # Output filenames
    finit = froot + '.init'
    fbfit = froot + '.fit'
    fxcm  = froot + '.xcm'
    fspec = froot + '.hdf5'

    #----------------------------------------------------------------------------------------------------
    # LOAD IN THE DATASET AND SET XSPEC PREFERENCES

    # Set the Xspec chatter level to zero
    XspecSettings.chatter = 0

    # Load in the data
    Spectrum(fdata)

    # Ignore bad data
    AllData.ignore('bad')

    #----------------------------------------------------------------------------------------------------
    # DEFINE THE MODEL

    # Define the model
    Model(theModel)
    m1 = AllModels(1)

    # Adjust the energy range for the model according to the input Emin, Emax, NE
    listEnergies = [E_min, E_max, NE, 'log', ]
    strEnergies  = [str(item) for item in listEnergies]
    setEnergies  = ' '.join(strEnergies)
    AllModels.setEnergies(setEnergies)

    # Choose Wilms et al. (2000) abundances if "TBabs" is a model component
    if ('TBabs' in m1.componentNames): Xset.abund = 'wilm'  # <-- Careful: 'wilms' will not work!

    #----------------------------------------------------------------------------------------------------
    # SET THE MODEL PARAMETERS

    # Initialize a list that will store thawed parameter indices
    thawed = []

    # Loop through components in the model
    compsAll  = [item[0] for item in setup]  # Component names from user input
    compNames = m1.componentNames            # Component names from Xspec
    Ncomps    = np.size(compNames)           # Number of components
    for i in range(Ncomps):
    
        # Component "i" object
        comp = m1.__getattribute__(compNames[i])
        
        # Locate the index of component "i" in compsAll
        icomp = compsAll.index(comp.name)

        # Loop through parameters in component "i"
        parsDict = setup[icomp][1]                          # Parameter names/vals dictionary for component "i"
        parNames = comp.__getattribute__('parameterNames')  # Parameter names for component "i" from Xspec
        Npars    = np.size(parNames)                        # Number of parameters for componet "i"
        for j in range(Npars):

            # Parameter "j" object
            par = comp.__getattribute__(parNames[j])

            # For Parameter "j", set [value, delta, min, bot, top, max]
            # Note: Delta sets whether the parameter is frozen (<0) or thawed (>0)
            par.values = parsDict[par.name]
        
            # Collect list of thawed parameter indices
            if (par.frozen is False): thawed.append(par.index)

    # Link/tie parameters to one another, if specified
    if (link is not None):
        # Loop through parameter links
        Nlinks = len(link)  # Number of links
        for i in range(Nlinks):
            plinkFrom = m1.__getattribute__(link[i][0]).__getattribute__(link[i][1])  # Parameter object to link from
            plinkTo   = m1.__getattribute__(link[i][2]).__getattribute__(link[i][3])  # Parameter object to link to
            plinkFrom.link = str(plinkTo.index)

    # Save the initial parameters file
    try: os.remove(finit)
    except OSError: pass
    Xset.save(finit, info='m')

    #----------------------------------------------------------------------------------------------------
    # PERFORM THE FIT

    # Fit the spectrum
    Fit.statMethod  = 'chi'
    Fit.statTest    = 'chi'
    Fit.nIterations = 1000
    Fit.query       = 'yes'
    Fit.perform()

    # Loop through thawed parameters and run the Xspec "error" command
    for ithaw in thawed:
        Fit.error('maximum 10 ' + str(dfstat) + ' ' + str(ithaw))

    # Collect fit results
    chisq  = Fit.statistic
    dof    = Fit.dof
    pvalue = Fit.nullhyp
    sigma  = sigma_chisq(pvalue=pvalue, df=1)
    covar  = Fit.covariance

    # Save the data + model + fit as an XCM file
    try: os.remove(fxcm)
    except OSError: pass
    Xset.save(fxcm, info='a')
    
    #----------------------------------------------------------------------------------------------------
    # SAVE THE BESTFIT: PARAMETERS, 1-SIGMA UNCERTAINTIES, FIT STATISTIC, DOF, P-VALUE, SIGMA-VALUE

    # Write to the best-fit results file
    f    = open(fbfit, 'w')
    hdr1 = 'Model               ' + theModel                    + '\n'
    hdr2 = 'Reduced Chi-Squared ' + '{0:.5g}'.format(chisq/dof) + '\n'
    hdr3 = 'Chi-Squared         ' + '{0:.5g}'.format(chisq)     + '\n'
    hdr4 = 'Degrees of Freedom  ' + '{0:.5g}'.format(dof)       + '\n'
    hdr5 = 'Null Hypo p-value   ' + '{0:.5g}'.format(pvalue)    + '\n'
    hdr6 = 'Conf Interval sigma ' + '{0:.5g}'.format(sigma)     + '\n' + '\n'
    hdr7 = ['Comp', 'Par', 'Unit', 'Value', '-Err', '+Err']
    f.write(hdr1)
    f.write(hdr2)
    f.write(hdr3)
    f.write(hdr4)
    f.write(hdr5)
    f.write(hdr6)
    f.write('\t'.join(hdr7)+'\n')

    # Loop through each parameter, collect its best-fit value (and 1-sigma uncertainty if thawed)
    for i in range(Ncomps):
        
        # Component "i" object
        comp = m1.__getattribute__(compNames[i])
        
        # Locate the index of component "i" in compsAll
        icomp = compsAll.index(comp.name)

        # Loop through parameters in component "i"
        parNames = comp.__getattribute__('parameterNames')  # Parameter names for component "i" from Xspec
        Npars    = np.size(parNames)                        # Number of parameters for componet "i"
        for j in range(Npars):

            # Parameter "j" object
            par = comp.__getattribute__(parNames[j])

            # Collect the component name, paramter name, parameter unit
            line = [comp.name, par.name, par.unit]

            # For Parameter "j", collect the best-fit value
            bfval = par.values[0]
            line.append('{0:.5g}'.format(bfval))
            
            # For Parameter "j", collect the 1-sigma uncertainties (if thawed)
            if (par.index in thawed):
                errm = par.error[0] - bfval  # (-)
                errp = par.error[1] - bfval  # (+)
                line.append('{0:.3g}'.format(errm))
                line.append('+' + '{0:.3g}'.format(errp))

            # Write the line to the best-fit results file
            f.write('\t'.join(line)+'\n')

    # Close the file
    f.close()

    #----------------------------------------------------------------------------------------------------
    # SAVE THE SPECTRUM: DATA, BESTFIT, RESIDUALS, RATIO, MODEL (PHYSICAL UNITS)

    # Set the x-axis to keV
    Plot.xAxis = 'keV'
    
    # Collect the data and bestfit [normalized counts s^-1 keV^-1]
    Plot('data')
    E_data     = np.array(Plot.x())      # [keV]
    FE_data    = np.array(Plot.y())      # [normalized counts s^-1 keV^-1]
    Eerr_data  = np.array(Plot.xErr())   # [keV]
    FEerr_data = np.array(Plot.yErr())   # [normalized counts s^-1 keV^-1]
    FE_bestfit = np.array(Plot.model())  # [normalized counts s^-1 keV^-1]

    # Collect the data - model residuals
    FE_resid    = FE_data - FE_bestfit  # [normalized counts s^-1 keV^-1]
    FEerr_resid = FEerr_data            # [normalized counts s^-1 keV^-1]
    
    # Collect the data / model ratio
    FE_ratio    = FE_data / FE_bestfit
    FEerr_ratio = FEerr_data / FE_bestfit

    # Collect the model in physical units [photons s^-1 cm^-2 keV^-1]
    Plot('model')
    E_model   = np.array(Plot.x())           # [keV]
    FE_model  = np.array(Plot.model())       # [photons s^-1 cm^-2 keV^-1]
    keV2erg   = 1.602176565e-09              # [erg keV^-1]
    EFE_model = E_model * FE_model * keV2erg # [erg s^-1 cm^-2 keV^-1]

    # Write results to an HDF5 file
    f = h5py.File(fspec, 'w')
    f.create_dataset('Energy',    data=E_data)       # [keV]
    f.create_dataset('Data',      data=FE_data)      # [normalized counts s^-1 keV^-1]
    f.create_dataset('EnergyErr', data=Eerr_data)    # [keV]
    f.create_dataset('DataErr',   data=FEerr_data)   # [normalized counts s^-1 keV^-1]
    f.create_dataset('BestFit',   data=FE_bestfit)   # [normalized counts s^-1 keV^-1]
    f.create_dataset('Resid',     data=FE_resid)     # [normalized counts s^-1 keV^-1]
    f.create_dataset('ResidErr',  data=FEerr_resid)  # [normalized counts s^-1 keV^-1]
    f.create_dataset('Ratio',     data=FE_ratio)     # [-]
    f.create_dataset('RatioErr',  data=FEerr_ratio)  # [-]
    f.create_dataset('Emodel',    data=E_model)      # [keV]
    f.create_dataset('FEmodel',   data=FE_model)     # [photons s^-1 cm^-2 keV^-1]
    f.create_dataset('EFEmodel',  data=EFE_model)    # [erg s^-1 cm^-2 keV^-1]
    f.close()
    
    # Clear the data and model <-- This is important!
    AllData.clear()
    AllModels.clear()

