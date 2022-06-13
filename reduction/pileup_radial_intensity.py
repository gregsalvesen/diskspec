import numpy as np
import h5py
from Swift_XRT_WTmode import *
from scipy.optimize import curve_fit
from astropy.io import fits

# PROBLEM:
# --------
# THE CALIBRATED KING MODEL PSF DOES A TERRIBLE JOB FITTING THE DATA.
# REPLACING beta WITH beta/2 GIVES A BETTER FIT, BUT THEN THE PSF CORE DOES NOT APPEAR PILED-UP, BUT IT SHOULD BE...
#
# IF YOU GET THIS WORKING, USE ENERGY BINS: [0.5,1], [1,2], [2,10] keV
#
# PLOT AND FIT THE INTENSITY PROFILE ON BOTH SIDES OF THE SOURCE, DO NOT COMBINE +/- RADII (WE ARE CURRENTLY DOUBLE COUNTING...NEED TO FIX THIS)

'''
PURPOSE:
--------
For a Swift/XRT observation in WT mode, create a radial intensity profile made up from source-centered annuli of 1-pixel width. Fit the outer wings of this observed profile with the model PSF (King function) and extrapolate to the PSF core region. Departures away from the model PSF in the inner regions will inform us where pile-up is an issue. (See Moretti et al. 2005 and Vaughan et al 2006).

NOTES:
------
- Must initialize HEAsoft
- These routines are run in the script: run_pileup_analysis.py
- XRT WT PSF coefficients: https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/data/swift/xrt/index.html
- Is it "better" to make images rather than spectra to retrieve the count rate?
- I confirmed with the Swift Help Desk that the PSF calibration coeeficients yield r_C in [pixels]
- Do we ever see an excess in intensity in the outer wings of the PSF? <-- JMM
'''

#----------------------------------------------------------------------------------------------------
# Return the radial intensity profile of the observed PSF...
# ...as well as the energy grid and off-axis angle needed to compute the model PSF for data fitting
def radial_intensity(fsrc_list, fbkg_list, Rin, Rout, phacuts, dE=None):
    '''
    fsrc_list -
    fbkg_list -
    Rin       -
    Rout      -
    phacuts   -
    dE        -  [keV]
    '''
    # Off-axis angle of the source
    theta = offaxis_angle(fitsfile=fsrc_list[0])  # [degrees]

    # Energy grid [keV]
    Elo   = phacuts[0] / 100.0         # [keV]
    Ehi   = phacuts[1] / 100.0         # [keV]
    if (dE is None): NE = 2
    else:            NE = int((Ehi - Elo) / dE + 1)
    Egrid = np.linspace(Elo, Ehi, NE)  # [keV]

    # PHA channel grid (1 XRT channel = 10 eV)
    Npha    = NE
    phagrid = np.linspace(phacuts[0], phacuts[1], Npha).astype(int)  # [channels]
    
    # Loop through each region (i.e., [Rin,Rout] values)...
    # ...and compute the radial intensity profile (w/ stdev) <-- i.e., count rate radial profile
    Nreg     = len(Rin)
    Rdata    = np.zeros(Nreg)
    PSFdata  = np.zeros(Nreg)
    PSFerr   = np.zeros(Nreg)
    Rdata    = np.zeros(Nreg)
    PSFdataE = np.zeros([Nreg, NE])
    PSFerrE  = np.zeros([Nreg, NE])
    for i in np.arange(Nreg):
        
        # Source and background spectral files for the current region
        fsrc = fsrc_list[i]
        fbkg = fbkg_list[i]
        
        # Radial mid-point location of the current region
        Rdata[i] = 0.5 * (Rin[i] + Rout[i])  # [pixels]

        # Bolometric (i.e., full range of phacuts) radial intensity profile of the observed PSF [counts/sec]
        PSFdata[i], PSFerr[i] = bkgsub_count_rate(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts, stdev=True)

        # Loop through each energy band
        for j in np.arange(NE-1):
            
            # PHA channel cuts for the current energy band
            phaband = [phagrid[j], phagrid[j+1]]  # [channels]
        
            # Build the radial intensity profile of the observed PSF [counts/sec] in current energy band
            PSFdataE[i,j], PSFerrE[i,j] = bkgsub_count_rate(fsrc=fsrc, fbkg=fbkg, phacuts=phaband, stdev=True)
            print i, j, PSFdataE[i,j]

    return Rdata, PSFdata, PSFerr, PSFdataE, PSFerrE, Egrid, theta

#----------------------------------------------------------------------------------------------------
# Fit a King function to the outer wings of the PSF profile and extrapolate to the core
def fitKing2PSF(Rdata, PSFdataE, Egrid, theta, PSFmin=10.0, PSFerrE=None):
    '''
    Rdata    - Radial locations corresponding to the PSF profile [pixels]
    PSFdataE  - Observed PSF radial profile (un-normalized) [counts/sec]
    Egrid    -
    theta    -
    cntrate0 - Count rate threshold below which to fit the PSF model to the observed PSF profile [counts/s]
    PSFerrE   - Uncertainties corresponding to the PSF profile (if None --> PSFerr = sqrt(PSF))
    '''
    # Use Poisson counting uncertainties if errors on the PSFdata are not supplied [counts/sec]
    if (PSFerrE is None): PSFerrE = np.sqrt(PSFdataE)

    # The inner core of the PSF cannot be used when fitting the PSF King model...
    # ...so find the radius interior to which we will ignore duing the fit...
    # ...which is decided to be where the count rate is less than some threshold cntrate0 [counts/sec].
    PSFdata = np.sum(PSFdataE, axis=1)  # <-- Collapsed the energy dimension [counts/sec]
    iR0     = np.where(PSFdata <= PSFmin)[0][0]
    Rwing   = Rdata[iR0]

    # Do not fit the data (Rdata, PSFdataE, PSFerrE) interior to Rwing --> i.e., only fit the PSF wing
    Rdata_wing    = Rdata[iR0:]
    PSFdataE_wing = PSFdataE[iR0:,:]
    PSFerrE_wing  = PSFerrE[iR0:,:]

    # Initialize the array that will store the best-fit model of the PSF radial profile
    NR     = len(Rdata)
    PSFfit = np.zeros(NR)

    # Loop through each energy bin and fit a King profile to the PSF
    NE = len(Egrid)
    for j in np.arange(NE-1):

        # Monochromatic energy bin (use the mid-point)
        Emono = 0.5 * (Egrid[j] + Egrid[j+1])

        # PSF King model core radius "Rc" and slope "beta"
        Rc0, beta0 = plane_function(E=Emono, theta=theta)
        
        # Amplitude initial guess for the King function model <-- choose A0 = PSF(R=Rwing)
        A0 = PSFdataE_wing[0,j] / ((1.0 + (Rwing / Rc0)**2)**(-1.0 * beta0))

        # King model for fitting, where the normalization is the only free parameter
        def kingfit(R, A):
            PSF = A * (1.0 + (R / Rc0)**2.0)**(-1.0 * beta0)
            return PSF
        
        # Fit a King profile to the PSF wing
        # --> The curve_fit() "f" parameter must be the name of a function...
        #     ...that takes the independent variable as the first argument...
        #     ...and the parameters to fit as separate remaining arguments
        xdata = Rdata_wing
        ydata = PSFdataE_wing[:,j]
        sigma = PSFerrE_wing[:,j]
        p0    = [A0]
        popt, pcov   = curve_fit(f=kingfit, xdata=xdata, ydata=ydata, p0=p0, sigma=sigma, bounds=(0, np.inf))
        perr         = np.sqrt(np.diag(pcov))
        PSFmodE      = king(R=Rdata, Rc=Rc0, beta=beta0, A=popt[0])
        PSFmodE_wing = PSFmodE[iR0:]
        chisq, dof   = chisqdof(ydata=ydata, ymod=PSFmodE_wing, Npars=len(p0))
        print "Reduced chi-square: ", chisq / dof, Rc0, beta0, popt[0]

        # Sum up the King profiles <-- DO WE NEED TO WEIGHT OR NORMALIZE THESE IN SOME WAY?
        PSFfit += PSFmodE

    return PSFfit, Rwing

#----------------------------------------------------------------------------------------------------
# King function
def king(R, Rc, beta, A=1.0):
    '''
    R    - Radius [pixels]
    Rc   - Core radius [pixels]
    beta - Slope [-]
    A    - Amplitude, or normalization factor
    '''
    PSF = A * (1.0 + (R / Rc)**2.0)**(-1.0 * beta)
    return PSF

#----------------------------------------------------------------------------------------------------
# Chi-square and degrees of freedom
def chisqdof(ydata, ymod, Npars):
    chisq = np.sum((ydata - ymod)**2 / ymod)
    dof   = ydata.size - 1 - Npars
    return chisq, dof

#----------------------------------------------------------------------------------------------------
# Return the core radius [pixels] and slope [-] of the model PSF as functions of an input energy and off-axis angle
def plane_function(E, theta, ffits='/Users/salvesen/research/fcolabs/caldb/swxpsf20010101v005.fits'):
    '''
    E     - Energy [keV]
    theta - Off-axis position [degrees]
    ffits - FITS file containing the XRT WT mode PSF coefficients (see Moretti et al. 2005)
    '''
    # Collect the values of the plane function coefficients
    hdul   = fits.open(ffits)
    coeffs = hdul['WT_PSF_COEF'].data
    hdul.close()
    a2, a3 = coeffs[2][1], coeffs[3][1]
    b2, b3 = coeffs[2][2], coeffs[3][2]
    c2, c3 = coeffs[2][3], coeffs[3][3]
    d2, d3 = coeffs[2][4], coeffs[3][4]

    # Unit conversions: E [keV / 0.1], theta [arcmin / 0.1]
    E     = E / 0.1             # [keV / 0.1]
    theta = theta * 60.0 / 0.1  # [arcmin / 0.1]
    
    # Calculate the core radius "Rc" and the slope "beta" for the input energy and position
    Rc   = a2 + (b2 * theta) + (c2 * E) + (d2 * E * theta)  # Rc(E, theta)   [pixels] <-- confirmed w/ Swift Help Desk
    beta = a3 + (b3 * theta) + (c3 * E) + (d3 * E * theta)  # beta(E, theta) [-]

    return Rc, beta

#----------------------------------------------------------------------------------------------------
# Determine the background-subtracted count rate for each region in a list of spectral files
def bksub_count_rate_regions(fsrc_list, fbkg_list, phacuts=None, stdev=True):
    '''
    fsrc_list - Source spectrum file list
    fbkg_list - Background spectrum file list
    phacuts   - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    stdev     - True/False flag to return the standard deviation in the count rate
    '''
    Nreg    = len(fsrc_list)
    cntrate = np.zeros(Nreg)
    stdrate = np.zeros(Nreg)
    for i in np.arange(Nreg):
        fsrc = fsrc_list[i]
        fbkg = fbkg_list[i]
        cntrate[i], stdrate[i] = bkgsub_count_rate(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts, stdev=stdev)
    return cntrate, stdrate

#----------------------------------------------------------------------------------------------------
# Write the radial intensity profiles for each pair of source/background files (all [Rin,Rout] values) to an HDF5 file
def write_radial_intensity(fout, fsrc_list, fbkg_list, Rin, Rout, PSFmin=10.0, phacuts=None, dE=None):
    '''
    fout      - HDF5 output filename containing the results of the radial intensity profile analysis
    fsrc_list - Source spectrum file list
    fbkg_list - Background spectrum file list
    Rin       - Array of annulis inner radii corresponding to the source/background lists [pixels]
    Rout      - Array of annulus outer radii corresponding to the source/background lists [pixels]
    phacuts   - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    '''
    # Collect the PSF radial profile [counts/sec]...
    # ...setup the energy grid [keV] for the PSF model fitting and get the off-axis angle [degrees]
    Rdata, PSFdata, PSFerr, PSFdataE, PSFerrE, Egrid, theta = \
        radial_intensity(fsrc_list=fsrc_list, fbkg_list=fbkg_list, Rin=Rin, Rout=Rout, phacuts=phacuts, dE=dE)

    # Fit the radial intensity profile with the in-flight calibrated PSF model
    PSFfit, Rwing = fitKing2PSF(Rdata=Rdata, PSFdataE=PSFdataE, PSFerrE=PSFerrE, Egrid=Egrid, theta=theta, PSFmin=PSFmin)

    # Write out the results to an HDF5 file
    f = h5py.File(fout, 'w')
    f.create_dataset('Rin',     data=Rin)      # [pixels]
    f.create_dataset('Rout',    data=Rout)     # [pixels]
    f.create_dataset('Rdata',   data=Rdata)    # [pixels]
    f.create_dataset('Rwing',   data=Rwing)    # [pixels]
    f.create_dataset('PSFmin',  data=PSFmin)   # [counts/second]
    f.create_dataset('PSFdata', data=PSFdata)  # [counts/second]
    f.create_dataset('PSFerr',  data=PSFerr)   # [counts/second]
    f.create_dataset('PSFfit',  data=PSFfit)   # [counts/second]
    f.close()


#----------------------------------------------------------------------------------------------------
