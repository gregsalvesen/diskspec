import numpy as np
import pylab as plt
import corner
import gs_stats
from astropy.io import fits

nwalk = 10000
pre   = "Nwalk_" + str(nwalk) + "_"

# Output filenames
dout       = "/Users/salvesen/research/fcolabs/results/"
fwalk_a    = dout + pre + "xspec_Walker_Paths_a.png"
fwalk_i    = dout + pre + "xspec_Walker_Paths_i.png"
fwalk_Mbh  = dout + pre + "xspec_Walker_Paths_Mbh.png"
fwalk_Mdd  = dout + pre + "xspec_Walker_Paths_Mdd.png"
fwalk_Dbh  = dout + pre + "xspec_Walker_Paths_Dbh.png"
fwalk_hd   = dout + pre + "xspec_Walker_Paths_hd.png"
fwalk_stat = dout + pre + "xspec_Walker_Paths_stat.png"
fpdf_a     = dout + pre + "xspec_Parameter_PDF_a.png"
fpdf_i     = dout + pre + "xspec_Parameter_PDF_i.png"
fpdf_Mbh   = dout + pre + "xspec_Parameter_PDF_Mbh.png"
fpdf_Mdd   = dout + pre + "xspec_Parameter_PDF_Mdd.png"
fpdf_Dbh   = dout + pre + "xspec_Parameter_PDF_Dbh.png"
fpdf_hd    = dout + pre + "xspec_Parameter_PDF_hd.png"

# Read in the MCMC chain output FITS file
rdir  = "/Users/salvesen/research/fcolabs/results/"
ffits = rdir + pre + "xspec_chain_burn_kerrbb.fits"
hdul  = fits.open(ffits)
hdr0  = hdul['primary'].header
hdr1  = hdul['chain'].header
data0 = hdul['primary'].data
data1 = hdul['chain'].data
hdul.close()

# Collect the chains
chains_a    = data1['a__3']
chains_i    = data1['i__4']
chains_Mbh  = data1['Mbh__5']
chains_Mdd  = data1['Mdd__6']
chains_Dbh  = data1['Dbh__7']
chains_hd   = data1['hd__8']
chains_stat = data1['FIT_STATISTIC']

# Number of walkers, number of steps, array of steps
Nwalkers    = hdr1['nwalkers']
Nsteps      = int(len(chains_stat) / Nwalkers)
steps       = np.arange(Nsteps) + 1

# Range of each parameter <-- GET THE MIN/MAX FROM THE PARAMETER FILE, BUT HARD CODE THEM IN FOR NOW
# ALSO GET THE INITIAL PARAMETERS FROM CHI^2 FITTING, OVERPLOT THESE, TOO!
a_min,   a_max   = 0.9, 0.9999
i_min,   i_max   = 0.0, 85.0
Mbh_min, Mbh_max = 0.0, 50.0
Mdd_min, Mdd_max = 0.0, 10.0
Dbh_min, Dbh_max = 0.0, 50.0
hd_min,  hd_max  = 1.0, 10.0
range_a   = (a_min, a_max)
range_i   = (i_min, i_max)
range_Mbh = (Mbh_min, Mbh_max)
range_Mdd = (Mdd_min, Mdd_max)
range_Dbh = (Dbh_min, Dbh_max)
range_hd  = (hd_min, hd_max)

# Mean / Sigma <-- GET THESE SOMEHOW INSTEAD OF HARDCODING
mean_i,   sigma_i   = 68.65, 0.91  # [degrees]
mean_Mbh, sigma_Mbh = 5.40,  0.18  # [Msun]
mean_Dbh, sigma_Dbh = 3.20,  0.2   # [kpc]

# Axes labels for each parameter
label_a    = r"$a$"
label_i    = r"$i_{\rm disk}\ \left[ {\rm degrees} \right]$"
label_Mbh  = r"$M\ \left[ M_{\odot} \right]$"
label_Mdd  = r"$\dot{M}\ \left[ {\rm 10^{18} g\ s^{-1}} \right]$"
label_Dbh  = r"$D\ \left[ {\rm kpc} \right]$"
label_hd   = r"$f_{\rm col}$"
label_stat = r"${\rm Fit\ Statistic}$"

#====================================================================================================
# FUNCTIONS

#----------------------------------------------------------------------------------------------------
# Extract the path taken by each walker from the output of XSPEC MCMC
def extract_walker_paths(chains, Nwalkers, Ncut=0):
    '''
    The XSPEC MCMC output chains are 1D arrays of length Nwalkers*Nsteps organized like so:
        Let w1, w2, w3, etc. represent walker # 1, 2, 3, etc.
        Let s1, s2, s3, etc. represent step # 1, 2, 3, etc.
        chains = [w1s1,w2s1,w3s1,...,w1s2,w2s2,w3s2,...,w1s3,w2s3,w3s3,...]

    This routine outputs a list of lists, each containing the path taken by a walker organized like so:
        walkers = [[w1s1,w1s2,w1s3,...], [w2s1,w2s2,w1s3,...], [w3s1,w3s2,w3s3,...], ...]

    "Ncut" is the (inclusive) number of steps to toss out (i.e., indices [0:Ncut] are removed).
    For example, specifying Ncut=10 will return:
        walkers = [[w1s11,w1s12,w1s13,...], [w2s11,w2s12,w1s13,...], [w3s11,w3s12,w3s13,...], ...]
    '''
    walkers = []
    for i in np.arange(Nwalkers):
        this_walker = chains[i::Nwalkers]   # Using [start:stop:step] syntax
        walkers.append(this_walker[Ncut:])  # Toss out the first "Ncut" steps
    return walkers

#====================================================================================================
# PLOT: Walker Paths

#----------------------------------------------------------------------------------------------------
# Plot the path taken by each walker
def plot_walker_paths(chains, fwalk, ylabel=''):

    # From the chains, unpack the path taken by each individual walker
    walkers = extract_walker_paths(chains=chains, Nwalkers=Nwalkers)

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"${\rm Iteration\ Step}$")
    ax.set_ylabel(ylabel)
    for i in np.arange(Nwalkers):
        ax.plot(steps, walkers[i], color='black', linewidth=1, alpha=0.1)
    fig.savefig(fwalk)
    plt.close()
#----------------------------------------------------------------------------------------------------

# Plot the walker paths for all parameters
plot_walker_paths(chains=chains_a,    fwalk=fwalk_a,    ylabel=label_a)
plot_walker_paths(chains=chains_i,    fwalk=fwalk_i,    ylabel=label_i)
plot_walker_paths(chains=chains_Mbh,  fwalk=fwalk_Mbh,  ylabel=label_Mbh)
plot_walker_paths(chains=chains_Mdd,  fwalk=fwalk_Mdd,  ylabel=label_Mdd)
plot_walker_paths(chains=chains_Dbh,  fwalk=fwalk_Dbh,  ylabel=label_Dbh)
plot_walker_paths(chains=chains_hd,   fwalk=fwalk_hd,   ylabel=label_hd)
plot_walker_paths(chains=chains_stat, fwalk=fwalk_stat, ylabel=label_stat)

#====================================================================================================
# PLOT: Parameter 1D PDFs

#----------------------------------------------------------------------------------------------------
def Gauss1D_pdf(x, mu=0.0, sigma=1.0, A=None):
    '''
    Inputs:
    -------
    x     - Value or array/list of bins at which to compute the PDF
    mu    - Mean value (i.e. peak x-value) of the Gaussian distribution
    sigma - Standard deviation from the mean (i.e. 1-sigma value)
    
    Output:
    -------
    Gaussian (normal) probability distribution function (PDF)
    '''
    if (A is None):
        A = 1.0 / np.sqrt(2.0 * np.pi * sigma**2)
    pdf = A * np.exp(-1.0 * (x - mu)**2 / (2.0 * sigma**2))
    return pdf

#----------------------------------------------------------------------------------------------------
def plot_parameter_PDF(chains, fpdf, Ncut=0, Nbins=100, xlabel='', prior=None, range=None):
    '''
    Inputs:
    -------
    
    Overplot the prior
    '''
    # From the chains, unpack the path taken by each individual walker
    walkers = extract_walker_paths(chains=chains, Nwalkers=Nwalkers, Ncut=Ncut)

    # Normalized probability density distribution of the walker paths
    xhist, xedges = np.histogram(walkers, bins=Nbins, range=range, density=True)
    xcenters = (np.roll(xedges,-1)[0:-1] + xedges[0:-1]) * 0.5
    
    # Parameter peak and uncertainties
    #xpeak, xm68, xp68 = gs_stats.error1D(x=xcenters, y=xhist)

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"${\rm PDF}$")
    
    # Plot the PDF
    ax.plot(xcenters, xhist, color='black', linestyle='solid', drawstyle='steps-mid')
    
    # Plot the peak and +/- uncertainties
    #ax.plot(xcenters, xhist, color='black', linestyle='solid')

    # Plot the prior
    if (prior is not None):
        ax.plot(prior[0], prior[1], color='blue', linestyle='dashed')
    
    fig.savefig(fpdf)
    plt.close()

#----------------------------------------------------------------------------------------------------

# Calculate the normalized priors for M, D, i
Nx    = 1000
x_i   = np.linspace(i_min,   i_max,   Nx)
x_Mbh = np.linspace(Mbh_min, Mbh_max, Nx)
x_Dbh = np.linspace(Dbh_min, Dbh_max, Nx)
pdf_i   = Gauss1D_pdf(x_i,   mu=mean_i,   sigma=sigma_i,   A=None)
pdf_Mbh = Gauss1D_pdf(x_Mbh, mu=mean_Mbh, sigma=sigma_Mbh, A=None)
pdf_Dbh = Gauss1D_pdf(x_Dbh, mu=mean_Dbh, sigma=sigma_Dbh, A=None)
prior_i   = (x_i,   pdf_i)
prior_Mbh = (x_Mbh, pdf_Mbh)
prior_Dbh = (x_Dbh, pdf_Dbh)

# Plot the 1D PDFs for all parameters
Nbins = 100
Ncut  = 500
plot_parameter_PDF(chains=chains_a,   fpdf=fpdf_a,   Nbins=Nbins, Ncut=Ncut, range=range_a,   xlabel=label_a)
plot_parameter_PDF(chains=chains_i,   fpdf=fpdf_i,   Nbins=Nbins, Ncut=Ncut, range=range_i,   xlabel=label_i,   prior=prior_i)
plot_parameter_PDF(chains=chains_Mbh, fpdf=fpdf_Mbh, Nbins=Nbins, Ncut=Ncut, range=range_Mbh, xlabel=label_Mbh, prior=prior_Mbh)
plot_parameter_PDF(chains=chains_Mdd, fpdf=fpdf_Mdd, Nbins=Nbins, Ncut=Ncut, range=range_Mdd, xlabel=label_Mdd)
plot_parameter_PDF(chains=chains_Dbh, fpdf=fpdf_Dbh, Nbins=Nbins, Ncut=Ncut, range=range_Dbh, xlabel=label_Dbh, prior=prior_Dbh)
plot_parameter_PDF(chains=chains_hd,  fpdf=fpdf_hd,  Nbins=Nbins, Ncut=Ncut, range=range_hd,  xlabel=label_hd)


#====================================================================================================
# PLOT: Parameter 2D PDFs


#====================================================================================================
# PLOT: Corner plot for all parameters

