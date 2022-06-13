import numpy as np
import pylab as plt
import corner
import gs_stats
import h5py
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter, MaxNLocator, NullFormatter
from matplotlib.colors import LogNorm

'''
PURPOSE:
--------

'''
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

def HalfGauss1D_pdf(x, mu=0.0, sigma=1.0, A=None):
    if (A is None):
        A = np.sqrt(2.0) / np.sqrt(np.pi * sigma**2)
    pdf = A * np.exp(-1.0 * (x - mu)**2 / (2.0 * sigma**2))
    return pdf

#====================================================================================================
# LOAD IN THE DATA

# Input filenames
tag   = "iFlat_"
obsID = "00030009021"
rdir  = "/Users/salvesen/research/fcolabs/results/" + obsID + "/"
fburn = rdir + tag + "chain_burn_" + obsID + ".hdf5"
fmcmc = rdir + tag + "chain_mcmc_" + obsID + ".hdf5"

# Output filenames
dout       = "/Users/salvesen/research/fcolabs/results/" + obsID + "/"
tag_walk   = tag + "Walker_Paths_"
tag_pdf    = tag + "Parameter_PDF_"
tag_2pdf   = tag + "Parameter_2PDF_"
tag_mcmc   = tag + "MCMC_"
fwalk_nH   = dout + tag_walk + "nH.png"
fwalk_a    = dout + tag_walk + "a.png"
fwalk_i    = dout + tag_walk + "i.png"
fwalk_Mbh  = dout + tag_walk + "Mbh.png"
fwalk_Mdd  = dout + tag_walk + "Mdd.png"
fwalk_Dbh  = dout + tag_walk + "Dbh.png"
fwalk_hd   = dout + tag_walk + "hd.png"
fwalk_norm = dout + tag_walk + "norm.png"
fwalk_prob = dout + tag_walk + "prob.png"
fpdf_nH    = dout + tag_pdf  + "nH.png"
fpdf_a     = dout + tag_pdf  + "a.png"
fpdf_i     = dout + tag_pdf  + "i.png"
fpdf_Mbh   = dout + tag_pdf  + "Mbh.png"
fpdf_Mdd   = dout + tag_pdf  + "Mdd.png"
fpdf_Dbh   = dout + tag_pdf  + "Dbh.png"
fpdf_hd    = dout + tag_pdf  + "hd.png"
fpdf_norm  = dout + tag_pdf  + "norm.png"
fcorner    = dout + tag_mcmc + "corner.png"
f2pdf_a_hd = dout + tag_2pdf + "a_hd.png"

# Range of each parameter <-- GET THE MIN/MAX FROM THE PARAMETER FILE, BUT HARD CODE THEM IN FOR NOW
# ALSO GET THE INITIAL PARAMETERS FROM CHI^2 FITTING, OVERPLOT THESE, TOO!
nH_min,   nH_max   = 0.0, 10.0
a_min,    a_max    = -1.0, 1.0#0.9999
i_min,    i_max    = 60.0, 90.0#0.0, 85.0
Mbh_min,  Mbh_max  = 0.0, 10.0
Mdd_min,  Mdd_max  = 0.0, 50.0
Dbh_min,  Dbh_max  = 0.0, 10.0
hd_min,   hd_max   = 1.0, 4.0#1.0, 10.0
norm_min, norm_max = 0.0, 10.0
range_nH   = (nH_min, nH_max)
range_a    = (a_min, a_max)
range_i    = (i_min, i_max)
range_Mbh  = (Mbh_min, Mbh_max)
range_Mdd  = (Mdd_min, Mdd_max)
range_Dbh  = (Dbh_min, Dbh_max)
range_hd   = (hd_min, hd_max)
range_norm = (norm_min, norm_max)
range_list = [range_nH, range_a, range_i, range_Mbh, range_Mdd, range_Dbh, range_hd, range_norm]

# PRIORS
# ------
# Mean / Sigma <-- GET THESE SOMEHOW INSTEAD OF HARDCODING
mean_i,   sigma_i   = 85.0,  2.0#68.65, 0.91   # [degrees]  (Alternative: 68.65, 0.91)
mean_Mbh, sigma_Mbh = 5.40,  0.18  # [Msun]
mean_Dbh, sigma_Dbh = 3.20,  0.2   # [kpc]
# Calculate the normalized priors for M, D, i
Nx        = 1000
x_i       = np.linspace(i_min,   85.0,   Nx)
x_Mbh     = np.linspace(Mbh_min, Mbh_max, Nx)
x_Dbh     = np.linspace(Dbh_min, Dbh_max, Nx)
pdf_i     = Gauss1D_pdf(x_i,   mu=mean_i,   sigma=sigma_i,   A=None)
#pdf_i     = HalfGauss1D_pdf(x_i,   mu=mean_i,   sigma=sigma_i,   A=None)
pdf_Mbh   = Gauss1D_pdf(x_Mbh, mu=mean_Mbh, sigma=sigma_Mbh, A=None)
pdf_Dbh   = Gauss1D_pdf(x_Dbh, mu=mean_Dbh, sigma=sigma_Dbh, A=None)
prior_i   = (x_i,   pdf_i)
prior_Mbh = (x_Mbh, pdf_Mbh)
prior_Dbh = (x_Dbh, pdf_Dbh)

# Read in the EMCEE chain from the burn-in
f       = h5py.File(fmcmc, 'r')
chains  = f['chain'][:,:,:]  # [walker, iteration, parameter]
lnprobs = f['lnprob'][:,:]   # [walker, iteration]
f.close()

# Collect chains for individual parameters
chains_nH   = chains[:,:,0]
chains_a    = chains[:,:,1]
chains_i    = chains[:,:,2]
chains_Mbh  = chains[:,:,3]
chains_Mdd  = chains[:,:,4]
chains_Dbh  = chains[:,:,5]
chains_hd   = chains[:,:,6]
chains_norm = chains[:,:,7]

# Number of walkers, number of steps, array of steps
Nwalkers = chains.shape[0]
Nsteps   = chains.shape[1]
Npars    = chains.shape[2]
steps    = np.arange(Nsteps) + 1

# Axes labels for each parameter
label_nH   = r"$n_{\rm H}\ \left[ 10^{22} / {\rm cm^{2}} \right]$"
label_a    = r"${\rm Black\ Hole\ Spin\ Parameter},\ a$"
label_i    = r"${\rm Inner\ Disk\ Inclination},\ i_{\rm disk}\ \left[ ^{\circ} \right]$"
label_Mbh  = r"$M\ \left[ M_{\odot} \right]$"
label_Mdd  = r"$\dot{M}\ \left[ {\rm 10^{18}~g/s} \right]$"
label_Dbh  = r"$D\ \left[ {\rm kpc} \right]$"
label_hd   = r"${\rm Color\ Correction\ Factor},\ f_{\rm col}$"
label_norm = r"$K$"
label_prob = r"$\ln \left( p \right)$"
label_list = [label_nH, label_a, label_i, label_Mbh, label_Mdd, label_Dbh, label_hd, label_norm]

#====================================================================================================

# Round to specified number of sigfigs.
def round_sigfigs(num, sig_figs):
    if num != 0:
        return round(num, -int(np.floor(np.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0

# https://pypi.org/project/corner/1.0.1/
def Gauss2D_sigma_levels(Nsigma):
    cdf = 1 - np.exp(-0.5 * (Nsigma)**2)
    return cdf


#====================================================================================================
# PLOT: Walker Paths

#----------------------------------------------------------------------------------------------------
# Plot the path taken by each walker
def plot_walker_paths(chains, fwalk, ylabel=''):

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"${\rm Iteration\ Step}$")
    ax.set_ylabel(ylabel)

    for i in np.arange(Nwalkers):
        ax.plot(steps, chains[i,:], color='black', linewidth=1, alpha=0.1)

    fig.savefig(fwalk)
    plt.close()

#----------------------------------------------------------------------------------------------------

'''
# Plot the walker paths for all parameters
plot_walker_paths(chains=chains_nH,   fwalk=fwalk_nH,   ylabel=label_nH)
plot_walker_paths(chains=chains_a,    fwalk=fwalk_a,    ylabel=label_a)
plot_walker_paths(chains=chains_i,    fwalk=fwalk_i,    ylabel=label_i)
plot_walker_paths(chains=chains_Mbh,  fwalk=fwalk_Mbh,  ylabel=label_Mbh)
plot_walker_paths(chains=chains_Mdd,  fwalk=fwalk_Mdd,  ylabel=label_Mdd)
plot_walker_paths(chains=chains_Dbh,  fwalk=fwalk_Dbh,  ylabel=label_Dbh)
plot_walker_paths(chains=chains_hd,   fwalk=fwalk_hd,   ylabel=label_hd)
plot_walker_paths(chains=chains_norm, fwalk=fwalk_norm, ylabel=label_norm)
plot_walker_paths(chains=lnprobs,     fwalk=fwalk_prob, ylabel=label_prob)
'''

#====================================================================================================
# PLOT: Parameter 1D PDFs

#----------------------------------------------------------------------------------------------------
def plot_parameter_PDF(chains, fpdf, Ncut=0, Nbins=100, xlabel='', prior=None, range=None, peakErr=False):
    '''
    Inputs:
    -------
    chains
    fpdf
    Ncut
    Nbins
    
    Overplot the prior
    '''
    # From the chains, unpack the path taken by each individual walker
    walkers = chains[:,Ncut:].flatten()

    # Normalized probability density distribution of the walker paths
    xhist, xedges = np.histogram(walkers, bins=Nbins, range=range, density=True)
    xcenters = (np.roll(xedges,-1)[0:-1] + xedges[0:-1]) * 0.5

    # Plotting preferences
    dpi   = 300
    xsize = 8.4
    ysize = 8.4
    fs, ls, lw, pad, tlmaj, tlmin = 32, 24, 2, 12, 10, 5
    left, bottom, xbox, ybox = 0.165, 0.15, 0.80, 0.80
    
    plt.rcParams['axes.linewidth'] = lw
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    ax = plt.axes([left, bottom, xbox, ybox])
    ax.set_xlabel(xlabel, fontsize=fs, labelpad=pad)
    ax.set_ylabel(r"${\rm Posterior\ PDF}$", fontsize=fs, labelpad=pad)
    
    ax.tick_params('both', direction='in', which='major', labelsize=ls, length=tlmaj, width=lw, pad=pad)
    ax.tick_params('x', direction='in', which='minor', length=tlmin, width=lw)
    ax.tick_params('y', direction='in', which='minor', length=tlmin, width=lw)

    # Plot the PDF
    ax.plot(xcenters, xhist, color='black', linestyle='solid', drawstyle='steps-mid', linewidth=lw)
    
    # Plot the peak and +/- uncertainties
    if (peakErr is True):
        xpeak, xm68, xp68 = gs_stats.error1D(x=xcenters, y=xhist)
        ax.axvline(x=xpeak,      ymin=0.0, ymax=1.0, color='black', linestyle='dashed', linewidth=lw)
        ax.axvline(x=xpeak-xm68, ymin=0.0, ymax=1.0, color='black', linestyle='dotted', linewidth=lw)
        ax.axvline(x=xpeak+xp68, ymin=0.0, ymax=1.0, color='black', linestyle='dotted', linewidth=lw)
        xpeak_str = str(round_sigfigs(xpeak, 3))
        xm68_str  = str(round_sigfigs(xm68, 1))
        xp68_str  = str(round_sigfigs(xp68, 1))
        #x_str     = label_list[i] + r"$\ = " + xpeak_str + "^{+" + xp68_str + "}_{-" + xm68_str + "}$"
        #plt.title(x_str, fontsize=ls)#, pad=pad)

    # Plot the prior
    if (prior is not None):
        ax.plot(prior[0], prior[1], color='blue', linestyle='dashed', linewidth=lw)

    xmin = range[0]
    xmax = range[1]
    ax.set_xlim(xmin, xmax)  # Set x-limits (must come after plotting)

    fig.savefig(fpdf, bbox_inches=0, dpi=dpi)
    plt.close()

#----------------------------------------------------------------------------------------------------

fpdf_i = "/Users/salvesen/Desktop/test_i.png"

# Plot the 1D PDFs for all parameters
Nbins = 50
Ncut  = 0
#plot_parameter_PDF(chains=chains_nH,   fpdf=fpdf_nH,   Nbins=Nbins, Ncut=Ncut, range=range_nH,   xlabel=label_nH)
#plot_parameter_PDF(chains=chains_a,    fpdf=fpdf_a,    Nbins=Nbins, Ncut=Ncut, range=range_a,    xlabel=label_a)#, peakErr=True)
plot_parameter_PDF(chains=chains_i,    fpdf=fpdf_i,    Nbins=Nbins, Ncut=Ncut, range=range_i,    xlabel=label_i,   prior=prior_i)
#plot_parameter_PDF(chains=chains_Mbh,  fpdf=fpdf_Mbh,  Nbins=Nbins, Ncut=Ncut, range=range_Mbh,  xlabel=label_Mbh, prior=prior_Mbh)
#plot_parameter_PDF(chains=chains_Mdd,  fpdf=fpdf_Mdd,  Nbins=Nbins, Ncut=Ncut, range=range_Mdd,  xlabel=label_Mdd)
#plot_parameter_PDF(chains=chains_Dbh,  fpdf=fpdf_Dbh,  Nbins=Nbins, Ncut=Ncut, range=range_Dbh,  xlabel=label_Dbh, prior=prior_Dbh)
#plot_parameter_PDF(chains=chains_hd,   fpdf=fpdf_hd,   Nbins=Nbins, Ncut=Ncut, range=range_hd,   xlabel=label_hd, peakErr=True)
#plot_parameter_PDF(chains=chains_norm, fpdf=fpdf_norm, Nbins=Nbins, Ncut=Ncut, range=range_norm, xlabel=label_norm)

quit()
#====================================================================================================
# PLOT: Corner plot for all parameters

#----------------------------------------------------------------------------------------------------
def plot_2Dpdf(fout, xchains, ychains, xlabel='', ylabel='', xrange=None, yrange=None, Nbins=100):

    x = xchains.flatten()
    y = ychains.flatten()
    
    # Plotting preferences
    dpi   = 300
    xsize = 8.4
    ysize = 8.4
    fs, ls, lw, pad, tlmaj, tlmin = 32, 24, 2, 12, 10, 5
    left, bottom, xbox, ybox = 0.165, 0.15, 0.80, 0.80
    
    plt.rcParams['axes.linewidth'] = lw
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    ax = plt.axes([left, bottom, xbox, ybox])
    ax.set_xlabel(xlabel, fontsize=fs, labelpad=pad)
    ax.set_ylabel(ylabel, fontsize=fs, labelpad=pad)
    
    ax.tick_params('both', direction='in', which='major', labelsize=ls, length=tlmaj, width=lw, pad=pad)
    ax.tick_params('x', direction='in', which='minor', length=tlmin, width=lw)
    ax.tick_params('y', direction='in', which='minor', length=tlmin, width=lw)

    xyhist, xedges, yedges, image = plt.hist2d(x, y, bins=Nbins, normed=True, cmap='Greys')
    #plt.contour(xyhist.transpose(), extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], norm=True, linewidths=lw, colors=['red','green','blue'], levels=levels, origin='lower')
    xmin, xmax = xrange[0], xrange[1]
    ymin, ymax = yrange[0], yrange[1]
    ax.set_xlim(xmin, xmax)  # Set x-limits (must come after plotting)
    ax.set_ylim(ymin, ymax)  # Set y-limits (must come after plotting)

    # Save the figure
    fig.savefig(fout, bbox_inches=0, dpi=dpi)
    plt.close()

# Plot the fcol-a parameter degeneracies
plot_2Dpdf(f2pdf_a_hd, xchains=chains_a, ychains=chains_hd, xlabel=label_a, ylabel=label_hd, xrange=range_a, yrange=range_hd)

#----------------------------------------------------------------------------------------------------
# STILL NEED TO OVERPLOT THE PRIORS
# LEVELS ARE NOT RIGHT https://people.ucsc.edu/~ianc/python/kdestats.html
# TILT THE TEXT ON THE AXES TO SQUEEZE IT IN BETTER
def plot_corner(fout, Nbins=50):

    # Flatter the "chains" 3D array [walker, iteration, parameter] to a 2D array [all iterations, parameter]
    #chains2D = chains.reshape(-1, chains.shape[-1])  # [iteration, parameter]

    # Plotting preferences
    dpi   = 300
    xsize = 8.4
    ysize = 8.4
    fs, ls, lw, pad, tlmaj, tlmin = 12, 10, 1, 5, 2, 1
    nullfmt = NullFormatter()

    sigmas = [1.0, 2.0, 3.0]
    levels = []
    for i in np.arange(len(sigmas)):
        levels.append(Gauss2D_sigma_levels(Nsigma=sigmas[i]))

    # FUUUUUUCK
    plt.rcParams['axes.linewidth'] = lw
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    Lgap, Rgap, Bgap, Tgap, Igap = 0.1, 0.05, 0.1, 0.05, 0.01  # Gaps left, right, bottom, top, in-between
    xbox = (1.0 - (Lgap + Rgap) - (Npars - 1) * Igap) / Npars
    ybox = (1.0 - (Bgap + Tgap) - (Npars - 1) * Igap) / Npars

    # Loop through rows
    bottom = Bgap
    for i in reversed(np.arange(Npars)):

        # Loop through columns
        left = Lgap
        for j in np.arange(i+1):
            
            # Define the axes object
            ax = plt.axes([left, bottom, xbox, ybox])

            # Find the min/max parameter ranges for plotting
            xmin, xmax = range_list[j]
            ymin, ymax = range_list[i]

            # Set tick properties
            ax.tick_params('both', direction='in', which='major', labelsize=ls, length=tlmaj, width=lw, pad=pad)
            ax.tick_params('x', direction='in', which='minor', length=tlmin, width=lw)
            ax.tick_params('y', direction='in', which='minor', length=tlmin, width=lw)

            # Label the x-axis
            if (i == (Npars-1)): ax.set_xlabel(label_list[j], fontsize=fs, labelpad=pad)
            else: ax.xaxis.set_major_formatter(nullfmt)  # No labels on the x-axis
            
            # Label the y-axis
            if ((j == 0) and (i != 0)): ax.set_ylabel(label_list[i], fontsize=fs, labelpad=pad)
            else: ax.yaxis.set_major_formatter(nullfmt)  # No labels on the x-axis

            # Plot the 2D contour
            if (j != i):
                x = chains[:,:,j].flatten()
                y = chains[:,:,i].flatten()
                xyhist, xedges, yedges, image = plt.hist2d(x, y, bins=Nbins, normed=True, cmap='Greys')
                #plt.contour(xyhist.transpose(), extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], norm=True, linewidths=lw, colors=['red','green','blue'], levels=levels, origin='lower')
                ax.set_xlim(xmin, xmax)  # Set x-limits (must come after plotting)
                ax.set_ylim(ymin, ymax)  # Set y-limits (must come after plotting)

            # Plot the 1D PDF
            if (j == i):
                x = chains[:,:,j].flatten()
                xhist, xedges = np.histogram(x, bins=Nbins, range=range_list[j], density=True)
                xcenters = (np.roll(xedges,-1)[0:-1] + xedges[0:-1]) * 0.5
                ax.plot(xcenters, xhist, color='black', linestyle='solid', drawstyle='steps-mid', linewidth=lw)
                xpeak, xm68, xp68 = gs_stats.error1D(x=xcenters, y=xhist)
                ax.axvline(x=xpeak,      ymin=0.0, ymax=1.0, color='black', linestyle='dashed', linewidth=lw)
                ax.axvline(x=xpeak-xm68, ymin=0.0, ymax=1.0, color='black', linestyle='dotted', linewidth=lw)
                ax.axvline(x=xpeak+xp68, ymin=0.0, ymax=1.0, color='black', linestyle='dotted', linewidth=lw)
                xpeak_str = str(round_sigfigs(xpeak, 3))
                xm68_str  = str(round_sigfigs(xm68, 1))
                xp68_str  = str(round_sigfigs(xp68, 1))
                x_str     = label_list[i] + r"$\ = " + xpeak_str + "^{+" + xp68_str + "}_{-" + xm68_str + "}$"
                plt.title(x_str, fontsize=ls)#, pad=pad)
                ax.set_xlim(xmin, xmax)  # Set x-limits (must come after plotting)

            # Increment the left location
            left = left + (xbox + Igap)

        # Increment the bottom location
        bottom = bottom + (ybox + Igap)

    # Save the figure
    fig.savefig(fout, bbox_inches=0, dpi=dpi)
    plt.close()

# Plot the corner plot
nH_min,   nH_max   = 0.9, 1.1
a_min,    a_max    = -1.0, 1.0
i_min,    i_max    = 65.0, 85.0
Mbh_min,  Mbh_max  = 4.8, 6.0
Mdd_min,  Mdd_max  = 1.0, 10.0
Dbh_min,  Dbh_max  = 2.5, 4.0
hd_min,   hd_max   = 1.0, 3.0
norm_min, norm_max = 0.5, 1.5
range_nH   = (nH_min, nH_max)
range_a    = (a_min, a_max)
range_i    = (i_min, i_max)
range_Mbh  = (Mbh_min, Mbh_max)
range_Mdd  = (Mdd_min, Mdd_max)
range_Dbh  = (Dbh_min, Dbh_max)
range_hd   = (hd_min, hd_max)
range_norm = (norm_min, norm_max)
range_list = [range_nH, range_a, range_i, range_Mbh, range_Mdd, range_Dbh, range_hd, range_norm]
label_nH   = r"$n_{\rm H}$"
label_a    = r"$a$"
label_i    = r"$i_{\rm disk}$"
label_Mbh  = r"$M$"
label_Mdd  = r"$\dot{M}$"
label_Dbh  = r"$D$"
label_hd   = r"$f_{\rm col}$"
label_norm = r"$K$"
label_prob = r"$\ln \left( p \right)$"
label_list = [label_nH, label_a, label_i, label_Mbh, label_Mdd, label_Dbh, label_hd, label_norm]
plot_corner(fout=fcorner)
