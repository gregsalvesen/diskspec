import numpy as np
import h5py
import pylab as plt
from matplotlib import colors, colorbar
from matplotlib.ticker import MultipleLocator, MaxNLocator

'''
PURPOSE:
--------
Produce plots from the pile-up analysis.
'''

#=============================================================================#
# SET THINGS UP

# Results directory for figure output, list of observation IDs, and tags for output filenames
resdir     = '/Users/salvesen/research/fcolabs/results/'
obsID_list = ['00030009021', '00030009022', '00030009023', '00030009025', '00030009026']
Nobs       = len(obsID_list)
tag_grades_ratio = 'pileup_grades_ratio_'
tag_radial_inten = 'pileup_radial_intensity_'
tag_spec_distort = 'pileup_spectral_distortions_'

# Plotting preferences
dpi   = 300
xsize = 8.4  # [inches]
ysize = 8.4  # [inches]
left, right, bottom, top = 0.15, 0.95, 0.125, 0.85
lw = 1.0

# Create a discrete color map
cmap       = plt.cm.nipy_spectral                        # Choose a color map to discretize
cmaplist   = [cmap(i) for i in range(cmap.N)]            # Extract all colors
cmap       = cmap.from_list('Custom', cmaplist, cmap.N)  # Create the new map
vmin, vmax = 0, 800                                      # Min/Max count rate [counts/second]
bounds     = np.linspace(vmin, vmax, 17)                 # Discrete color grid (every 50 cts/s)
norm       = colors.BoundaryNorm(bounds, cmap.N)         # Normalize
cbticks    = np.linspace(vmin, vmax, 9)                  # Colorbar tick marks (every 100 cts/s)

#=============================================================================#
# GRADES RATIO

# Plot the grades ratio as a function of Rcut for a given Swift/XRT observation
def plot_grades_ratio(fh5, fout):

    # Collect the data
    f = h5py.File(fh5, 'r')
    Rcut    = f['Rcut'][:]                     # [pixels]
    frac    = f['grades_ratio'][:]             # [-]
    fstd    = f['grades_ratio_stdev'][:]       # [-]
    cntrate = f['bkgsub_count_rate'][:]        # [counts / second]
    stdrate = f['bkgsub_count_rate_stdev'][:]  # [counts / second]
    f.close()

    # Plot the grades ratio as a function of Rcut
    ymin, ymax = 0.84, 0.95
    fig = plt.figure(figsize=(xsize,ysize), dpi=dpi)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$R_{\rm cut}\ {\rm [pixels]}$")
    ax.set_ylabel(r"${\rm Counts\ Ratio\ (Grade\ 0 / Grades\ 0}$"+r"--"+r"${\rm 2)}$")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-1, np.max(Rcut)+1)
    scat = ax.scatter(Rcut, frac, c=cntrate, cmap=cmap, norm=norm, marker='o', s=75, edgecolors='black', linewidths=1.0, zorder=2)
    ax.errorbar(x=Rcut, y=frac, yerr=fstd, ecolor='black', elinewidth=lw, linestyle='None', zorder=1)

    # Colorbar axis
    ax2 = fig.add_axes([0.15, 0.875, 0.80, 0.033])
    cb  = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal', spacing='proportional', ticks=cbticks, boundaries=bounds)
    ax2.set_xlabel(r"${\rm Count\ Rate\ [counts\ second^{-1}]}$", size=16)
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=16)
    cb.ax.tick_params(which='minor', axis='both', bottom=False, top=False, left=False, right=False)

    # Save the figure
    fig.savefig(fout, bbox_inches=0, dpi=dpi)
    plt.close()

#=============================================================================#
# RADIAL INTENSITY PROFILE <-- JUST A SKELETON RIGHT NOW!!! NEEDS WORK!!!

# Plot the radial intensity profile and the PSF model fit for a given Swift/XRT observation
def plot_radial_intensity(fh5, fout):

    # Collect the data
    f = h5py.File(fh5, 'r')
    Rin     = f['Rin'][:]        # [pixels]
    Rout    = f['Rout'][:]       # [pixels]
    Rdata   = f['Rdata'][:]      # [pixels]
    Rwing   = f['Rwing'].value   # [pixels]
    PSFmin  = f['PSFmin'].value  # [counts/sec]
    PSFdata = f['PSFdata'][:]    # [counts/sec]
    PSFerr  = f['PSFerr'][:]     # [counts/sec]
    PSFfit  = f['PSFfit'][:]     # [counts/sec]
    f.close()
    
    # Index of the "wing"
    iR0 = np.argmin(np.abs(Rdata - Rwing))

    # Plot the radial intensity profile and the (extrapolated) PSF model fit
    fig = plt.figure(figsize=(xsize,ysize), dpi=dpi)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$R\ {\rm [pixels]}$")
    ax.set_ylabel(r"${\rm Count\ Rate\ [counts\ s^{-1}\ pixel^{-1}]}$")
    ax.set_yscale('log')
    scat = ax.scatter(Rdata, PSFdata, marker='o', s=75, edgecolors='black', linewidths=1.0, zorder=2)#, c=cntrate, cmap=cmap, norm=norm)
    ax.errorbar(x=Rdata, y=PSFdata, yerr=PSFerr, ecolor='black', elinewidth=lw, linestyle='None', zorder=1)
    
    # Solid line = fit; Dashed line = extrapolated
    ax.plot(Rdata[0:iR0+1], PSFfit[0:iR0+1], 'k--', linewidth=lw)
    ax.plot(Rdata[iR0:],    PSFfit[iR0:],    'k-',  linewidth=lw)

    # Data/model ratio
    #print PSFdata[iR0:] / PSFfit[iR0:]

    # Save the figure
    fig.savefig(fout, bbox_inches=0, dpi=dpi)
    plt.close()

#=============================================================================#
# SPECTRAL DISTORTIONS

# Plot the diskbb best-fit Tin as a function of Rcut for a given Swift/XRT observation
def plot_stectral_distortion(fh5, fout):

    # Collect the data
    f = h5py.File(fh5, 'r')
    Rcut      = f['Rcut'][:]                     # [pixels]
    nH        = f['nH'][:]                       # [atoms cm^-2]
    Tin       = f['Tin'][:]                      # [keV]
    TinErr    = f['TinErr'][:]                   # [keV]
    normbb    = f['normbb'][:]                   # [...]
    normbbErr = f['normbbErr'][:]                # [...]
    chisq     = f['chisq'][:]                    # [-]
    dof       = f['dof'][:]                      # [-]
    cntrate   = f['bkgsub_count_rate'][:]        # [counts/second]
    stdrate   = f['bkgsub_count_rate_stdev'][:]  # [counts/second]
    f.close()

    # Plot the diskbb best-fit Tin as a function of Rcut
    ymin, ymax = 0, 2
    fig = plt.figure(figsize=(xsize,ysize), dpi=dpi)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$R_{\rm cut}\ {\rm [pixels]}$")
    ax.set_ylabel(r"$T_{\rm in}\ [{\rm keV}]}$")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-1, np.max(Rcut)+1)
    scat = ax.scatter(Rcut, Tin, c=cntrate, cmap=cmap, norm=norm, marker='o', s=75, edgecolors='black', linewidths=1.0, zorder=2)
    ax.errorbar(x=Rcut, y=Tin, yerr=TinErr, ecolor='black', elinewidth=lw, linestyle='None', zorder=1)
    
    # Colorbar axis
    ax2 = fig.add_axes([0.15, 0.875, 0.80, 0.033])
    cb  = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal', spacing='proportional', ticks=cbticks, boundaries=bounds)
    ax2.set_xlabel(r"${\rm Count\ Rate\ [counts\ second^{-1}]}$", size=16)
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=16)
    cb.ax.tick_params(which='minor', axis='both', bottom=False, top=False, left=False, right=False)

    # Save the figure
    fig.savefig(fout, bbox_inches=0, dpi=dpi)
    plt.close()

#=============================================================================#
# MAKE THE PLOTS

print ""
print "PILE-UP PLOTS"

# Loop through each observation ID
for i in np.arange(Nobs):

    # Current observation ID
    obsID = obsID_list[i]
    print "...Making plots for ObsID: ", obsID

    # Grades Ratio
    fh5_grades_ratio  = resdir + obsID + '/' + tag_grades_ratio + obsID + '.h5'
    fout_grades_ratio = resdir + obsID + '/' + tag_grades_ratio + obsID + '.png'
    #plot_grades_ratio(fh5=fh5_grades_ratio, fout=fout_grades_ratio)

    # Radial Intensity Profile
    fh5_radial_inten  = resdir + obsID + '/' + tag_radial_inten + obsID + '.h5'
    fout_radial_inten = resdir + obsID + '/' + tag_radial_inten + obsID + '.png'
    #plot_radial_intensity(fh5=fh5_radial_inten, fout=fout_radial_inten)

    # Spectral Distortions
    fh5_spec_distort  = resdir + obsID + '/' + tag_spec_distort + obsID + '.h5'
    fout_spec_distort = resdir + obsID + '/' + tag_spec_distort + obsID + '.png'
    #plot_spectral_distortion(fh5=fh5_spec_distort, fout=fout_spec_distort)

#=============================================================================#
