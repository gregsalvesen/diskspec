import numpy as np
import pylab as plt
import corner
import h5py
from matplotlib.ticker import NullFormatter, MultipleLocator, MaxNLocator, AutoMinorLocator, FormatStrFormatter

# Output file
fout = "test.png"

# Read in the data
fmcmc = "emcee.hdf5"
f = h5py.File(fmcmc, 'r')
chains  = f['chain']   # (walker, iteration, param) <-- the burn-in is discarded
lnprobs = f['lnprob']  # (walker, iteration) <-- the burn-in is discarded


# Collect individual parameter chains
chain_NH   = chains[:,:,0].flatten()  # 10^22 [atoms cm^-2]
chain_Mdot = chains[:,:,1].flatten()  # 10^18 [g s^-1]
chain_fcol = chains[:,:,2].flatten()  # [-]

print chain_NH
quit()
f.close()

'''
#====================================================================================================
# Plotting preferences
dpi = 600
alpha = 0.25
#fs, fs_sm, lw, pad, tlmaj, tlmin = 32, 20, 2, 10, 10, 5
fs, fs_sm, lw, pad, tlmaj, tlmin = 10, 7.5, 1.0, 5, 5, 2.5
left, right, bottom, top = 0.15, 0.95, 0.175, 0.975
#xsize, ysize = 8.4, 8.4
xsize, ysize = 2.8, 2.8
xinc_maj, xinc_min = 0.1, 0.025
yinc_maj, yinc_min = 0.05, 0.01
xmajorLocator = MultipleLocator(xinc_maj)
xminorLocator = MultipleLocator(xinc_min)
ymajorLocator = MultipleLocator(yinc_maj)
yminorLocator = MultipleLocator(yinc_min)

# Draw contour levels at 0.5, 1.0, 1.5, 2.0 sigma.
levels = 1.0 - np.exp(-0.5 * np.array([0.5,1.0,1.5,2.0])**2)

# Bins for the filled contours (use the (mu,sigma)-grid resolution of 100x100)
bins = [100,100]

# Men
#xmin, xmax = -6.4, -5.9
#ymin, ymax = 0.0, 0.25
#plt.rcParams['axes.linewidth'] = lw
fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
ax = fig.add_subplot(111)
ax.set_xlabel(r"$f_{\rm col}$", fontsize=fs*1.25, labelpad=pad*0.5)
ax.set_ylabel(r"$Mdot$", fontsize=fs*1.25, labelpad=pad*0.5)
#ax.set_xlim(xmin, xmax)
#ax.set_ylim(ymin, ymax)
#ax.xaxis.set_major_locator(xmajorLocator)
#ax.xaxis.set_minor_locator(xminorLocator)
#ax.yaxis.set_major_locator(ymajorLocator)
#ax.yaxis.set_minor_locator(yminorLocator)
ax.tick_params('x', direction='in', labelsize=fs_sm, length=tlmaj, width=lw, which='major', pad=pad)
ax.tick_params('y', direction='in', labelsize=fs_sm, length=tlmaj, width=lw, which='major', pad=pad*0.5)
ax.tick_params('both', direction='in', length=tlmin, width=lw, which='minor')
#plt.xticks(np.arange(xmin, xmax+xinc_maj, xinc_maj))
#plt.yticks(np.arange(ymin, ymax+yinc_maj, yinc_maj))

corner.hist2d(x=chain_fcol, y=chain_Mdot, range=[(xmin,xmax),(ymin,ymax)], ax=ax, bins=bins, levels=levels, plot_datapoints=True, plot_density=True, plot_contours=True, fill_contours=False, contour_kwargs={'linewidths':1.0}, data_kwargs={'alpha':0.01,'markersize':1})

#ax.plot([mu_fit_M, mu_fit_M], [ymin,ymax], 'b-', linewidth=lw)
#ax.plot([xmin,xmax], [sigma_fit_M,sigma_fit_M], 'b-', linewidth=lw)
#ax.text(0.075, 0.925, r"${\rm BLAH!}$", color='b', transform=ax.transAxes, ha='left', va='top', fontsize=fs)

fig.savefig(fout_M, bbox_inches=0, dpi=dpi)
plt.close()
'''
