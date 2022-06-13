import numpy as np
import pylab as plt
import h5py
from matplotlib.ticker import NullFormatter, FuncFormatter

'''
NAME:
-----
plot_spectrum

PURPOSE:
--------
Plot the spectrum that is saved to an HDF5 file after running xspec_modelfit.py

EXAMPLE:
--------
>> from plot_spectrum import *
>> fspec = "/Users/salvesen/research/fcolabs/results/junk/TBabs_diskbb.hdf5"
>> fout  = "/Users/salvesen/research/fcolabs/results/junk/TBabs_diskbb_spec.png"
>> inst  = plot_spectrum(fspec=fspec)
>> inst.plot_bestfit(fout=fout)

NOTES:
------
- Some bullshit having to do with matplotlib and tick marks:
  https://github.com/matplotlib/matplotlib/issues/8436
  plt.rc('text', usetex=True)
- Need to work on plot_model()
'''
#====================================================================================================
class plot_spectrum:

    #----------------------------------------------------------------------------------------------------
    # Initialization
    def __init__(self, fspec, E_min=0.5, E_max=10.0):
        self.E_min = E_min
        self.E_max = E_max

        # Read in the data, bestfit, and model spectra
        f = h5py.File(fspec, 'r')
        self.xdata     = f['Energy'][:]
        self.ydata     = f['Data'][:]
        self.xdataErr  = f['EnergyErr'][:]
        self.ydataErr  = f['DataErr'][:]
        self.ybfit     = f['BestFit'][:]
        self.yresid    = f['Resid'][:]
        self.yresidErr = f['ResidErr'][:]
        self.yratio    = f['Ratio'][:]
        self.yratioErr = f['RatioErr'][:]
        self.xmodel    = f['Emodel'][:]
        self.ymodel    = f['FEmodel'][:]
        self.xymodel   = f['EFEmodel'][:]
        f.close()

    #----------------------------------------------------------------------------------------------------
    # HELPER FUNCTIONS

    def myFormatter(self, y, pos):
        formatstring = '{:g}'.format(y)
        return formatstring.format(y)

    #----------------------------------------------------------------------------------------------------
    # PLOT: SPECTRUM, BESTFIT MODEL, RADIO, RESIDUALS
    def plot_bestfit(self, fout):

        # Settings
        lwthin    = 0.5
        lwnormal  = 1.0
        lwthick   = 2.0
        padLeft   = 0.15
        padRight  = 0.05
        padBot    = 0.125
        padTop    = 0.05
        divide    = 0.5
        width     = 1.0 - padLeft - padRight
        heightTop = 1.0 - divide - padTop
        heightMid = (divide - padBot) * 0.5
        heightBot = heightMid
        rectTop   = [padLeft, divide, width, heightTop]
        rectMid   = [padLeft, padBot+heightMid, width, heightMid]
        rectBot   = [padLeft, padBot, width, heightBot]
        
        # Define the figure and subplot axes objects
        fig = plt.figure()#figsize=(xsize, ysize), dpi=dpi)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        axTop = plt.axes(rectTop)
        axMid = plt.axes(rectMid)
        axBot = plt.axes(rectBot)

        # Plot the data and best-fit model
        axTop.set_ylabel(r"${\rm normalized\ counts\ s^{-1}\ keV^{-1}}$")
        axTop.set_xlim(self.E_min, self.E_max)
        axTop.set_xscale('log')
        axTop.set_yscale('log')
        axTop.get_yaxis().set_major_formatter(FuncFormatter(self.myFormatter))
        axTop.errorbar(x=self.xdata, y=self.ydata, xerr=self.xdataErr, yerr=self.ydataErr, \
                       fmt='none', color='k', linewidth=lwthin, zorder=0)
        axTop.plot(self.xdata, self.ybfit, color='C1', drawstyle='steps-mid', zorder=1)

        # Plot the data/model ratio
        axMid.set_ylabel(r"${\rm Ratio}$")
        axMid.set_xlim(self.E_min, self.E_max)
        axMid.set_xscale('log')
        axMid.get_yaxis().set_major_formatter(FuncFormatter(self.myFormatter))
        axMid.errorbar(x=self.xdata, y=self.yratio, xerr=self.xdataErr, yerr=self.yratioErr, \
                       fmt='none', color='k', linewidth=lwthin, zorder=0)
        axMid.plot([self.E_min,self.E_max], [1,1], color='C0', linewidth=1.0, zorder=1)

        # Plot the data - model residuals
        axBot.set_xlabel(r"${\rm Energy\ (keV)}$")
        axBot.set_ylabel(r"${\rm Residuals}$")
        axBot.set_xlim(self.E_min, self.E_max)
        axBot.set_xscale('log')
        axBot.get_yaxis().set_major_formatter(FuncFormatter(self.myFormatter))
        axBot.errorbar(x=self.xdata, y=self.yresid, xerr=self.xdataErr, yerr=self.yresidErr, \
                       fmt='none', color='k', linewidth=lwthin, zorder=0)
        axBot.plot([self.E_min,self.E_max], [0,0], color='C0', linewidth=lwnormal, zorder=1)

        # No labels on the x-axis <-- This must come AFTER set_xscale('log') calls, which is fucking stupid
        nullfmt = NullFormatter()  # No labels
        axTop.xaxis.set_major_formatter(nullfmt)
        axTop.xaxis.set_minor_formatter(nullfmt)
        axMid.xaxis.set_major_formatter(nullfmt)
        axMid.xaxis.set_minor_formatter(nullfmt)
    
        # Custom tick labels <-- Don't get me started on how fucking stupid this is
        axBot.set_xticks([1, 10], minor=False)
        axBot.set_xticklabels(['1', '10'], minor=False)
        axBot.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 2, 3, 4, 5, 6, 7, 8, 9], minor=True)
        axBot.set_xticklabels(['0.5', '', '', '', '', '2', '', '', '5', '', '', '', ''], minor=True)

        # Save the figure
        fig.savefig(fout)
        plt.close()

    #----------------------------------------------------------------------------------------------------
    # PLOT: MODEL (???ENERGIES EXTEND???UNFOLDED???)
    def plot_model(self, fout):
    
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.set_xlabel(r"${\rm Energy\ (keV)}$")
        ax.set_ylabel(r"${\rm erg\ s^{-1}\ cm^{-2}\ keV^{-1}}$")
        ax.set_xlim(self.E_min, self.E_max)
        ax.set_xscale('log')
        ax.set_yscale('log')
    
        ax.plot(self.xmodel, self.xymodel, color='C1', zorder=1)
    
        fig.savefig(fout)
        plt.close()

    #----------------------------------------------------------------------------------------------------
