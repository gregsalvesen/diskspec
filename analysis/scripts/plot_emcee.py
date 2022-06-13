import numpy as np
import pylab as plt
import h5py

'''
PURPOSE:
--------
Plot the one-dimensional posterior for the disk inclination.

INPUTS:
-------
fmcmc - List of HDF5 filenames containing the MCMC chains and lnprobs output by xspec_emcee
par   - Attributes: index, name, minval, maxval
'''

#====================================================================================================
class plot_emcee:

    #------------------------------------------------------------------------------------------------
    # INITIALIZATION
    def __init__(self, fmcmc, par):

        # Parameter of interest
        self.par = par

        # Initialize lists that will contain...
        self.chains   = []  # MCMC chains
        self.lnprobs  = []  # MCMC lnprobs
        self.Nwalkers = []  # Number of walkers
        self.Nsteps   = []  # Number of steps/iterations
        self.Npars    = []  # Number of parameters
        self.steps    = []  # Array of steps/iterations

        # Loop through the input HDF5 files
        self.Nmcmc = len(fmcmc)
        for i in range(self.Nmcmc):
        
            # Extract the EMCEE chains and lnprobs from the current HDF5 file
            f = h5py.File(fmcmc[i], 'r')
            self.chains.append(f['chain'][:,:,:])  # [walker, iteration, parameter]
            self.lnprobs.append(f['lnprob'][:,:])  # [walker, iteration]
            f.close()

            # Number of walkers, number of steps, array of steps
            self.Nwalkers.append(self.chains[i].shape[0])
            self.Nsteps.append(self.chains[i].shape[1])
            self.Npars.append(self.chains[i].shape[2])
            self.steps.append(np.arange(self.Nsteps[i]) + 1)

    #------------------------------------------------------------------------------------------------
    # HELPER FUNCTIONS

    #................................................................................................
    def TopHat_pdf(self, x, xmin=0.0, xmax=1.0):
        '''
        Inputs:
        -------
        x    - Value or array/list of bins at which to compute the PDF
        xmin - Minimum value of the top hot distribution (inclusive)
        xmax - Maximum value of the top hat distribution (inclusive)
    
        Output:
        -------
        Top hat probability distribution function (PDF)
        '''
        if ((x >= xmin) and (x <= xmax)):
            pdf = 1.0 / (xmax - xmin)
        else:
            pdf = 0.0
        return pdf
        
    #................................................................................................
    def Gaussian_pdf(self, x, mu=0.0, sigma=1.0, A=None):
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

    #------------------------------------------------------------------------------------------------
    # PLOT: Acceptance Fraction (or something similar)
    #def acceptance_fraction(self, fout):

    #------------------------------------------------------------------------------------------------
    # PLOT: Walker Paths
    def walker_paths(self, fout):

        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.set_xlabel(r'${\rm Iteration\ Step}$')
        ax.set_ylabel(self.par.label)

        for i in np.arange(self.Nmcmc):
            for j in np.arange(self.Nwalkers[i]):
                thisChain = self.chains[i][j,:,self.par.index]
                ax.plot(self.steps[i], thisChain, color='black', linewidth=1, alpha=0.1)
        
        fig.savefig(fout)
        plt.close()

    #------------------------------------------------------------------------------------------------
    # PLOT: One-Dimensional Parameter Posterior Marginalized Over All Other Parameters
    #def posterior_1D(self, fout, paridx):

    #------------------------------------------------------------------------------------------------
    # PLOT: Two-Dimensional Parameter Posterior Marginalized Over All Other Parameters
    #def posterior_2D(self, fout, paridx1, paridx2):

    #------------------------------------------------------------------------------------------------
    # PLOT: All Parameter Posteriors (2D and 1D)
    #def posterior_mcmc(self, fout):

