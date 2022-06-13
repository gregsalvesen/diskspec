import numpy as np
from plot_emcee import *

'''
'''

# Input filenames
rdir  = "/Users/salvesen/research/fcolabs/results/chains/"
fmcmc = rdir + "iFlat_iJet_chain_emcee_kerrbb.hdf5"#"xspec_chain_emcee_kerrbb.hdf5"

# Output filenames
dout     = "/Users/salvesen/research/fcolabs/results/plots/"
tag_walk = "iFlat_iJet_Walker_Paths_"#"xspec_Walker_Paths_"
fwalk_a  = dout + tag_walk + "a.png"


# Make the parameter object
class Parameter(object):
    def __init__(self, label, index):
        self.label = label
        self.index = index
def make_parameter(label, index):
    parameter = Parameter(label, index)
    return parameter


# Instance of the plot_emcee class
par_a = make_parameter(label=r'$a$', index=0)
inst  = plot_emcee(fmcmc=[fmcmc], par=par_a)

# Plots
inst.walker_paths(fout=fwalk_a)
