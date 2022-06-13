#!/usr/bin/env python

"""
Use EMCEE to do MCMC in Xspec.
Jeremy Sanders 2012-2016

Requires Python 2.7+, numpy, scipy, h5py and emcee
"""

from __future__ import print_function, division

import sys
import argparse
import time
import re
import itertools

import h5py
import numpy as N
import emcee

import xspec_pool

def get_initial_parameters(parameters, nwalkers):
    """Construct list of initial parameter values for each walker."""
    p0 = []
    for walker in xrange(nwalkers):
        pwalker = []
        # for each walker, use initial parameters based on parameter
        # and delta parameter
        for par in parameters:
            width = par.delta
            swidth = par.sigma*0.1
            if swidth > 0 and swidth < width:
                # use sigma if delta is badly adjusted
                width = swidth
            
            v = N.random.normal(par.initval, width)
            # clip to hard range
            v = N.clip(v, par.minval, par.maxval)
            pwalker.append(v)
        p0.append( N.array(pwalker) )
    return N.array(p0)

def expand_systems(systems):
    """Allow system*N syntax in systems."""
    out = []
    for s in systems:
        m = re.match(r'([A-Za-z0-9]+)\*([0-9]+)', s)
        if m:
            out += [m.group(1)]*int(m.group(2))
        else:
            out.append(s)
    return out

def do_mcmc(xcm, nwalkers=100, nburn=100, niters=1000, systems = ['localhost'],
            outchain='out.dat', outhdf5='out.hdf5', debug=False,
            continuerun=False, autosave=True,
            nochdir=False, initialparameters=None,
            lognorm=False, chunksize=4, parameterpriors=None):
    """Do the actual MCMC process."""

    print("Loading XCM file")
    xspec = xspec_pool.Xspec(
        xcm, expand_systems(systems), debug=debug, nochdir=nochdir)

    if lognorm:
        print("Using prior equivalent to log parameter")
        xspec.log_norms_priors()

    if not initialparameters:
        print("Getting initial parameters")
        p0 = get_initial_parameters(xspec.thawedpars, nwalkers)
    else:
        print("Loading initial parameters from", initialparameters)
        p0 = N.loadtxt(initialparameters)
        
    #----------------------------------------------------------------------------------------------------
    if not parameterpriors:
        print("Using uniform priors for all thawed parameters")
    else:
        print("Using Gaussian priors for parameters specifed in", parameterpriors)
        xspec.parameter_priors(parameterpriors)
    #----------------------------------------------------------------------------------------------------

    ndims = p0.shape[1]
    pool = xspec_pool.XspecPool(xspec)

    # sample the mcmc
    sampler = emcee.EnsembleSampler(nwalkers, ndims, None, pool=pool)

    print("Starting MCMC")
    if not continuerun and nburn > 0:
        # burn in
        print("Burn in period started")
        pos, prob, state = sampler.run_mcmc(p0, nburn)
        sampler.reset()
        print("Burn in period finished")
    else:
        # no burn in
        state = None
        pos = p0

    if not continuerun:
        # create new datasets, extensible along number of iterations
        hdf5file = h5py.File(outhdf5, "w")
        chain = hdf5file.create_dataset(
            "chain",
            (nwalkers, niters, ndims),
            maxshape=(nwalkers, None, ndims))
        lnprob = hdf5file.create_dataset(
            "lnprob",
            (nwalkers, niters),
            maxshape=(nwalkers, None))
        start = 0

    else:
        print("Continuing from existing chain in", outhdf5)

        hdf5file = h5py.File(outhdf5, "r+")
        chain = hdf5file["chain"]
        lnprob = hdf5file["lnprob"]

        start = chain.attrs["count"]
        pos = N.array(chain[:, start-1, :])
        print("Restarting at iteration", start)

        chain.resize((nwalkers, niters, ndims))
        lnprob.resize((nwalkers, niters))

    # iterator interface allows us to trap ctrl+c and know where we are
    lastsave = time.time()
    index = start
    try:
        for p, l, s in sampler.sample(
            pos, rstate0=state, storechain=False,
            iterations=niters-start):

            chain[:, index, :] = p
            lnprob[:, index] = l
            index += 1

            if autosave and time.time() - lastsave > 60*10:
                chain.attrs["count"] = index
                hdf5file.flush()

    except KeyboardInterrupt:
        chain.attrs["count"] = index
        print("Ctrl+C pressed - ending")

    else:
        chain.attrs["count"] = index
        print("Writing chain", outchain)
        with open(outchain, "w") as chainf:
            write_xspec_chain(
                chainf, chain, lnprob,
                xspec.thawedpars, xspec.xspec_thawed_idxs(),
                nwalkers)

    hdf5file.close()

def write_xspec_chain(chainf, chain, lnprob, params, paridxs, nwalkers):
    """Write an xspec text chain file to file object chainf.
    """

    chainf.write('! Markov chain file generated by xspec "chain" command.\n')
    chainf.write('!    Do not modify, else file may not reload properly.\n')

    nwalkers, niters, ndims = chain.shape

    # real length could be shorter
    niters = chain.attrs["count"]

    chainf.write('!Length: %i  Width: %i\n' % (niters*nwalkers, ndims+1))
    chainf.write('!Type: GoodmanWeare\n')
    chainf.write('!NWalkers: %i\n' % nwalkers)

    # header for contents of file
    hdr = []
    for par, idx in itertools.izip(params, paridxs):
        hdr.append("%s %s %s" % (
                idx, par.name,
                par.unit if par.unit else "0"))
    hdr.append("Likelihood")
    chainf.write('!%s\n' % ' '.join(hdr))

    fmt = '\t'.join(['%g']*ndims)

    # then each walker separately
    for wi in xrange(nwalkers):
        chainw = chain[wi, :, :]
        statw = lnprob[wi, :]

        # write output
        for pars, stat in itertools.izip(chainw, statw):
            line = fmt % tuple(pars) + '\t' + '%g' % stat + '\n'
            chainf.write(line)

def main():
    """Main program."""

    p = argparse.ArgumentParser(
        description="Xspec MCMC with EMCEE. Jeremy Sanders 2012-2016.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("xcm", metavar="XCM",
                   help="Input XCM file")
    p.add_argument("--niters", metavar="N", type=int, default=5000,
                   help="Number of iterations")
    p.add_argument("--nburn", metavar="N", type=int, default=500,
                   help="Number of burn iterations")
    p.add_argument("--nwalkers", metavar="N", type=int, default=50,
                   help="Number of walkers")
    p.add_argument("--systems", default="localhost", metavar="LIST",
                   help="Space separated list of systems to run on")
    p.add_argument("--output-hdf5", default="emcee.hdf5", metavar="FILE",
                   help="Output HDF5 file")
    p.add_argument("--output-chain", default="emcee.chain", metavar="FILE",
                   help="Output text file")
    p.add_argument("--continue-run",  action="store_true", default=False,
                   help="Continue from an existing chain (in HDF5)")
    p.add_argument("--debug", action="store_true", default=False,
                   help="Create xspec log files")
    p.add_argument("--no-chdir", action="store_true", default=False,
                   help="Do not chdir to xcm file directory before execution")
    p.add_argument("--initial-parameters", metavar="FILE",
                   help="Provide initial parameters")
    p.add_argument("--log-norm", action="store_true", default=False,
                   help="log norm values during MCMC")
    p.add_argument('--chunk-size', metavar='N', type=int, default=4,
                   help='Number of sets of parameters to pass to xspec')
    p.add_argument("--parameter-priors", metavar="FILE",
                   help="Provide parameter priors")

    args = p.parse_args()

    sampler = do_mcmc(
        args.xcm,
        systems = args.systems.split(),
        nwalkers = args.nwalkers,
        nburn = args.nburn,
        niters = args.niters,
        outchain = args.output_chain,
        outhdf5 = args.output_hdf5,
        continuerun = args.continue_run,
        debug = args.debug,
        nochdir = args.no_chdir,
        initialparameters = args.initial_parameters,
        lognorm = args.log_norm,
        chunksize = args.chunk_size,
        parameterpriors = args.parameter_priors,
    )

    print("Done")

if __name__ == '__main__':
    main()
