import glob
import os
from Swift_XRT_WTmode import *

#=============================================================================#
'''
PURPOSE:
--------
Perform data reduction tasks for a Swift/XRT observation in WT mode.

NOTES:
------
- Must initialize HEAsoft
- File naming conventions are imposed
'''
#=============================================================================#
def reduce_data(obsID, regsrc, regbkg,
        ra='"16 54 00.137"', dec='"-39 50 44.90"',
        tag='', grade='0-2', grpmin=20, syserr=0.03, phacuts=None,
        datdir='/Users/salvesen/research/fcolabs/data/',
        wrtdir='/Users/salvesen/research/fcolabs/data/clean/',
        runxrt=True, runimg=True, runsrc=True, runbkg=True,
        runarf=True, runrmf=True, rungrp=True):
    """
    Inputs:
    -------
    obsID   : Swift observation identification number (string)
    regsrc  : Region file for source extraction (ds9 output in WCS format)
    regbkg  : Region file for background extraction (ds9 output in WCS format)
    ra      : Source RA (string, '"hh mm ss.s"')
    dec     : Source DEC (string, '"deg arcmin arcsec"')
    tag     : Additional tag on output filenames (string, placed b4 obsID tag)
    grade   : Grades to extract from the Swift/XRT WT mode data (string)
    grpmin  : Group the spectrum to have this minimum number of counts per bin
    syserr  : Add this fractional systematic error to the data
    phacuts : Range of photon energy channels to keep [phalcut, phahcut]
              (1 XRT channel = 10 eV)
    datdir  : Directory containing the raw Swift data (obsID will be appended)
    wrtdir  : Directory to write the output files to (obsID will be appended)

    True/False Flags:
    -----------------
    runxrt : Run xrtpipline() <-- time intensive
    runimg : Run extract_image() on the full event-lists
    runsrc : Run extract_spectrum() on the source region
    runbkg : Run extract_spectrum() on the background region
    runarf : Run xrtmkarf()
    runrmf : Run downloadrmf()
    rungrp : Run grppha()
    """
    # Swift observation ID number
    print("Dataset: ", obsID)

    # Directories: Data, Write
    datdir = datdir + obsID + '/'
    wrtdir = wrtdir + obsID + '/'

    # Run xrtpipeline (also create exposure map needed to build the ARF later)
    if (runxrt is True):
        logxrt = wrtdir + 'xrtpipeline_' + tag + obsID + '.log'
        print("...Cleaning data...")
        xrtpipeline(indir         = datdir,
                    steminputs    = 'sw' + obsID,
                    outdir        = wrtdir,
                    srcra         = ra,
                    srcdec        = dec,
                    createexpomap = 'yes',
                    cleanup       = 'yes',
                    logfile       = logxrt)

    # Extract image
    if (runimg is True):
        print("...Extracting image...")
        logimg = wrtdir + 'image_' + tag + obsID + '.log'
        extract_image(infile  = 'sw' + obsID + 'xwtw2po_cl.evt',
                      indir   = wrtdir,
                      outfile = wrtdir + 'image_' + tag + obsID + '.img',
                      xcofile = wrtdir + 'image_' + tag + obsID + '.xco',
                      region  = regsrc,
                      phacuts = phacuts,
                      logfile = logimg)

    # Extract source spectrum
    if (runsrc is True):
        outfile  = wrtdir + 'source_' + tag + obsID + '.pha'
        xcofile  = wrtdir + 'source_' + tag + obsID + '.xco'
        logsrc   = wrtdir + 'source_' + tag + obsID + '.log'
        print("...Extracting source spectrum...")
        extract_spectrum(infile     = 'sw' + obsID + 'xwtw2po_cl.evt',
                         indir      = wrtdir,
                         outfile    = outfile,
                         xcofile    = xcofile,
                         region     = regsrc,
                         phacuts    = phacuts,
                         grade      = grade,
                         logfile    = logsrc,
                         background = False)

    # Extract background spectrum
    if (runbkg is True):
        infile   = 'sw' + obsID + 'xwtw2po_cl.evt'
        indir    = wrtdir
        outfile  = wrtdir + 'background_' + tag + obsID + '.pha'
        xcofile  = wrtdir + 'background_' + tag + obsID + '.xco'
        logbkg   = wrtdir + 'background_' + tag + obsID + '.log'
        print("...Extracting background spectrum...")
        extract_spectrum(infile     = infile,
                         indir      = indir,
                         outfile    = outfile,
                         xcofile    = xcofile,
                         region     = regbkg,
                         phacuts    = phacuts,
                         grade      = grade,
                         logfile    = logbkg,
                         background = True)

    # Build ancillary response file
    # (Taking care to use the exposure map created by xrtpipeline)
    if (runarf is True):
        phafile  = wrtdir + 'source_' + tag + obsID + '.pha'
        outfile  = wrtdir + 'source_' + tag + obsID + '.arf'
        expofile = wrtdir + 'sw' + obsID + 'xwtw2po_ex.img'
        srcx     = -1
        srcy     = -1
        psfflag  = 'yes'
        psffile  = 'CALDB'
        rmffile  = 'CALDB'
        logarf   = wrtdir + 'arf_' + tag + obsID + '.log'
        print("...Building ancillary response file...")
        xrtmkarf(phafile,
                 outfile,
                 expofile,
                 srcx    = srcx,
                 srcy    = srcy,
                 psfflag = psfflag,
                 psffile = psffile,
                 rmffile = rmffile,
                 logfile = logarf)

    # Download the proper response matrix file
    if (runrmf is True):
        print("...Downloading response matrix file...")
        downloadrmf(logfile = logarf,
                    outdir  = wtrdir)

    # Group the source spectrum
    if (rungrp is True):
        grptag   = 'grp' + str(grpmin) + '_' 'serr' + str(int(round(syserr*100))) + '_'
        infile   = wrtdir + 'source_' + tag + obsID + '.pha'
        outfile  = wrtdir + 'source_' + grptag + tag + obsID + '.pha'
        backfile = wrtdir + 'background_' + tag + obsID + '.pha'
        respfile = glob.glob(wrtdir + '*.rmf')[0]
        ancrfile = wrtdir + 'source_' + tag + obsID + '.arf'
        loggrp   = wrtdir + 'source_' + grptag + tag + obsID + '.log'
        print("...Grouping the source spectrum...")
        grppha(infile   = infile,
               outfile  = outfile,
               backfile = backfile,
               respfile = respfile,
               ancrfile = ancrfile,
               grpmin   = grpmin,
               syserr   = syserr,
               phacuts  = phacuts,
               logfile  = loggrp)

'''
# NEEDS WORK! Maybe just dump the output of grep *.log to a file
# Search for ERRORS in the log files
print("...Searching log files for errors...")
logfiles = [logxrt, logimg, logsrc, logbkg, logarf, loggrp]
for file in logfiles:
    print("......Log file: ", file)
    a = os.popen("grep -i error" + file)
    print(a.size())
'''
#=============================================================================#
