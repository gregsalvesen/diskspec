import numpy as np
import subprocess
import os
from astropy.io import fits

#=============================================================================#
'''
PURPORSE:
---------
Perform various data reduction tasks for Swift/XRT in Windowed Timing (WT) mode.

FUNCTIONS:
----------
Data Reduction Functions
- xrtpipeline(indir, steminputs, outdir, srcra, srcdec, [createexpomap, cleanup, logfile])
- xrtmkarf(phafile, outfile, expofile, [srcx, srcy, psfflag, psffile, rmffile, logfile])
- downloadrmf(logfile, outdir)
- extract_image(infile, indir, outfile, xcofile, [region, logfile])
- extract_spectrum(infile, indir, outfile, xcofile, region, [grade, logfile])
- grppha(infile, outfile, backfile, respfile, ancrfile, [grpmin, syserr, logfile])

Helper Functions
- fkeypar(fitsfile, keyword, [index])
- fparkey(value, fitsfile, keyword)
- write_xco(xcofile, lines)
- counts(fitsfile, [phacuts])
- offaxis_angle(fitsfile)
- goodtime(fitsfile)
- count_rate(fitsfile, [phacuts])
- bkgsub_counts(fsrc, fbkg, [phacuts, stdev])
- bkgsub_count_rate(fsrc, fbkg, [phacuts, stdev])
- optimize_region()

NOTES:
------
- Must initialize HEAsoft
- Swift XRT Data Reduction Guide: http://www.swift.ac.uk/analysis/xrt/files/xrt_swguide_v1_2.pdf
- Xselect User Guide: https://heasarc.gsfc.nasa.gov/ftools/xselect/node1.html
- Xselect Keys: http://www.darts.isas.ac.jp/pub/legacy.gsfc.nasa.gov/software/lheasoft/lheasoft6.0/linux_mac/heasarc/headas/powerpc-unknown-linux-gnu/bin/xselect.key
'''
#=============================================================================#
# DATA REDUCTION FUNCTIONS

#-----------------------------------------------------------------------------#
def xrtpipeline(indir, steminputs, outdir, srcra, srcdec, \
                createexpomap='yes', cleanup='no', logfile=None):
    """
    Produce XRT Level 2 cleaned event-files from Level 1 event-files
    (using xrtpipeline defaults)
    
    indir          : Directory containing the Level 1 data
    steminputs     : Root filename (will be 'sw' + the observation ID number)
    outdir         : Output directory for the (many) generated Level 2 files
    srcra          : Source RA (string, format: '"hh mm ss.s"')
    srcdec         : Source DEC (string, format: '"hh mm ss.s"')
    createexpomap  : Flag to create an exposure map (yes/no), which is needed by the ARF later
    cleanup        : Flag to cleanup the (many) secondary files that are generated (yes/no)
    logfile        : Output log filename (w/ path)
    """
    cmd = 'xrtpipeline' + \
          ' indir=' + indir + ' steminputs=' + steminputs + ' outdir=' + outdir + \
          ' srcra=' + srcra + ' srcdec=' + srcdec + \
          ' createexpomap=' + createexpomap + ' cleanup=' + cleanup + ' clobber=yes'
    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def xrtmkarf(phafile, outfile, expofile, srcx=-1, srcy=-1, psfflag='yes', \
    psffile='CALB', rmffile='CALDB', logfile=None):
    """
    Build an ancillary response file (ARF)

    phafile  : Source spectrum (ungrouped) filename (w/ path)
    outfile  : Output ARF filename (w/ path)
    expofile : Exposure map filename (w/ path) produced from xrtpipeline()
    srcx     : Source x-position in detector coordinates (use -1 to use coordinates in PHA input file)
    srcy     : Source x-position in detector coordinates (use -1 to use coordinates in PHA input file)
    psfflag  : Flag to correct the PSF for point sources (yes/no)
    psffile  : PSF filename (use 'CALDB')
    rmffile  : RMF filename (use 'CALDB')
    logfile  : Output log filename (w/ path)
    """
    cmd = 'xrtmkarf' + \
          ' phafile=' + phafile + ' outfile=' + outfile + ' expofile=' + expofile + \
          ' srcx=' + str(srcx) + ' srcy=' + str(srcy) + \
          ' psfflag=' + psfflag + ' rmffile=' + rmffile + ' clobber=yes'
    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def downloadrmf(logfile, outdir):
    """
    Download the correct response matrix file (RMF) from the output in the xrtarf logfile
    (Alternatively, could use quzcif: https://heasarc.gsfc.nasa.gov/ftools/caldb/quzcif.html)

    logfile : Log filename (w/ path) output by xrtmkarf() and containing the RMF info
    outdir  : Output directory for the downloaded RMF
    """
    rmfftp  = os.popen("grep .rmf " + logfile).read().split("'")[1]
    rmffile = rmfftp.split('/')[-1]
    cmd     = 'wget ' + rmfftp + ' -O ' + outdir + '/' + rmffile + ' -nc'
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def extract_spectrum(infile, indir, outfile, xcofile, region, grade='0-2', phacuts=None, background=False, logfile=None):
    """
    Extract a spectrum by writing an XCO file and executing it in Xselect

    infile     : Cleaned Level 2 event filename (_cl.evt) w/o the path
    indir      : Path to the infile
    outfile    : Output spectrum filename (w/ path)
    xcofile    : Output XCO filename (w/ path) used to run Xselect commands
    region     : Region file to extract (ds9 output in WCS format)
    grade      : Grades to extract in WT mode (string, choose '0-2' or '0')
    phacuts    : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    logfile    : Output log filename (w/ path)
    background : Flag for if this is a background spectrum (if set to True, subtract 1 from BACKSCAL)
    """
    # Initialize the list of lines to write to an XCO file
    lines = []
    # Name the Xselect session
    lines.append('xsel')
    # Read in the events file
    lines.append('read' + ' read_what=event' + ' infiles=' + infile + ' data_dir=' + indir + \
                 ' reset_inst=yes' + ' reset_miss=yes' + ' reset_datamode=yes')
    # Region to extract
    lines.append('filter' + ' filter_what=region' + ' region=' + region)
    # Grade to extract
    lines.append('filter' + ' filter_what=grade' + ' grade=' + grade)
    # PHA Channels to extract
    if (phacuts is not None):
        lines.append('filter' + ' filter_what=pha_cutoff' + \
                     ' phalcut=' + str(phacuts[0]) + ' phahcut=' + str(phacuts[1]))
    # Extract the spectrum
    lines.append('extract' + ' bin_what=spectrum')
    # Save the spectrum
    lines.append('save' + ' save_what=spectrum' + ' save_str=' + outfile + ' clobberit=yes')
    # Exit session and quit Xselect
    lines.append('exit' + ' save_session=no')
    lines.append('quit')
    # Write the lines above to an XCO file
    write_xco(xcofile=xcofile, lines=lines)

    # Use Xselect to run the commands in the XCO file
    cmd = 'xselect @' + xcofile
    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

    # Change the BACKSCAL keyword (if specificed by the user)
    arcsec2pix  = 1.0 / 2.357  # 1 XRT pixel = 2.357"
    rin_arcsec  = float(os.popen("grep annulus " + region).read().split(',')[-2].split('"')[0])
    rout_arcsec = float(os.popen("grep annulus " + region).read().split(',')[-1].split('"')[0])
    rin_pix     = rin_arcsec * arcsec2pix
    rout_pix    = rout_arcsec * arcsec2pix
    value       = rout_arcsec - rin_arcsec
    if (background is True): value = value - 1.0
    fparkey(value=value, fitsfile=outfile, keyword='backscal')

#-----------------------------------------------------------------------------#
def grppha(infile, outfile, backfile, respfile, ancrfile, grpmin=20, syserr=0, phacuts=None, logfile=None):
    """
    Group the counts in a source spectrum

    infile   : Input source spectrum filename (w/ path)
    outfile  : Output grouped source spectrum filename (w/ path)
    backfile : Background filename (w/ path)
    respfile : Response matrix filename (w/ path)
    ancrfile : Acnillary response filename (w/ path)
    grpmin   : Minimum counts/bin grouping to apply (20 is standard for chi-square fitting)
    syserr   : Fractional systematic error to include
    phacuts  : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    logfile  : Output log filename (w/ path)
    """
    cmd = 'grppha' + \
          ' infile=' + infile + ' outfile=' + outfile + ' clobber=true' + ' comm="'

    if (phacuts is None):
        # Flag channels 0-29 as bad and add systematic errors to all channels
        cmd += 'bad 0-29' + ' & ' + 'systematics ' + '0-1023 ' + str(syserr) + ' & '
    else:
        # Flag channels phacuts[0]-29 as bad
        if (phacuts[0] < 29):
            cmd += 'bad ' + str(phacuts[0]) + '-29' + ' & '
        # Add systematic errors to phacuts[0]-phacuts[1]
        cmd += 'systematics ' + str(phacuts[0]) + '-' + str(phacuts[1]) + ' ' + str(syserr) + ' & '

    cmd += 'group min ' + str(grpmin) + ' & ' + \
          'chkey backfile ' + backfile + ' & ' + \
          'chkey ancrfile ' + ancrfile + ' & ' + \
          'chkey respfile ' + respfile + ' & ' + \
          'write !' + outfile + ' & ' + \
          'quit' + '"'

    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def extract_image(infile, indir, outfile, xcofile, region=None, phacuts=None, logfile=None):
    """
    Extract an image by writing an XCO file and executing it in Xselect

    infile   : Cleaned Level 2 event filename (_cl.evt) w/o the path
    indir    : Path to the infile
    outfile  : Output image filename (w/ path)
    xcofile  : Output XCO filename (w/ path) used to run Xselect commands
    region   : Region file to extract (ds9 output in WCS format)
    phacuts  : Photon energy bounds to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    logfile  : Output log filename (w/ path)
    """
    # Initialize the list of lines to write to an XCO file
    lines = []
    # Name the Xselect session
    lines.append('xsel')
    # Read in the events file
    lines.append('read' + ' read_what=event' + ' infiles=' + infile + ' data_dir=' + indir + \
                 ' reset_inst=yes' + ' reset_miss=yes' + ' reset_datamode=yes')
    # Apply region filter
    if (region is not None):
        lines.append('filter' + ' filter_what=region' + ' region=' + region)
    # Apply energy filter
    if (phacuts is not None):
        lines.append('filter' + ' filter_what=pha_cutoff' + \
                     ' phalcut=' + str(phacuts[0]) + ' phahcut=' + str(phacuts[1]))
    # Extract the image
    lines.append('extract' + ' bin_what=image')
    # Save the image
    lines.append('save' + ' save_what=image' + ' save_str=' + outfile + ' clobberit=yes')
    # Exit session and quit Xselect
    lines.append('exit' + ' save_session=no')
    lines.append('quit')
    # Write the lines above to an XCO file
    write_xco(xcofile=xcofile, lines=lines)

    # Use Xselect to run the commands in the XCO file
    cmd = 'xselect @' + xcofile
    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def extract_curve(infile, indir, outfile, xcofile, region, grade='0-2', phacuts=None, background=False, logfile=None):
    """
    Extract a light curve by writing an XCO file and executing it in Xselect
    ??? CAN WE CHANGE THE BINSIZE ???

    infile     : Cleaned Level 2 event filename (_cl.evt) w/o the path
    indir      : Path to the infile
    outfile    : Output spectrum filename (w/ path)
    xcofile    : Output XCO filename (w/ path) used to run Xselect commands
    region     : Region file to extract (ds9 output in WCS format)
    grade      : Grades to extract in WT mode (string, choose '0-2' or '0')
    phacuts    : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    logfile    : Output log filename (w/ path)
    background : Flag for if this is a background spectrum (if set to True, subtract 1 from BACKSCAL)
    """
    # Initialize the list of lines to write to an XCO file
    lines = []
    # Name the Xselect session
    lines.append('xsel')
    # Read in the events file
    lines.append('read' + ' read_what=event' + ' infiles=' + infile + ' data_dir=' + indir + \
                 ' reset_inst=yes' + ' reset_miss=yes' + ' reset_datamode=yes')
    # Region to extract
    lines.append('filter' + ' filter_what=region' + ' region=' + region)
    # Grade to extract
    lines.append('filter' + ' filter_what=grade' + ' grade=' + grade)
    # PHA Channels to extract
    if (phacuts is not None):
        lines.append('filter' + ' filter_what=pha_cutoff' + \
                     ' phalcut=' + str(phacuts[0]) + ' phahcut=' + str(phacuts[1]))
    # Extract the light curve
    lines.append('extract' + ' bin_what=curve')
    # Save the light curve
    lines.append('save' + ' save_what=curve' + ' save_str=' + outfile + ' clobberit=yes')
    # Exit session and quit Xselect
    lines.append('exit' + ' save_session=no')
    lines.append('quit')
    # Write the lines above to an XCO file
    write_xco(xcofile=xcofile, lines=lines)

    # Use Xselect to run the commands in the XCO file
    cmd = 'xselect @' + xcofile
    if (logfile is not None): cmd += '&> ' + logfile
    subprocess.call(cmd, shell=True)

    # Change the BACKSCAL keyword (if specificed by the user)
    arcsec2pix  = 1.0 / 2.357  # 1 XRT pixel = 2.357"
    rin_arcsec  = float(os.popen("grep annulus " + region).read().split(',')[-2].split('"')[0])
    rout_arcsec = float(os.popen("grep annulus " + region).read().split(',')[-1].split('"')[0])
    rin_pix     = rin_arcsec * arcsec2pix
    rout_pix    = rout_arcsec * arcsec2pix
    value       = rout_arcsec - rin_arcsec
    if (background is True): value = value - 1.0
    print("VALUE = ", value)
    fparkey(value=value, fitsfile=outfile, keyword='backscal')

    # Exposure correction...not sure I fully understand this yet
    #http://www.swift.ac.uk/analysis/xrt/lccorr.php

#-----------------------------------------------------------------------------#
#def xrtlccorr():


#=============================================================================#
# HELPER FUNCTIONS

#-----------------------------------------------------------------------------#
def fkeypar(fitsfile, keyword, index=0):
    """
    Return the value (string) of a keyword in a FITS file
    
    fitsfile : FITS filename (w/ path)
    keyword  : Keyword name in the FITS header to alter
    """
    cmd = 'fkeypar' + \
          ' fitsfile=' + fitsfile + '[' + str(index) + ']' + ' keyword=' + keyword
    subprocess.call(cmd, shell=True)
    cmd = 'pget fkeypar value'
    value = os.popen(cmd).read().split('\n')[0]
    return value

#-----------------------------------------------------------------------------#
def fparkey(value, fitsfile, keyword):
    """
    Change the BACKSCAL keyword for a source or background spectrum
    - For source, set value = rout - rin
    - For background, value = rout - rin - 1 <-- because the end of window pixel is flagged as bad
    
    value    : Value to set for the FITS keyword
    fitsfile : FITS filename (w/ path)
    keyword  : Keyword name in the FITS header to alter
    """
    cmd = 'fparkey' + \
          ' value=' + str(value) + ' fitsfile=' + fitsfile + ' keyword=' + keyword
    subprocess.call(cmd, shell=True)

#-----------------------------------------------------------------------------#
def write_xco(xcofile, lines):
    """
    Write an XCO file line-by-line
    
    xcofile : Write to this filename (w/ path)
    lines   : List of lines to write
    """
    f = open(xcofile, 'w')
    for item in lines:
        f.write(item + '\n')
    f.close()

#-----------------------------------------------------------------------------#
def offaxis_angle(fitsfile):
    """
    Calculate the off-axis angle [degrees]
    
    fitsfile : FITS filename (w/ path)
    """
    ra_obj  = float(fkeypar(fitsfile=fitsfile, keyword='RA_OBJ', index=0))   # [degrees]
    dec_obj = float(fkeypar(fitsfile=fitsfile, keyword='DEC_OBJ', index=0))  # [degrees]
    ra_pnt  = float(fkeypar(fitsfile=fitsfile, keyword='RA_PNT', index=0))   # [degrees]
    dec_pnt = float(fkeypar(fitsfile=fitsfile, keyword='DEC_PNT', index=0))  # [degrees]
    offaxis = np.sqrt((ra_obj - ra_pnt)**2 + (dec_obj - dec_pnt)**2)         # [degrees]
    return offaxis

#-----------------------------------------------------------------------------#
def goodtime(fitsfile):
    """
    Return the total "good time interval" (float) from an extracted spectrum

    fitsfile : FITS filename (w/ path)
    """
    f      = fits.open(fitsfile)
    gti    = f['GTI']
    tstart = gti.data['START']
    tstop  = gti.data['STOP']
    tgood  = np.sum(tstop - tstart)
    f.close()
    return tgood

#-----------------------------------------------------------------------------#
def counts(fitsfile, phacuts=None, stdev=True):
    """
    Return the number of counts (int) from an extracted spectrum with optional energy filter
    (Optionally return the standard deviation in the number of counts)

    fitsfile : FITS filename (w/ path)
    phacuts  : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    stdev    : True/False flag to return the standard deviation in the total counts
    """
    f        = fits.open(fitsfile)
    spectrum = f['SPECTRUM']
    channels = spectrum.data['CHANNEL']
    counts   = spectrum.data['COUNTS']
    f.close()
    if (phacuts is not None):
        ilo    = np.where(channels >= phacuts[0])[0][0]
        ihi    = np.where(channels <= phacuts[1])[0][-1]
        counts = counts[ilo:ihi+1]

    # Get the total counts and standard deviation
    Ncounts = np.sum(counts)
    stdcnts = np.std(counts)
    if (stdev is True): return Ncounts, stdcnts
    if (stdev is False): return Ncounts

#-----------------------------------------------------------------------------#
def count_rate(fitsfile, phacuts=None, stdev=True):
    """
    Return the count rate (counts/second) w/ stdev from an extracted spectrum with optional energy filter

    fitsfile : FITS filename (w/ path)
    phacuts  : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    stdev    : True/False flag to return the standard deviation in the count rate
    """
    # Calculate the total and standard deviation for the counts
    Ncounts, stdcnts = counts(fitsfile=fitsfile, phacuts=phacuts, stdev=stdev)
    # Calculate the total and standard deviation for the count rate
    tgood    = goodtime(fitsfile=fitsfile)
    cnt_rate = float(Ncounts) / tgood
    std_rate = float(stdcnts) / tgood
    if (stdev is True): return cnt_rate, std_rate
    if (stdev is False): return cnt_rate

#-----------------------------------------------------------------------------#
def bkgsub_counts(fsrc, fbkg, phacuts=None, stdev=True):
    """
    Calculate the background-subtracted counts (and standard deviation)

    fsrc    : Source spectrum file
    fbkg    : Background spectrum file
    phacuts : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    stdev   : True/False flag to return the standard deviation in the count rate
    """
    # Collect the counts (and standard deviation) for the source region and background region
    cnts_src, stds_src = counts(fitsfile=fsrc,  phacuts=phacuts, stdev=True)
    cnts_bkg, stds_bkg = counts(fitsfile=fbkg,  phacuts=phacuts, stdev=True)

    # Collect the BACKSCAL keywords
    backscal_src = float(fkeypar(fitsfile=fsrc, keyword='BACKSCAL', index=1))
    backscal_bkg = float(fkeypar(fitsfile=fbkg, keyword='BACKSCAL', index=1))

    # Background subtracted counts (and standard deviation)
    bkgsub_cnts = cnts_src - (cnts_bkg * (backscal_src / backscal_bkg))  # [counts / second]
    bkgsub_stds = stds_src - (stds_bkg * (backscal_src / backscal_bkg))  # [counts / second]
    if (stdev is True): return bkgsub_cnts, bkgsub_stds
    if (stdev is False): return bkgsub_cnts

#-----------------------------------------------------------------------------#
def bkgsub_count_rate(fsrc, fbkg, phacuts=None, stdev=True):
    """
    Calculate the background-subtracted count rate (and standard deviation)

    fsrc    : Source spectrum file
    fbkg    : Background spectrum file
    phacuts : Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    stdev   : True/False flag to return the standard deviation in the count rate
    """
    # Collect the count rate (and standard deviation) for the source region and background region
    cntrate_src, stdrate_src = count_rate(fitsfile=fsrc,  phacuts=phacuts, stdev=True)
    cntrate_bkg, stdrate_bkg = count_rate(fitsfile=fbkg,  phacuts=phacuts, stdev=True)

    # Collect the BACKSCAL keywords
    backscal_src = float(fkeypar(fitsfile=fsrc, keyword='BACKSCAL', index=1))
    backscal_bkg = float(fkeypar(fitsfile=fbkg, keyword='BACKSCAL', index=1))

    # Background subtracted count rate (and standard deviation)
    bkgsub_cntrate = cntrate_src - (cntrate_bkg * (backscal_src / backscal_bkg))  # [counts / second]
    bkgsub_stdrate = stdrate_src - (stdrate_bkg * (backscal_src / backscal_bkg))  # [counts / second]
    if (stdev is True): return bkgsub_cntrate, bkgsub_stdrate
    if (stdev is False): return bkgsub_cntrate

#-----------------------------------------------------------------------------#
def optimize_region(cr_min, cr_max):
    """
    !!! UNDER CONSTRUCTION !!!
    Iteratively draw regions until the desired count rate is obtained

    count_rate : Desired count rate in the region (counts/sec)
    MAKE SURE YOU ARE CONSIDERING THE COUNT RATE IN THE RELEVANT ENERGY RANGE 0.5-10 keV
    """
    # Extract an image filtered by energy and count rate (check that extracting an image and spectrum gives same count rate)
    
    return

#=============================================================================#
