from Swift_XRT_WTmode import *

'''
PURPOSE:
--------
Retrieve some basic information about a Swift/XRT observation:

EXAMPLE:
--------
datadir = '/Users/salvesen/research/fcolabs/data/pileup/00030009021/'
fsrc    = datadir + "source_Rcut0_grd0-2_00030009021.pha"
fbkg    = datadir + "background_Rcut0_grd0-2_00030009021.pha"
phacuts = [50, 1000]
getinfo(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts)

NOTES:
------
- Must initialize HEAsoft
'''

def getinfo(fsrc, fbkg=None, phacuts=None):
    '''
    fsrc    - Spectrum file
    fbkg    - Background file
    phacuts - Range of photon energy channels to keep [phalcut, phahcut] (1 XRT channel = 10 eV)
    '''
    # Retrieve some keyword values from the FITS header
    source   = fkeypar(fitsfile=fsrc, keyword='OBJECT')
    obsID    = fkeypar(fitsfile=fsrc, keyword='OBS_ID')
    mjd      = fkeypar(fitsfile=fsrc, keyword='MJD-OBS')
    date     = fkeypar(fitsfile=fsrc, keyword='DATE-OBS')
    nevents  = fkeypar(fitsfile=fsrc, keyword='NEVENTS')  # <-- Not working?
    exposure = fkeypar(fitsfile=fsrc, keyword='EXPOSURE')
    tgood    = goodtime(fitsfile=fsrc)
    
    # Calculate the count rate and its standard deviation (background subtracted if fbkg is provided)
    if (fbkg is None): cntrate, stdrate = count_rate(fitsfile=fsrc, phacuts=phacuts, stdev=True)
    if (fbkg is not None): cntrate, stdrate = bkgsub_count_rate(fsrc=fsrc, fbkg=fbkg, phacuts=phacuts, stdev=True)
    
    # Off-axis angle [arcmin]
    deg2arcmin = 60.0  # [arcmin / degree]
    offaxis    = offaxis_angle(fitsfile=fsrc) * deg2arcmin

    print ""
    print "Source:                            ", source
    print "Observation ID:                    ", obsID
    print "Start Date [MJD; J2000]:           ", mjd
    print "Start Date [yyyy-mm-dd; hh:mm:ss]: ", date
    print "Number of Events:                  ", nevents
    print "Net Exposure Time [seconds]:       ", exposure
    print "Goodtime Interval [seconds]:       ", tgood
    print "Count Rate +/- Stdev [cnts/sec]:   ", cntrate, stdrate
    print "Off-Axis Angle [arcmin]:           ", offaxis
    print ""

