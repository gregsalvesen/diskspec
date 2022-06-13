from xspec_models import *
from plot_spectrum import *

'''
PURPOSE:
--------
Perform lots of spectral fits with a traditional chi-squared approach.
'''

# Swift observation ID number and its corresponding: clean data directory "clndir", results directory "outdir"
obsID  = '00030009021'
clndir = '/Users/salvesen/research/fcolabs/data/clean/' + obsID + '/'
outdir = '/Users/salvesen/research/fcolabs/results/' + obsID + '/'

# Data file
fdata = clndir + 'source_grp20_' + obsID + '.pha'

# Energy range and number of bins for the model
E_min = 0.5
E_max = 10.0
NE    = 1000

# GRO J1655-40 Parameters
nH   = 0.74   # Galactic column density [10^22 atoms cm^-2] (Dickey & Lockman 1990)
aFe  = 0.98   # Black hole spin parameter [-] (Reis et al. 2009)
aCF  = 0.70   # Black hole spin parameter [-] (Shafee et al. 2006)
iBin = 68.65  # Binary orbital inclination angle [degrees] (???)
iJet = 85.0   # Jet inclination angle [degrees] (???)
Mbh  = 5.4    # Black hole mass [M_sun] (???)
Dbh  = 3.2    # Distance [kpc] (???)
fcol = 1.7    # Color correction factor [-] (initial guess)

# bhspec specifics
Mbhspec = np.log10(Mbh)
cosiBin = np.cos(iBin * np.pi/180.0)
cosiJet = np.cos(iJet * np.pi/180.0)
Dbhspec = (10.0 / Dbh)**2

# Freeze a parameter by setting its delta < 0
freeze = -1

# Thaw a parameter by setting its delta > 0
def thaw(x, delta=0.0001):
    dx = np.abs(x * delta)  # Absolute value to avoid freezing parameters that can have negative values
    return dx

# Create an instance of the xspec_models class
xspec_models = xspec_models(fdata=fdata, E_min=E_min, E_max=E_max, NE=NE)

#----------------------------------------------------------------------------------------------------

# Parameter: [value, delta, min, bot, top, max]

#----------------------------------------------------------------------------------------------------
# Demonstrate that there is intrinsic absorption

# TBabs * diskbb  <-- w/ nH frozen
froot       = outdir + "TBabs_diskbb_nHfroz_" + obsID
vals_nH     = [nH,  freeze   ]
vals_Tin    = [1.0, thaw(1.0)]
vals_normbb = [1.0, thaw(1.0)]
xspec_models.TBabs_diskbb(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * diskbb
froot       = outdir + "TBabs_diskbb_" + obsID
vals_nH     = [nH,  thaw(nH) ]
vals_Tin    = [1.0, thaw(1.0)]
vals_normbb = [1.0, thaw(1.0)]
xspec_models.TBabs_diskbb(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

#----------------------------------------------------------------------------------------------------
# Demonstrate that a high-energy component is not necessary

# TBabs * (diskbb + powerlaw)
froot         = outdir + "TBabs_diskbb_powerlaw_" + obsID
vals_nH       = [ nH,  thaw(nH) ]
vals_Tin      = [ 1.0, thaw(1.0)]
vals_normbb   = [ 1.0, thaw(1.0)]
vals_PhoIndex = [-2.3, freeze   ]
vals_normpl   = [ 1.0, thaw(1.0)]
xspec_models.TBabs_diskbb_powerlaw(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb, \
    PhoIndex=vals_PhoIndex, normpl=vals_normpl)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * simpl * (diskbb)
froot         = outdir + "TBabs_diskbb_simpl_" + obsID
vals_nH       = [nH,   thaw(nH)  ]
vals_Tin      = [1.0,  thaw(1.0) ]
vals_normbb   = [1.0,  thaw(1.0) ]
vals_Gamma    = [2.3,  freeze    ]
vals_FracSctr = [0.05, thaw(0.05)]
vals_UpScOnly = [0.0,  freeze    ]
xspec_models.TBabs_diskbb_simpl(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb, \
    Gamma=vals_Gamma, FracSctr=vals_FracSctr, UpScOnly=vals_UpScOnly)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * (diskbb + compTT)
froot         = outdir + "TBabs_diskbb_compTT_" + obsID
vals_nH       = [nH,   thaw(nH)  ]
vals_Tin      = [1.0,  thaw(1.0) ]
vals_normbb   = [1.0,  thaw(1.0) ]
vals_Redshift = [0.0,  freeze    ]
vals_T0       = [1.0,  thaw(1.0) ]  # <-- Links to the diskbb Tin (see xspec_models.py)
vals_kT       = [20.0, freeze    ]
vals_taup     = [1.0,  thaw(1.0) ]
vals_approx   = [1.0,  freeze    ]
vals_normTT   = [1.0,  thaw(1.0) ]
xspec_models.TBabs_diskbb_compTT(froot=froot, nH=vals_nH, Tin=vals_Tin, normbb=vals_normbb, \
    Redshift=vals_Redshift, T0=vals_T0, kT=vals_kT, taup=vals_taup, approx=vals_approx, normTT=vals_normTT)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

#----------------------------------------------------------------------------------------------------
# Demonstrate that a physical disk model is a good fit, but how to choose a spin and inclination?

# TBabs * kerrbb
vals_nH     = [nH,   thaw(nH) ]
vals_eta    = [0.0,  freeze   ]
vals_Mbh    = [Mbh,  freeze   ]
vals_Mdd    = [1.0,  thaw(1.0)]
vals_Dbh    = [Dbh,  freeze   ]
vals_hd     = [1.7,  thaw(1.7)]
vals_rflag  = [1.0,  freeze   ]
vals_lflag  = [1.0,  freeze   ]
vals_normbb = [1.0,  freeze   ]  # <-- kerrbb norm is frozen if Mbh, Dbh, i are all frozen

# a = free; i = iBin
froot  = outdir + "TBabs_kerrbb_iBin_" + obsID
vals_a = [0.0,  thaw(1.0)]
vals_i = [iBin, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# a = free; i = iJet
froot  = outdir + "TBabs_kerrbb_iJet_" + obsID
vals_a = [0.0,  thaw(1.0)]
vals_i = [iJet, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# a = free; hd = 1.7; i = iBin
froot   = outdir + "TBabs_kerrbb_hd_iBin_" + obsID
vals_a  = [0.0,  thaw(1.0)]
vals_i  = [iBin, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=[1.7, freeze], rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# a = free; hd = 1.7; i = iJet
froot   = outdir + "TBabs_kerrbb_hd_iJet_" + obsID
vals_a  = [0.0,  thaw(1.0)]
vals_i  = [iJet, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=[1.7, freeze], rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# a = aFe; i = iBin
froot  = outdir + "TBabs_kerrbb_aFe_iBin_" + obsID
vals_a = [aFe,  freeze]
vals_i = [iBin, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * kerrbb
# a = aFe; i = iJet
froot  = outdir + "TBabs_kerrbb_aFe_iJet_" + obsID
vals_a = [aFe,  freeze]
vals_i = [iJet, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * kerrbb
# a = aCF; i = iBin
froot  = outdir + "TBabs_kerrbb_aCF_iBin_" + obsID
vals_a = [aCF,  freeze]
vals_i = [iBin, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * kerrbb
# a = aCF; i = iJet
froot  = outdir + "TBabs_kerrbb_aCF_iJet_" + obsID
vals_a = [aCF,  freeze]
vals_i = [iJet, freeze]
xspec_models.TBabs_kerrbb(froot=froot, nH=vals_nH, eta=vals_eta, a=vals_a, i=vals_i, Mbh=vals_Mbh, \
    Mdd=vals_Mdd, Dbh=vals_Dbh, hd=vals_hd, rflag=vals_rflag, lflag=vals_lflag, normbb=vals_normbb)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

#----------------------------------------------------------------------------------------------------
# Explore other disk models as a way of getting at systematic errors

# TBabs * diskpn
froot       = outdir + "TBabs_diskpn_" + obsID
vals_nH     = [nH,  thaw(nH) ]
vals_T_max  = [1.0, thaw(1.0)]
vals_R_in   = [6.0, freeze   ]
vals_normpn = [1.0, thaw(1.0)]
xspec_models.TBabs_diskpn(froot=froot, nH=vals_nH, T_max=vals_T_max, R_in=vals_R_in, normpn=vals_normpn)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * ezdiskbb
froot       = outdir + "TBabs_ezdiskbb_" + obsID
vals_nH     = [nH,  thaw(nH) ]
vals_T_max  = [1.0, thaw(1.0)]
vals_normez = [1.0, thaw(1.0)]
xspec_models.TBabs_ezdiskbb(froot=froot, nH=vals_nH, T_max=vals_T_max, normez=vals_normez)
plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

# TBabs * bhspec
froot          = outdir + "TBabs_bhspec_" + obsID
vals_log_mass  = [Mbhspec, freeze    ]
vals_log_lumin = [-1.0,    thaw(-1.0)]
vals_inc       = [cosiBin, freeze    ]
vals_spin      = [aFe,     freeze    ]
vals_normbh    = [Dbhspec, freeze    ]
#xspec_models.TBabs_bhspec(froot=froot, nH=vals_nH, log_mass=vals_log_mass, log_lumin=vals_log_lumin, \
#    inc=vals_inc, spin=vals_spin, normbh=vals_normbh)
#plot_spectrum(fspec=froot+'.hdf5').plot_bestfit(fout=froot+'_spectrum.png')

