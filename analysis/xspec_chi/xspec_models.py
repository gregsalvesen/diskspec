import numpy as np
from xspec_modelfit import *

'''
'''

#====================================================================================================
class xspec_models:

    #----------------------------------------------------------------------------------------------------
    # Initialization
    def __init__(self, fdata, E_min=0.5, E_max=10.0, NE=100):
        self.fdata = fdata
        self.E_min = E_min
        self.E_max = E_max
        self.NE    = NE

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*diskbb
    # ######
    def TBabs_diskbb(self, froot, nH, Tin, normbb):
        model = "TBabs*diskbb"
        setup = [('TBabs',  {'nH':nH}), \
                 ('diskbb', {'Tin':Tin, 'norm':normbb})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*(diskkbb+powerlaw)
    # ######
    def TBabs_diskbb_powerlaw(self, froot, nH, Tin, normbb, PhoIndex, normpl):
        model = "TBabs*(diskbb+powerlaw)"
        setup = [('TBabs',    {'nH':nH}), \
                 ('diskbb',   {'Tin':Tin, 'norm':normbb}), \
                 ('powerlaw', {'PhoIndex':PhoIndex, 'norm':normpl})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)
    
    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*simpl*(diskkbb)
    # ######
    def TBabs_diskbb_simpl(self, froot, nH, Tin, normbb, Gamma, FracSctr, UpScOnly):
        model = "TBabs*simpl*(diskbb)"
        setup = [('TBabs',  {'nH':nH}), \
                 ('diskbb', {'Tin':Tin, 'norm':normbb}), \
                 ('simpl',  {'Gamma':Gamma, 'FracSctr':FracSctr, 'UpScOnly':UpScOnly})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*(diskkbb+compTT)
    # ######
    def TBabs_diskbb_compTT(self, froot, nH, Tin, normbb, Redshift, T0, kT, taup, approx, normTT):
        model = "TBabs*(diskbb+compTT)"
        setup = [('TBabs',  {'nH':nH}), \
                 ('diskbb', {'Tin':Tin, 'norm':normbb}), \
                 ('compTT', {'Redshift':Redshift, 'T0':T0, 'kT':kT, 'taup':taup, 'approx':approx, 'norm':normTT})]
        link  = [('compTT', 'T0', 'diskbb', 'Tin')]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata, link=link)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*kerrbb
    # ######
    def TBabs_kerrbb(self, froot, nH, eta, a, i, Mbh, Mdd, Dbh, hd, rflag, lflag, normbb):
        model = "TBabs*kerrbb"
        setup = [('TBabs',  {'nH':nH}), \
                 ('kerrbb', {'eta':eta, 'a':a, 'i':i, 'Mbh':Mbh, 'Mdd':Mdd, 'Dbh':Dbh, 'hd':hd, \
                             'rflag':rflag, 'lflag':lflag, 'norm':normbb})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*diskpn
    # ######
    def TBabs_diskpn(self, froot, nH, T_max, R_in, normpn):
        model = "TBabs*diskpn"
        setup = [('TBabs',  {'nH':nH}), \
                 ('diskpn', {'T_max':T_max, 'R_in':R_in, 'norm':normpn})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*ezdiskbb
    # ######
    def TBabs_ezdiskbb(self, froot, nH, T_max, normez):
        model = "TBabs*ezdiskbb"
        setup = [('TBabs',    {'nH':nH}), \
                 ('ezdiskbb', {'T_max':T_max, 'norm':normez})]
        xspec_modelfit(theModel=model, setup=setup, froot=froot, fdata=self.fdata)

    #----------------------------------------------------------------------------------------------------
    # ######
    # MODEL: TBabs*bhspec
    # ######

    '''
    # BHSPEC specifics <-- might be worth treating this completely separate?
    log_mass  = np.log10(Mbh)
    cos_iBin  = np.cos(iBin * np.pi/180.0)
    Dbhspec   = (10.0 / Dbh)**2
    # Load in the (appropriate) BHSPEC model with default parameters
    if (alpha == 0.01): fbhspec = "bhspec_spin_0.01.fits"
    if (alpha == 0.1):
        if (a < 0.0):  fbhspec = "bhspec_spin_0.1.fits"
        if (a >= 0.0): fbhspec = "bhspec_spin2_0.1.fits"
    Model("atable{"+dir+fbhspec+"}")
    
    val_log_mass,  min_log_mass,  max_log_mass  = log_mass,  0.477, 1.477   # Log10(black hole mass) [Msun]
    val_log_lumin, min_log_lumin, max_log_lumin = -1.0,     -2.0,   0.0     # Log10(disk luminosity) [LEdd]
    val_inc,       min_inc,       max_inc       = cos_iBin,  0.0,   1.0     # Cos(disk inclination) [-]
    val_spin,      min_spin,      max_spin      = aFeK,     -1.0,   np.nan  # Black hole spin parameter [-] <-- Max value?
    val_norm,      min_norm,      max_norm      = Dbhspec,   0.0,   1.0e24  # "Normalization", but really the distance [???]
    '''
