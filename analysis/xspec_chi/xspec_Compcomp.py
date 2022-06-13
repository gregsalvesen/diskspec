import numpy as np
from xspec import *
from pvalue2sigma import *

'''
PURPOSE:
--------
Show that for a high-energy Comptonizing component (powerlaw, simpl, compTT) does not statistically improve the TBabs*diskbb fit. Do this by calculating the null hypothesis p-value (and corresponding sigma-value) resulting from an F-test.
'''

rdir   = "/Users/salvesen/research/fcolabs/results/"
obsIDs = ["00030009021", "00030009022", "00030009023", "00030009025", "00030009026"]

# Return the degrees of freedom from an output .fit file
def return_dof(ffit):
    with open(ffit, "r") as f:
        for line in f:
            if line.startswith("Degrees of Freedom"):
                dof = int(line.split(" ")[-1])
                return dof

# Return the chi-squared value from an output .fit file
def return_chisq(ffit):
    with open(ffit, "r") as f:
        for line in f:
            if line.startswith("Chi-Squared"):
                chisq = float(line.split(" ")[-1])
                return chisq

# Loop through each observation ID
for obsID in obsIDs:

    # Files containing the fit results
    f_diskbb = rdir + obsID + "/TBabs_diskbb_"          + obsID + ".fit"
    f_powlaw = rdir + obsID + "/TBabs_diskbb_powerlaw_" + obsID + ".fit"
    f_simpl  = rdir + obsID + "/TBabs_diskbb_simpl_"    + obsID + ".fit"
    f_compTT = rdir + obsID + "/TBabs_diskbb_compTT_"   + obsID + ".fit"

    # Degrees of Freedom
    dof_diskbb = return_dof(ffit=f_diskbb)
    dof_powlaw = return_dof(ffit=f_powlaw)
    dof_simpl  = return_dof(ffit=f_simpl)
    dof_compTT = return_dof(ffit=f_compTT)

    # Chi-squared Statistic
    chisq_diskbb = return_chisq(ffit=f_diskbb)
    chisq_powlaw = return_chisq(ffit=f_powlaw)
    chisq_simpl  = return_chisq(ffit=f_simpl)
    chisq_compTT = return_chisq(ffit=f_compTT)

    # Delta Chi-squared
    Dchisq_powlaw = chisq_powlaw - chisq_diskbb
    Dchisq_simpl  = chisq_simpl  - chisq_diskbb
    Dchisq_compTT = chisq_compTT - chisq_diskbb

    # F-test p-value
    pvalue_powlaw = Fit.ftest(chisq2=chisq_powlaw, dof2=dof_powlaw, chisq1=chisq_diskbb, dof1=dof_diskbb)
    pvalue_simpl  = Fit.ftest(chisq2=chisq_simpl,  dof2=dof_simpl,  chisq1=chisq_diskbb, dof1=dof_diskbb)
    pvalue_compTT = Fit.ftest(chisq2=chisq_compTT, dof2=dof_compTT, chisq1=chisq_diskbb, dof1=dof_diskbb)

    # F-test sigma-value
    #sigma_powlaw = sigma_fstat(pvalue=pvalue_powlaw, dfn=dof_powlaw, dfd=dof_diskbb)
    #sigma_simpl  = sigma_fstat(pvalue=pvalue_simpl,  dfn=dof_simpl,  dfd=dof_diskbb)
    #sigma_compTT = sigma_fstat(pvalue=pvalue_compTT, dfn=dof_compTT, dfd=dof_diskbb)

    # Print out the results
    print "Observation ID: ", obsID
    print "        powlaw: ", Dchisq_powlaw, pvalue_powlaw#, sigma_powlaw
    print "         simpl: ", Dchisq_simpl,  pvalue_simpl#,  sigma_simpl
    print "        comptt: ", Dchisq_compTT, pvalue_compTT#, sigma_compTT
    print ""
