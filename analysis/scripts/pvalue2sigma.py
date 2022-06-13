import numpy as np
import scipy.stats

# Convert a null hypothesis p-value into a sigma value.
# http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html

#----------------------------------------------------------------------------------------------------
# chi-square statistic distribution

# Calculate the p-value
def pvalue_chisq(chisq, dof):
    pvalue = scipy.stats.chi2.sf(x=chisq, df=dof)  # <-- This is just (1 - cdf)
    return pvalue

# Calculate the sigma-value
def sigma_chisq(pvalue, df=1):
    confint = 1.0 - pvalue
    sigma   = np.sqrt(scipy.stats.chi2.ppf(confint, df=df))
    return sigma

#----------------------------------------------------------------------------------------------------
# F-statistic distribution

# Calculate the sigma-value <-- IS THIS RIGHT?
def sigma_fstat(pvalue, dfn, dfd):
    confint = 1.0 - pvalue
    sigma = np.sqrt(scipy.stats.f.ppf(confint, dfn=dfn, dfd=dfd))
    return sigma
