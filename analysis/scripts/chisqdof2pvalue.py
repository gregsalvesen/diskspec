import numpy as np
import decimal
import scipy.integrate as integrate

'''
PURPOSE:
--------
Given a chi^2 value and number of degrees of freedom, calculate the null hypothesis p-value.

NOTES:
------
This routine cannot handle a large number of degrees of freedom, which is a problem.
Look into scipy.special.gamma()
Hmm, I don't think the integral in Gamma2() below can be done with a standard gamma function?

REFERENCES:
-----------
http://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
'''
def calc_pvalue(chisq, dof):

    half    = decimal.Decimal(0.5)
    one     = decimal.Decimal(1.0)
    two     = decimal.Decimal(2.0)
    halfdof = decimal.Decimal(0.5 * dof)

    def G1Igrand(t, x):
        x      = decimal.Decimal(x)
        t      = decimal.Decimal(t)
        Igrand = t**(x - one) * np.exp(-t)
        return Igrand
    
    def Gamma1(x, lolim, uplim):
        G1 = integrate.quad(func=G1Igrand, a=lolim, b=uplim, args=(x,))
        return G1
    
    # Factor
    G1Tuple = Gamma1(x=halfdof, lolim=0.0, uplim=np.inf)
    G1      = decimal.Decimal(G1Tuple[0])
    factor  = (two**(halfdof) * G1)**(-one)

    def G2Igrand(t, x):
        x      = decimal.Decimal(x)
        t      = decimal.Decimal(t)
        Igrand = t**(x - one) * np.exp(-t * half)
        return Igrand
    
    def Gamma2(x, lolim, uplim):
        G2 = integrate.quad(func=G2Igrand, a=lolim, b=uplim, args=(x,))
        return G2

    # p-value
    G2Tuple = Gamma2(x=halfdof, lolim=chisq, uplim=np.inf)
    G2      = decimal.Decimal(G2Tuple[0])
    pvalue  = factor * G2
    return pvalue
