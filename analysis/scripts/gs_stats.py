import numpy as np
from scipy.integrate import quad, simps
from scipy.special import erf, erfc
from scipy.optimize import curve_fit

'''
This is a collection of handy statistics tools:

interpolated_median -
confidence_interval -
error1D -

'''
#================================================================================
def interpolated_median(x, y):
    '''
        Find the (x,y) point corresponding to the interpolated median of the distribution
        '''
    N = len(x)
    area = 0.0
    area50 = 0.5
    i = 0
    while (area < area50):
        area = np.trapz(y=y[0:i+1], x=x[0:i+1])
        i = i+1
    
    i = i-1  # <-- Account for +1 on last iteration
    iL, iR = i, i-1
    xR, yR = x[iR], y[iR]
    xL, yL = x[iL], y[iL]
    areaR  = np.trapz(y=y[0:iR+1], x=x[0:iR+1])
    areaL  = np.trapz(y=y[0:iL+1], x=x[0:iL+1])
    
    # Interpolation factor
    f = (area50 - areaL) / (areaR - areaL)
    
    x_ML = f * (xR - xL) + xL
    y_ML = f * (yR - yL) + yL
    
    return x_ML, y_ML


#================================================================================
def confidence_interval(x, y, sigma=0.682689492137, thold=0.001):
    
    x_ML, y_ML = interpolated_median(x=x, y=y)
    
    # Split the x-array into left and right of the x-median
    iL_arr = np.where(x <= x_ML)[0]
    iR_arr = np.where(x > x_ML)[0]
    xL_arr = x[iL_arr]
    xR_arr = x[iR_arr]
    yL_arr = y[iL_arr]
    yR_arr = y[iR_arr]
    
    NiL = len(iL_arr)
    NiR = len(iR_arr)
    
    # Start at y = 0.5*y_ML and area = 0
    y_now    = 0.5 * y_ML
    area_now = 0.0
    y_top    = y_ML
    y_bot    = 0.0
    # Iterate until the enclosed area is within the desired threshold
    while (np.abs(area_now - sigma) > thold):
        
        # Find the xL value corresponding to y_now
        yL = y_ML
        iL = NiL-1
        while (yL > y_now):
            yL = yL_arr[iL]
            iL = iL-1
        iL = iL+1
        xL = xL_arr[iL]
        yL = yL_arr[iL]
        iLp1 = iL+1
        xLp1 = xL_arr[iLp1]
        yLp1 = yL_arr[iLp1]
        # Interpolate to get xL_now
        fL_now = (y_now - yL) / (yLp1 - yL)
        xL_now = fL_now * (xLp1 - xL) + xL
        
        # Find the xR value corresponding to y_now
        yR = y_ML
        iR = 0
        while (yR > y_now):
            yR = yR_arr[iR]
            iR = iR+1
        iR = iR-1
        xR = xR_arr[iR]
        yR = yR_arr[iR]
        iRm1 = iR-1
        xRm1 = xR_arr[iRm1]
        yRm1 = yR_arr[iRm1]
        # Interpolate to get xR_now
        fR_now = (y_now - yRm1) / (yR - yRm1)
        xR_now = fR_now * (xR - xRm1) + xRm1
        
        # Integrate the distribution over xL_now --> xR_now
        # (Cut off the area beyond the interpolation)
        area = np.trapz(y=y[iL:NiL+iR+1], x=x[iL:NiL+iR+1])
        # How much area are we missing to the left [xL --> x_now]?
        area_xtraL = np.trapz(y=[yL,y_now], x=[xL,xL_now])
        # How much area are we missing to the right [x_now --> xR]?
        area_xtraR = np.trapz(y=[y_now,yR], x=[xR_now,xR])
        # Find the proper enclosed area
        area_now = area - (area_xtraL + area_xtraR)
        
        # Does area_now = the desired sigma to within our threshold?
        if (area_now > sigma):  # Raise y_now to contain smaller area
            if (y_now > y_bot): y_bot = y_now
            y_now = 0.5 * (y_now + y_top)
        if (area_now < sigma):  # Lower y_now to contain larger area
            if (y_now < y_top): y_top = y_now
            y_now = 0.5 * (y_now + y_bot)
    
    x_valm = xL_now
    x_valp = xR_now
    y_valm = y_now
    y_valp = y_now
    
    return x_ML, y_ML, x_valm, y_valm, x_valp, y_valp


#================================================================================
def error1D(x, y, sigma=0.682689492137, normal=False):
    '''
    Purpose:
    --------
    Compute the 1D (possibly asymmetric) errorbar at constant probability 
    value on the input PDF.
        
    Inputs:
    -------
    x     - Bins corresponding to the PDF
    y     - Probability density function (PDF)
    sigma - Integrate out in the PDF until a fraction "sigma" of the likelihood is enclosed.
    
    Output:
    -------
    xpeak, errm, errp = [PDF peak, -error, +error]
    '''
    # Index of the maximum likelihood value.
    iML = np.argmax(y)
    
    # Beginning and ending indices.
    iB = 0
    iE = y.size-1
    
    # Initialize the starting index.
    iR = iML+1
    if (iR > iE):
        iR = iE
    
    # Total enclosed area.
    #area_tot = simps(y, x=x) <-- GIVES BAD RESULT!
    area_tot = np.trapz(y=y, x=x)
    
    # Initialize the enclosed area.
    area = 0.0
    
    # Compute the area under the PDF until sigma-% is achieved. <-- This is very slow, can we do it smarter???
    xL = np.interp(y[iR], y[0:iML-1], x[0:iML-1])
    iL = int(np.argmin(np.abs(xL - np.array(x))))
    while (area < sigma):
        #xL = np.interp(y[iR], y[0:iML-1], x[0:iML-1])
        #iL = int(np.argmin(np.abs(xL - np.array(x))))

        # Handle very skewed distributions
        if (iL > iR):
            iL = iB
            iR = iE
        if ((iL != iB) and (iR != iE)):
            #area = simps(y[iL:iR+1], x=x[iL:iR+1]) / area_tot <-- GIVES BAD RESULT!
            area = np.trapz(y[iL:iR+1], x=x[iL:iR+1]) / area_tot
            iR += 1
            xL = np.interp(y[iR], y[0:iML-1], x[0:iML-1])
            iL = int(np.argmin(np.abs(xL - np.array(x))))
        if ((iL == iB) and (iR != iE)):
            area = np.trapz(y[iB:iR+1], x=x[iB:iR+1]) / area_tot
            iR += 1
            xL = np.interp(y[iR], y[0:iML-1], x[0:iML-1])
            iL = int(np.argmin(np.abs(xL - np.array(x))))
        if ((iL != iB) and (iR == iE)):
            area = np.trapz(y[iL:iE+1], x=x[iL:iE+1]) / area_tot
            xL = np.interp(y[iE], y[0:iML-1], x[0:iML-1])
            iL = int(np.argmin(np.abs(xL - np.array(x))))
            iML -= 1
            iL -= 1
        if ((iL == iB) and (iR == iE)):
            area = np.trapz(y[iB:iE+1], x=x[iB:iE+1]) / area_tot
        #else:
        #    xpeak, errm, errp = pdf_peakerr(x, y, sigma=sigma)  # Compute xpeak, errm, errp with pdf_peakerr
        #    area = sigma + 1.0  # Pop out of the loop.
        #iR += 1
        
    # Account for the extra iR += 1 performed on the last while loop iteration.
    iR -= 1
    
    # Collect the results
    #if ((iL != iB) and (iR != iE)):
    #    xpeak = x[iML]
    #    errm  = x[iML] - x[iL]
    #    errp  = x[iR] - x[iML]
    
    xpeak = x[iML]
    errm  = x[iML] - x[iL]
    errp  = x[iR] - x[iML]
    
    return xpeak, errm, errp


#===============================================================================
def pdf_peakerr(x, y, sigma=0.682689492137):
    '''
    Purpose:
    --------
    Find the peak of a probability distribution function and compute the 1D
    (asymmetric) errorbar.
            
    Inputs:
    ------
    x     - Bins corresponding to the PDF
    y     - Probability distribution function (PDF)
    sigma - Integrate out in the PDF until a fraction "sigma" of the
            likelihood is enclosed (default is 1-sigma = 68%).
        
    Outputs:
    -------
    xpeak, errm, errp = [PDF peak, -error, +error]
    '''
        
    # Index of the peak (i.e., the maximum likelihood index).
    ipeak = np.argmax(y)
    
    # Value corresponding to the peak of the PDF.
    xpeak = x[ipeak]
        
    # Beginning and ending indices.
    iB = 0
    iE = y.size-1
        
    # Initialize the indices that will be rightward and leftward of the peak.
    iR = ipeak
    iL = ipeak
        
    # Total area contained in the PDF.
    #area_tot = simps(y, x=x) <-- GIVES BAD RESULT!
    area_tot = np.trapz(y, x=x)
    
    # Initialize the enclosed area to the peak of the distribution.
    #area = 0.0
    dx = x[1] - x[0]
    area = y[ipeak] * dx
    
    # Step to the left and right incrementally until you reach the desired error bar or run out of bins. <-- Use a more sophisticated integrator, this just sums stuff
    while ((area / area_tot) < sigma):
        if (iR != iE):
            if (np.sum(y[iR::]) != 0.0):
                iR += 1
                area += y[iR] * dx
        if (iL != iB):
            if (np.sum(y[0:iL+1]) != 0.0):
                iL -= 1
                area += y[iL] * dx
        if ((iL == iB) and (iR == iE)):  # Escape this awful loop that needs looking into...
            area = area_tot
        
    errm = x[ipeak] - x[iL]
    errp = x[iR] - x[ipeak]

    return xpeak, errm, errp


#===============================================================================
def distfunc(y, Nbins, xmin, xmax):
    '''
    Purpose:
    --------
    Compute the binned distribution from an unbinned distribution determined
    from a Monte Carlo simulation.
            
    Inputs:
    ------
    y     - Unbinned distribution from a Monte Carlo simulation (array)
    Nbins - Number of bins
    xmin  - Minimum value for the histogram bin range
    xmax  - Maximum value for the histogram bin range
            
    Output:
    -------
    Binned distribution, bin center values
    '''
    bincenters = []
    dist, binedges = np.histogram(y, bins=Nbins, range=(xmin,xmax))
    for i in xrange(Nbins):
        bincenters.append(0.5 * (binedges[i] + binedges[i+1]))
    return dist, bincenters


#===============================================================================            
def histbound(x, x_peak, errm, errp):
    '''
    Purpose:
    --------
    Given a distribution, remove the parts that do not fall within the
    specified +/- range.
            
    Inputs:
    ------
    x      - Unbinned distribution from a Monte Carlo simulation (array)
    x_peak - Peak as determined in the binned distribution
    errm   - Minus error (lower bound = x_peak - errm)
    errp   - Plus error (upper bound = x_peak + errp)
            
    Output:
    -------
    Input distribution with values lying outside of the given bounds
    filtered out.
    '''
    x_errm = x_peak - errm
    x_errp = x_peak + errp
    iin    = range(x.size)                     # Initialize list of bounded ("in") indices
    iout_m = np.where(x < x_errm)[0].tolist()  # Leftward out of bounds ("out") indices
    iout_p = np.where(x > x_errp)[0].tolist()  # Rightward out of bounds ("out") indices
    iout   = list(set(iout_m + iout_p))        # List of out of bounds indices
    for i in xrange(np.array(iout).size):
        iin.remove(iout[i])                    # Remove "out" indices from list of "in" indices
    return x[iin]

    
#===============================================================================
def sigma_search(x_peak, err, prob, x_min, x_max, A=1.0, acc=1e-2):
    '''
    Purpose:
    --------
    Given a confidence interval that is not 1-sigma (e.g. 90%), find the
    corresponding 1-sigma confidence interval assuming the errors are
    Gaussian distributed.
            
    Inputs:
    -------
    x_peak - x-value corresponding to the peak of the Gaussian.
             (This is the data point measured value in practice)
    err    - Error bar on x_peak, the measured data point
    prob   - Confidence associated with err (e.g. 0.90 for 90% confidence)
    x_min  - Minimum allowable value for the data point (e.g. +1 for black hole spin)
    x_max  - Maximum allowable value for the data point (e.g. -1 for black hole spin)
    A      - Amplitude of the peak of the Gaussian (i.e. where CDF = 0.5)
    acc    - Absolute (not relative) accuracy to which the 1-sigma value is computed
             (also used to compute the number of x-bins to use)
            
    Output:
    -------
    Standard deviation (sigma) corresponding to the Gaussian distribution that
    encloses the confidence region 'prob' within the x-range, (x_peak +/- err).
    '''
        
    sigma_68 = 0.682689492137  # Confidence interval corresponding to 1-sigma
    sigma_guess = err          # Initial guess for the confidence interval
    x_err = x_peak + err       # x-location of the error bar extent

    # Compute the magic factor for converting between sigma and prob
    dprob = acc + 1.0
    while (np.abs(dprob) > acc):
        # Compute the CDF probability given the guessed sigma
        prob_guess = 0.5 * (1.0 + erf((x_err - x_peak) / (2.0 * sigma_guess**2.0)))
        # Modify sigma_guess based on if prob_guess is > or < prob
        dprob = prob_guess - prob
        if (prob_guess >= prob):
            dsigma_guess = 0.5 * dprob # Need to increase sigma_guess
        else:
            dsigma_guess = 0.5 * dprob # Need to decrease sigma_guess
        # Increment the sigma guess
        sigma_guess += dsigma_guess

    # This technique below should work too, but it takes forever.
    '''
    x_err = x_peak - err  # x-location of the error bar extent
    dx = acc + 1.0        # Initialize the dx value to be > the accuracy requested.
        
    # Make *sure* there are adequate bins to resolve the 1-sigma value.
    Nbins = int(10.0 * np.abs(x_max - x_min) / acc)

    # Iteratively compute Gaussian PDFs until converging on the proper sigma.
    x = np.linspace(x_min, x_max, Nbins)
    while (np.abs(dx) > acc):
        # Compute the PDF
        pdf = Gauss1D_pdf(x, A=A, mu=x_peak, sigma=sigma_guess)
        
        # Given the PDF just constructed, find the x-values that enclose a
        # fraction 'prob' of the distribution.
        x_peak_guess, errm_guess, errp_guess = error1D(x, pdf, sigma=prob)
        
        # +/- value from the peak of the distribution where a fraction 'prob' of the distribution is enclosed.
        x_err_guess = x_peak - errm_guess  # x-location
        if (np.allclose(errm_guess, errp_guess, rtol=acc, atol=acc) == False):
            print "ERRORS ARE ASYMMETRIC -- THIS IS NOT A GAUSSIAN!"
            sys.exit(1)
        
        # Modify sigma_guess based on if x_prob is > or < x_err
        dx = x_err_guess - x_err
        x_sigma_guess = x_peak - sigma_guess
        if (x_err_guess >= x_err):
            dsigma_guess = 0.5 * (x_err_guess - x_err) # Need to increase sigma_guess
        else:
            dsigma_guess = 0.5 * (x_err_guess - x_err) # Need to decrease sigma_guess
        sigma_guess += dsigma_guess
    '''
    
    # Account for the extra +dsigma made on the final iteration.
    sigma_guess -= dsigma_guess

    return sigma_guess





#================================================================================
def dblGauss_fit(xdata, ydata):
    '''
    Purpose:
    --------
    
    '''
    # Define the Gaussian fitting function.
    def func_Gauss(x, mu, sigma, A):
        return A * np.exp(-1.0*(x-mu)**2.0 / (2.0*sigma**2.0))
    
    # Define the double matched Gaussian fitting function.
    def func_dblGauss(x, mu, sigmaL, sigmaR, A):
        pdf = np.zeros(len(x))
        i_peak = int((np.abs(x-mu)).argmin())
        pdf[0:i_peak] = A * np.exp(-1.0*(x[0:i_peak]-mu)**2.0 / (2.0*sigmaL**2.0))
        pdf[i_peak::] = A * np.exp(-1.0*(x[i_peak::]-mu)**2.0 / (2.0*sigmaR**2.0))
        return pdf
    
    # Generate initial parameter guesses by fitting the distribution with a single Gaussian.
    pGauss_opt, pGauss_cov = curve_fit(func_Gauss, xdata, ydata)
    p0_guess = [pGauss_opt[0], pGauss_opt[1], pGauss_opt[1], pGauss_opt[2]]
    
    # Fit the distribution with a double matched Gaussian.
    p_opt, p_cov = curve_fit(func_dblGauss, xdata, ydata, p0=p0_guess)
    x_mu = p_opt[0]
    
    # Compute the best-fit double matched Gaussian distribution
    y_bestfit = func_dblGauss(xdata, p_opt[0], p_opt[1], p_opt[2], p_opt[3])
    
    return x_mu

