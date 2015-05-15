#!/usr/bin/env python

import kepmsg
import numpy, scipy, math, random
from math import *
from scipy import stats, linalg
from scipy.linalg import pinv
from numpy import * 

# -----------------------------------------------------------
# calculate sum of array

def sum(a):
    return a.sum()

# -----------------------------------------------------------
# calculate sum of array

def sumerr(a):

    work = 0.0
    for item in a:
	work += item**2
    err = sqrt(work)
    return err

# -----------------------------------------------------------
# calculate mean of numeric list

def mean(list):

    try:
        mean = scipy.stats.nanmean(list)
    except:
        total = 0.0
        for item in list:
            total += item
        mean = total / len(list)
    return mean

# -----------------------------------------------------------
# calculate error on mean of numeric list

def mean_err(list):

    total = 0.0
    for item in list:
	total = total + item**2
    err = sqrt(total) / len(list)
    return err

# -----------------------------------------------------------
# calculate median of numeric list

def median(list,logfile):

    list.sort()
    n = len(list)
    if (n == 0):
	message = 'ERROR -- KEPSTAT.MEDIAN: Supplied list is empty'
	status = kepmsg.err(logfile,message)
	median = None
    elif (n < 3):
	median = mean(list)
    else:
	median = list[n/2]
    return median
	
# -----------------------------------------------------------
# minimum of array

def min(array):

    minm = array[0]
    for i in range(1,len(array)):
	if (array[i] < minm): minm = array[i]
    return minm

# -----------------------------------------------------------
# maximum of array

def max(array):

    maxm = array[0]
    for i in range(1,len(array)):
	if (array[i] > maxm): maxm = array[i]
    return maxm

# -----------------------------------------------------------
# minimum of array with error

def mine(array,error):

    minm = array[0] - error[0]
    for i in range(1,len(array)):
	if (array[i] - error[i] < minm): minm = array[i] - error[i]
    return minm

# -----------------------------------------------------------
# maximum of array with error

def maxe(array,error):

    maxm = array[0] + error[0]
    for i in range(1,len(array)):
	if (array[i] + error[i] > maxm): maxm = array[i] + error[i]
    return maxm

# -----------------------------------------------------------
# root mean square of two arrays

def rms(array1,array2,logfile,verbose):

    sigma = 0
    status = 0
    if (len(array1) != len(array2)):
	message  = 'ERROR -- KEPSTAT.RMS: Arrays have unequal sizes - '
	message += 'Array1 = ' + str(len(array1)) + ', array2 = ' + str(len(array2))
	status = kepmsg.err(logfile,message,verbose)
    if (status == 0):
	for i in range(len(array1)):
	    sigma += (array1[i] - array2[i])**2
	sigma = math.sqrt(sigma / len(array1))
    return sigma, status

# -----------------------------------------------------------
# root mean square of two 2D arrays

def rms2d(array1,array2):

    sigma = 0
    n = 0
    array3 = (array1 - array2)**2
    for a in array3:
	sigma += a
	n += 1
    sigma /= n
    return sigma

# -----------------------------------------------------------
# standard deviation of an array

def stdev(array):

    sigma = 0.0
    average = mean(array)
    for i in range(len(array)):
	sigma += (array[i] - average)**2
    sigma = math.sqrt(sigma / len(array))
    return average, sigma

# -----------------------------------------------------------
# adjust data relative to random number constrained by error bars

def randarray(signal,err):

    random.seed()
    out = numpy.zeros(len(signal),dtype='float32')
    for i in range(len(signal)):
        out[i] = signal[i] + err[i] * inv_normal_cummulative_function(random.random())
    return out

# -----------------------------------------------------------
# remove infinities from light curve data

def removeinfinlc(x, cols):

    for j in range(len(cols)):
        work = []
        datatype = cols[j].dtype
	for i in range(len(x)):
            if numpy.isfinite(x[i]):
                work.append(cols[j][i])
        cols[j] = numpy.array(work,dtype=datatype)
    return cols

# -----------------------------------------------------------
# filter on data within time ranges

def filterOnRange(intime,tstart,tstop):
    
    status = 0
    outlist = []
    for i in range(len(intime)):
        for j in range(len(tstart)):
            if intime[i] > tstart[j] and intime[i] < tstop[j]:
                if len(outlist) == 0:
                    outlist.append(i)
                elif i > outlist[-1]:
                    outlist.append(i)
                    
    return outlist, status

#------------------------------------------------------------------
# inverse normal cummulative function

def inv_normal_cummulative_function(p):

# Lower tail quantile for standard normal distribution function.
#
# This function returns an apprximation of the inverse cumulative
# standard normal distribution function.  i.e., given P, it returns
# an approximation to the X satisfying P = Pr{Z <= X} where Z is a
# random variable from the standard normal distribution.
#
# The algorithm uses a minimax approximation by rational functions
# and the result has a relative error whose absolute value is less
# than 1.15e-9.
#
# Author:      Peter J. Acklam
# Time-stamp:  2000-07-19 18:26:14
# E-mail:      pjacklam@online.no
# WWW URL:     http://home.online.no/~pjacklam

    if p == 0.0: p = 0.0000000001
    if p == 1.0: p = 0.9999999999

# coefficients in rational approximations

    a = [-3.969683028665376e1,  2.209460984245205e2,
          -2.759285104469687e2,  1.383577518672690e2,
          -3.066479806614716e1,  2.506628277459239]
    b = [-5.447609879822406e1,  1.615858368580409e2,
          -1.556989798598866e2,  6.680131188771972e1,
          -1.328068155288572e1]
    c = [-7.784894002430293e-3, -3.223964580411365e-1,
          -2.400758277161838, -2.549732539343734,
          4.374664141464968,  2.938163982698783]
    d = [7.784695709041462e-3,  3.224671290700398e-1,
         2.445134137142996,  3.754408661907416]

# define break-points

    plow = 0.02425
    phigh = 1.0 - plow

# rational approximation for lower region

    if p < plow:
        q = sqrt(-2.0 * log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)

# rational approximation for upper region

    if phigh < p:
        q = sqrt(-2.0 * log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)

# rational approximation for central region

    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0)

# -----------------------------------------------------------
# bit map decoding

def bitInBitmap(bitmap,bit):
    
    flag = False
    for i in range(10,-1,-1):
        if (bitmap - 2**i >= 0):
            bitmap = bitmap - 2**i
            if 2**i == bit:
                flag = True
        else:
            continue

    return flag

# -----------------------------------------------------------
# principal components

def princomp(A):

 """ performs principal components analysis 
     (PCA) on the n-by-p data matrix A
     Rows of A correspond to observations, columns to variables. 

 Returns :  
  coeff :
    is a p-by-p matrix, each column containing coefficients 
    for one principal component.
  score : 
    the principal component scores; that is, the representation 
    of A in the principal component space. Rows of SCORE 
    correspond to observations, columns to components.

  latent : 
    a vector containing the eigenvalues 
    of the covariance matrix of A.
 """

# computing eigenvalues and eigenvectors of covariance matrix

 M = (A - numpy.mean(A.T,axis=1)).T # subtract the mean (along columns)
 [latent,coeff] = numpy.linalg.eig(numpy.cov(M))
 score = numpy.dot(coeff.T,M) # projection of the data in the new space

 return coeff,score,latent

# -----------------------------------------------------------
def savitzky_golay(y,window_size,order,deriv=0):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """

    status = 0
    sg = 0.0

# test the input attributes

    try:
        window_size = abs(int(window_size))
        order = abs(int(order))
    except:
        message = 'ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size and order must be of type int'
        status = kepmsg.err(None,message,True)
    if window_size % 2 != 1 or window_size < 1:
        message = 'ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size size must be a positive odd number'
        status = kepmsg.err(None,message,True)
    if window_size < order + 2:
        message = 'ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size is too small for the polynomials order'
        status = kepmsg.err(None,message,True)

# precompute fit coefficients

    if status == 0:
        order_range = range(order + 1)
        half_window = (window_size - 1) / 2
        b = mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = linalg.pinv(b).A[deriv]

# pad the signal at the extremes with values taken from the signal itself

        firstvals = y[0] - abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + abs(y[-half_window-1:-1][::-1] - y[-1])
        y = concatenate((firstvals, y, lastvals))
        sg = convolve(m,y,mode='valid')

    return sg, status


# -----------------------------------------------------------
def running_frac_std(time,flux,wid,sig=None):

# calculate running fractional standard deviation across the array flux within a window of width wid

    hwid = wid / 2
    runstd = zeros(len(flux))
    for i in range(len(time)):
        valsinwid = flux[logical_and(time < time[i] + hwid, time > time[i] - hwid)]
        if sig is None:
            runstd[i] = std(valsinwid) / mean(valsinwid) 
        else:
            runstd[i] = std(sig_clip(valsinwid,sig)) / mean(valsinwid)

    return array(runstd)

