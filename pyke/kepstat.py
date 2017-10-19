import math
import random
import numpy as np
from scipy import linalg
from . import kepmsg


def mean_err(array):
    return math.sqrt(np.nansum(array * array)) / len(array)


def rms(array1, array2, logfile, verbose):
    """root mean square of two arrays"""
    if len(array1) != len(array2):
        message  = ("ERROR -- KEPSTAT.RMS: Arrays have unequal sizes - "
                    "array1 = {0}, array2 = {1}".format(len(array1),
                                                        len(array2)))
        kepmsg.err(logfile, message, verbose)

    return math.sqrt(np.nanmean((array1 - array2) ** 2))


def randarray(signal, err):
    """adjust data relative to random number constrained by error bars"""
    random.seed()
    out = np.zeros(len(signal), dtype='float32')
    for i in range(len(signal)):
        out[i] = signal[i] + err[i] * inv_normal_cummulative_function(random.random())
    return out


def removeinfinlc(x, cols):
    """remove infinities from light curve data"""
    for j in range(len(cols)):
        work = []
        datatype = cols[j].dtype
        for i in range(len(x)):
            if np.isfinite(x[i]):
                work.append(cols[j][i])
        cols[j] = np.array(work, dtype=datatype)
    return cols


def filterOnRange(intime,tstart,tstop):
    """filter on data within time ranges"""
    outlist = []
    for i in range(len(intime)):
        for j in range(len(tstart)):
            if intime[i] > tstart[j] and intime[i] < tstop[j]:
                if len(outlist) == 0:
                    outlist.append(i)
                elif i > outlist[-1]:
                    outlist.append(i)
    return outlist


def inv_normal_cummulative_function(p):
    """
    Lower tail quantile for standard normal distribution function.

    This function returns an apprximation of the inverse cumulative
    standard normal distribution function.  i.e., given P, it returns
    an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    random variable from the standard normal distribution.

    The algorithm uses a minimax approximation by rational functions
    and the result has a relative error whose absolute value is less
    than 1.15e-9.

    Author:      Peter J. Acklam
    Time-stamp:  2000-07-19 18:26:14
    E-mail:      pjacklam@online.no
    WWW URL:     http://home.online.no/~pjacklam
    """

    if p == 0.0:
        p = 0.0000000001
    if p == 1.0:
        p = 0.9999999999

# coefficients in rational approximations
    a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
         1.383577518672690e2, -3.066479806614716e1, 2.506628277459239]
    b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
         6.680131188771972e1, -1.328068155288572e1]
    c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838,
         -2.549732539343734, 4.374664141464968, 2.938163982698783]
    d = [7.784695709041462e-3, 3.224671290700398e-1,
         2.445134137142996, 3.754408661907416]
# define break-points
    plow = 0.02425
    phigh = 1.0 - plow
# rational approximation for lower region
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return ((((((c[0] * q + c[1]) * q + c[2]) * q+c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0))
# rational approximation for upper region
    if phigh < p:
        q = math.sqrt(-2.0 * math.log(1 - p))
        return -((((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                 ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0))
# rational approximation for central region
    q = p - 0.5
    r = q * q
    return ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0))

def bitInBitmap(bitmap, bit):
    """bit map decoding"""
    flag = False
    for i in range(10, -1,- 1):
        if bitmap - 2**i >= 0:
            bitmap = bitmap - 2**i
            if 2**i == bit:
                flag = True
        else:
            continue
    return flag

def savitzky_golay(y, window_size, order, deriv=0):
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

    sg = 0.0

    try:
        window_size = abs(int(window_size))
        order = abs(int(order))
    except:
        message = ("ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size and order "
                   "must be of type int")
        kepmsg.err(None, message, True)
    if window_size % 2 != 1 or window_size < 1:
        message = ("ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size size must be "
                   "a positive odd number")
        kepmsg.err(None, message, True)
    if window_size < order + 2:
        message = ("ERROR -- KEPSTAT.SAVITZKY_GOLAY: window_size is too small "
                   "for the polynomials order")
        kepmsg.err(None, message, True)

    # precompute fit coefficients
    order_range = range(order + 1)
    half_window = (window_size - 1) / 2
    b = mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = linalg.pinv(b).A[deriv]

    # pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    sg = np.convolve(m, y, mode='valid')

    return sg


def running_frac_std(time, flux, wid):
    """calculate running fractional standard deviation across the array flux
       within a window of width wid
    """

    hwid = wid / 2
    runstd = np.zeros(len(flux))
    for i in range(len(time)):
        valsinwid = flux[np.logical_and(time < time[i] + hwid, time > time[i] - hwid)]
        runstd[i] = np.std(valsinwid) / np.mean(valsinwid)

    return np.array(runstd)
