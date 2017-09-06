from .keparray import rebin2D
import math
import numpy as np
import sys
import os
import glob
from math import modf, cos, sin, radians, exp
from scipy import ndimage, interpolate
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate
from scipy.interpolate import RectBivariateSpline, interp2d
from . import kepio, kepmsg


def poly0(p, x):
    return p[0] + 0.0 * x

def poly1(p, x):
    return p[0] + p[1] * x

def poly2(p, x):
    return p[0] + p[1] * x + p[2] * x * x

def poly3(p, x):
    return p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3

def poly4(p, x):
    return p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4

def poly5(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x**5)

def poly6(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x ** 5 + p[6] * x ** 6)

def poly7(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x ** 5 + p[6] * x ** 6 + p[7] * x ** 7)

def poly8(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x ** 5 + p[6] * x ** 6 + p[7] * x ** 7 + p[8] * x ** 8)

def poly9(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x ** 5 + p[6] * x ** 6 + p[7] * x ** 7 + p[8] * x ** 8
            + p[9] * x**9)

def poly10(p, x):
    return (p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4
            + p[5] * x ** 5 + p[6] * x ** 6 + p[7] * x ** 7 + p[8] * x ** 8
            + p[9] * x ** 9 + p[10] * x ** 10)

def poly1con(p, x):
    return p[0] + x

def gauss(p, x):
    return p[0] * np.exp(-(x - p[1]) ** 2 / (2.0 * p[2] ** 2))

def gauss0(p, x):
    return p[0] * np.exp(-x**2 / (2.0 * p[1] ** 2))

def congauss(p, x):
    return p[0] + p[1] * np.exp(-(x - p[2]) ** 2 / (2.0 * p[3] ** 2))

def moffat0(p, x):
    return p[0] / (1.0 + (x / p[1]) ** 2) ** p[2]

def conmoffat(p, x):
    return p[0] + p[1] / (1.0 + ((x - p[2]) / p[3]) ** 2) ** p[4]

def sine(p, x):
    return p[0] * np.sin(2.0 * np.pi * x / p[1] - p[2])

def smooth(x, window_len=10, window='hanning'):
    """Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
     x: the input signal
     window_len: the dimension of the smoothing window
     window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
         flat window will produce a moving average smoothing.

    output:
     the smoothed signal

    example:
    t = linspace(-2, 2, 0.1)
    x = sin(t) + randn(len(t)) * 0.1
    y = smooth(x)

    see also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """
    window_len = int(window_len)

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")

    s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval("np.{0}(window_len)".format(window))

    y = np.convolve(w / w.sum(), s, mode='same')

    return y[window_len-1:-window_len+1]

def pei(law,wave,ebmv,rv,a_i,lambda_i,b_i,n_i):
    """redden a spectrum using Pei Y.C., 1992 ApJ, 395, 130
    Rv = 3.08 : Milky Way (1)
    Rv = 3.16 : LMC       (2)
    Rv = 2.93 : SMC       (3)
    """

    # extinction at B (a_b)
    a_b = ebmv * (1. + rv)
    # convert Angstroms to microns
    wave = wave / 1e4
    # build function
    xi = 0.
    for i in range(6):
        term  = math.pow((wave / lambda_i[law, i]), n_i[law, i])
        term += math.pow((lambda_i[law, i] / wave), n_i[law, i])
        term += b_i[law, i]
        term  = a_i[law, i] / term
        xi   += term
    # remove a_b normalization on the extinction curve
    a_lambda = a_b * xi
    if (wave < 0.08):
        a_lambda = 0.
    # linearize extinction factor
    return 10.**(-a_lambda / 2.512)

def pei_paramters():
    """
    Data from Pei Y.C., 1992 ApJ, 395, 130 (Table 4).
    Rv = 3.08 : Milky Way (1)
    Rv = 3.16 : LMC       (2)
    Rv = 2.93 : SMC       (3)
    """
    a_i = np.zeros([4, 6])
    lambda_i = np.zeros([4, 6])
    b_i = np.zeros([4, 6])
    n_i = np.zeros([4, 6])

    # Milky Way Extinction Law
    a_i[1, 0] = 165.        # BKG
    a_i[1, 1] = 14.         # FUV
    a_i[1, 2] = 0.045       # 2175 AA
    a_i[1, 3] = 0.002       # 9.7 um
    a_i[1, 4] = 0.002       # 18 um
    a_i[1, 5] = 0.012       # FIR

    lambda_i[1, 0] =  0.047 # BKG
    lambda_i[1, 1] =  0.08  # FUV
    lambda_i[1, 2] =  0.22  # 2175 AA
    lambda_i[1, 3] =  9.7   # 9.7 um
    lambda_i[1, 4] =  18.   # 18 um
    lambda_i[1, 5] =  25.   # FIR

    b_i[1, 0] = 90.         # BKG
    b_i[1, 1] =  4.         # FUV
    b_i[1, 2] = -1.95       # 2175 AA
    b_i[1, 3] = -1.95       # 9.7 um
    b_i[1, 4] = -1.8        # 18 um
    b_i[1, 5] =  0.         # FIR

    n_i[1, 0] = 2.          # BKG
    n_i[1, 1] = 6.5         # FUV
    n_i[1, 2] = 2.          # 2175 AA
    n_i[1, 3] = 2.          # 9.7 um
    n_i[1, 4] = 2.          # 18 um
    n_i[1, 5] = 2.          # FIR

    return a_i, lambda_i, b_i, n_i

def polyval(x, c, tensor=True):
    """1-d polynomial interpolation"""

    c = np.array(c, ndmin=1, copy=0)
    if c.dtype.char in '?bBhHiIlLqQpP':
        c = c + 0.0
    if isinstance(x, (tuple, list)):
        x = np.asarray(x)
    if isinstance(x, np.ndarray) and tensor:
        c = c.reshape(c.shape + (1,) * x.ndim)

    c0 = c[-1] + x * 0
    for i in range(2, len(c) + 1) :
        c0 = c[-i] + c0 * x

    return c0

def polyval2d(x,y,c):
    """2-d polynomial interpolation"""
    try:
        x, y = np.array((x, y), copy=0)
    except:
        raise ValueError('x, y are incompatible')
    c = polyval(x, c)
    c = polyval(y, c, tensor=False)

    return c

def PRFgauss2d(params, *args):
    """
    2-d Gaussian interpolation

    notation from Vanderburg and Johnson (2014)
    pos: x, y position to extrapolate Gaussian to
    cen: x, y center of Gaussian
    A: amplitude of Gaussian
    sigma: x, y width of Gassian
    B: amplitude of rotation term
    D: background
    """

    # parameters
    cen = [params[0], params[1]]
    A = params[2]
    sigma = [params[3], params[4]]
    B = params[5]
    D = params[6]

    # arguments

    posx = args[0]
    posy = args[1]
    flux = args[2]

    dx = posx - cen[0]
    dy = posy - cen[1]
    z = np.square(dx) / sigma[0] ** 2 + np.square(dy) / sigma[1] ** 2
    g = A * np.exp(-z - B * dx * dy) + D

    res = np.square(flux - g)

    return res


def PRF2DET(flux, OBJx, OBJy, DATx, DATy, wx, wy, a, splineInterpolation):
    """
    PRF interpolation function
    """

    # trigonometry
    cosa = np.cos(radians(a))
    sina = np.sin(radians(a))

    # where in the pixel is the source position?
    PRFfit = np.zeros((np.size(DATy), np.size(DATx)))
    for i in range(len(flux)):
        FRCx, INTx = modf(OBJx[i])
        FRCy, INTy = modf(OBJy[i])
        if FRCx > 0.5:
            FRCx -= 1.0
            INTx += 1.0
        if FRCy > 0.5:
            FRCy -= 1.0
            INTy += 1.0
        FRCx = -FRCx
        FRCy = -FRCy

        # constuct model PRF in detector coordinates
        for (j, y) in enumerate(DATy):
            for (k, x) in enumerate(DATx):
                xx = x - INTx + FRCx
                yy = y - INTy + FRCy
                dx = xx * cosa - yy * sina
                dy = xx * sina + yy * cosa
                PRFfit[j, k] = PRFfit[j, k] + splineInterpolation(dy * wy, dx * wx) * flux[i]

    return PRFfit

def PRF(params, *args):
    """
    PRF model
    """
    # arguments
    DATx = args[0]
    DATy = args[1]
    DATimg = args[2]
    DATerr = args[3]
    nsrc = args[4]
    splineInterpolation = args[5]
    col = args[6]
    row = args[7]

    # parameters
    f = np.empty((nsrc))
    x = np.empty((nsrc))
    y = np.empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc + i]
        y[i] = params[nsrc * 2 + i]

    # calculate PRF model binned to the detector pixel size
    PRFfit = PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)
    # calculate the sum squared difference between data and model
    PRFres = np.nansum(np.square(DATimg - PRFfit))
    # keep the fit centered
    if max(abs(col - x[0]), abs(row - y[0])) > 10.0:
        PRFres = 1.0e300
    return PRFres

def PRFwithBackground(params, *args):
    """
    PRF model with variable background
    """
    # arguments
    DATx = args[0]
    DATy = args[1]
    DATimg = args[2]
    DATerr = args[3]
    nsrc = args[4]
    bterms = args[5] + 1
    bx = args[6]
    by = args[7]
    splineInterpolation = args[8]
    col = args[9]
    row = args[10]

    # parameters
    f = np.empty((nsrc))
    x = np.empty((nsrc))
    y = np.empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc + i]
        y[i] = params[nsrc * 2 + i]
    b = np.array([params[nsrc * 3:nsrc * 3 + bterms],
               params[nsrc * 3+ bterms:nsrc * 3 + bterms * 2]])

    # calculate PRF model binned to the detector pixel size
    PRFfit = PRF2DET(f, x, y, DATx, DATy, 1.0, 1.0, 0.0, splineInterpolation)

    # add background
    if bterms == 1:
        PRFfit += params[nsrc * 3]
    else:
        PRFfit += polyval2d(bx, by, b)

    # calculate the sum squared difference between data and model
    PRFres = np.nansum(np.square(DATimg - PRFfit) / np.square(DATerr))

    # keep the fit centered
    if max(abs(col - x[0]), abs(row - y[0])) > 5.0:
        PRFres = 1.0e300

    return PRFres

def PRFwithFocusAndBackground(params, *args):
    """PRF model with variable focus and background"""
    DATx = args[0]
    DATy = args[1]
    DATimg = args[2]
    DATerr = args[3]
    nsrc = args[4]
    bterms = args[5] + 1
    bx = args[6]
    by = args[7]
    splineInterpolation = args[8]
    col = args[9]
    row = args[10]

    # parameters
    f = np.empty((nsrc))
    x = np.empty((nsrc))
    y = np.empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc + i]
        y[i] = params[nsrc * 2 + i]
    if bterms == 1:
        b = params[nsrc * 3]
    else:
        b = np.array([params[nsrc * 3:nsrc * 3 + bterms],
                   params[nsrc * 3 + bterms:nsrc * 3 + bterms * 2]])
    wx = params[-3]
    wy = params[-2]
    a = params[-1]

    try:
        PRFfit = PRF2DET(f, x, y, DATx, DATy, wx, wy, a, splineInterpolation)

        # add background
        if bterms == 1:
            PRFfit = PRFfit + b
        else:
            PRFfit = PRFfit + polyval2d(bx, by, b)

        # calculate the sum squared difference between data and model
        PRFres = np.nansum(np.square(DATimg - PRFfit) / np.square(DATerr))
    except:
        PRFres = 1.0e30

    # keep the fit centered
    if max(abs(col - x[0]),abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres

def PRFwithFocus(params, *args):
    """PRF model with variable focus"""

    # arguments
    DATx = args[0]
    DATy = args[1]
    DATimg = args[2]
    DATerr = args[3]
    nsrc = args[4]
    splineInterpolation = args[5]
    col = args[6]
    row = args[7]

    # parameters
    f = np.empty((nsrc))
    x = np.empty((nsrc))
    y = np.empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]
    wx = params[-3]
    wy = params[-2]
    a = params[-1]

    # iterate over sources
    try:
        PRFfit = PRF2DET(f, x, y, DATx, DATy, wx, wy, 0.0, splineInterpolation)

        # calculate the sum squared difference between data and model
        PRFres = np.nansum(np.square(DATimg - PRFfit) / np.square(DATerr))
    except:
        PRFres = 1.0e30

    # keep the fit centered
    if max(abs(col - x[0]), abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres

def kepler_prf_2d(params, *args):
    """the residual between pixel data and 2D Kepler PRF model"""
    data = args[0]
    prf = args[1]
    prfDelY = args[2]
    prfDelX = args[3]
    prfDimY = args[4]
    prfDimX = args[5]
    prfY0 = args[6]
    prfX0 = args[7]
    interpolation = args[8]
    verbose = args[9]
    f, y, x = params

    # interpolate PRF centroid to new pixel position
    model = shift(prf, [y, x], order=1, mode='constant')

    # extract the PRF model within the data limits
    model = model[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX]

    # rebin the PRF image to the same size and dimension of the data image
    model = rebin2D(model, [np.shape(data)[0], np.shape(data)[1]],
                    interpolation, True, False)
    model = model / prfDelY / prfDelX

    # calculate the sum squared difference between data and model
    residual = np.nansum(np.square(data - model * f))

    # write out parameters
    if verbose:
        txt = ("\rPearson\'s chisq = {0} for {1} dof"
               .format(np.nansum(np.square(data - model * f) / np.absolute(data)),
                       (np.shape(data)[0] * np.shape(data)[1] - len(params))))
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual

def kepler_multi_prf_2d(params, *args):
    """the residual between pixel data and 2D Kepler multiple PRF model"""
    # arguments
    data = args[0]
    prf = args[1]
    prfDelY = args[2]
    prfDelX = args[3]
    prfDimY = args[4]
    prfDimX = args[5]
    prfY0 = args[6]
    prfX0 = args[7]
    interpolation = args[8]
    verbose = args[9]

    # parameters
    nsrc = len(params) // 3
    f = np.empty((nsrc))
    y = np.empty((nsrc))
    x = np.empty((nsrc))
    model = np.zeros((prfDimY + 1, prfDimX + 1))
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc + i]
        x[i] = params[nsrc * 2 + i]

        # interpolate PRF centroid to new pixel position
        tmp = shift(prf, [y[i], x[i]], order=1, mode='constant')

        # extract the PRF model within the data limits
        model = model + tmp[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX] * f[i]

    # rebin the PRF image to the same size and dimension of the data image
    model = rebin2D(model, [np.shape(data)[0], np.shape(data)[1]], interpolation,
                    True, False)
    model = model / prfDelY / prfDelX

    # calculate the sum squared difference between data and model
    residual = np.nansum(np.square(data - model))

    # write out parameters
    if verbose:
        txt = ("\rPearson\'s chisq = {0} for {1} dof"
               .format(np.nansum(square(data - model) / np.absolute(data)),
                       (np.shape(data)[0] * np.shape(data)[1] - len(params))))
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual

def kepler_bkg_multi_prf_2d(params, *args):
    """the residual between pixel data and 2D Kepler multiple PRF model with background"""

    # arguments
    data = args[0]
    prf = args[1]
    prfDelY = args[2]
    prfDelX = args[3]
    prfDimY = args[4]
    prfDimX = args[5]
    prfY0 = args[6]
    prfX0 = args[7]
    interpolation = args[8]
    verbose = args[9]

    # parameters
    nsrc = (len(params) - 1) // 3
    f = np.empty((nsrc))
    y = np.empty((nsrc))
    x = np.empty((nsrc))
    b = params[nsrc * 3]
    model = np.zeros((prfDimY + 1,prfDimX + 1))
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc + i]
        x[i] = params[nsrc * 2 + i]

        # interpolate PRF centroid to new pixel position
        tmp = shift(prf, [y[i], x[i]], order=1, mode='constant')

        # extract the PRF model within the data limits
        model = model + tmp[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX] * f[i]

    # rebin the PRF image to the same size and dimension of the data image
    model = rebin2D(model, [np.shape(data)[0], np.shape(data)[1]], interpolation, True,
                    False)
    model = model / prfDelY / prfDelX
    model = model + b

    # calculate the sum squared difference between data and model
    residual = np.nansum(np.square(data - model))

    # write out parameters
    if verbose:
        txt = ("\rPearson\'s chisq = {0} for {1} dof"
               .format(int(np.nansum(np.square(data - model) / data)),
                       (np.shape(data)[0] * np.shape(data)[1] - len(params))))
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual

def kepler_focus_multi_prf_2d(params,*args):
    """the residual between pixel data and 2D Kepler multiple PRF model with background and focus"""

    # arguments
    data = args[0]
    prf = args[1]
    prfDelY = args[2]
    prfDelX = args[3]
    datDimX = args[4]
    datDimY = args[5]
    interpolation = args[6]
    verbose = args[7]

    # parameters
    nsrc = (len(params) - 2) // 3
    f = np.empty((nsrc))
    y = np.empty((nsrc))
    x = np.empty((nsrc))
    b = params[nsrc * 3]
    w = params[nsrc * 3 + 1]
    if w > 1.5:
        w = 1.5
    elif w < 1.0:
        w = 1.0
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc + i]
        x[i] = params[nsrc * 2 + i]

    # dimensions of data image if it had PRF-sized pixels
    prfDimY = datDimY / prfDelY / w
    prfDimX = datDimX / prfDelX / w
    print(w, prfDimY, prfDimX)

    # location of the data image centered on the PRF image (in PRF pixel units)
    prfY0 = (np.shape(prf)[0] - prfDimY) / 2
    prfX0 = (np.shape(prf)[1] - prfDimX) / 2

    # iterate over each source identified in the mask, build model
    DY = 0.0; DX = 0.0
    if int(prfDimY) % 2 == 0: DY = 1.0
    if int(prfDimX) % 2 == 0: DX = 1.0
    model = np.zeros((prfDimY+DY,prfDimX+DX))

    for i in range(nsrc):
        # interpolate PRF centroid to new pixel position
        tmp = shift(prf,[y[i]/w,x[i]/w],order=1,mode='constant')
        # extract the PRF model within the data limits
        model = model + tmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX] * f[i]
        # rebin the PRF image to the same size and dimension of the data image
    model = rebin2D(model,[np.shape(data)[0],np.shape(data)[1]],interpolation,True,False)
    model = model / prfDelY / prfDelX / w / w
    # add background to model
    model = model + b
    # calculate the sum squared difference between data and model
    residual = np.nansum(np.square(data - model))

    # write out parameters
    if verbose:
        txt = ("\rPearson\'s chisq = {0} for {1} dof"
               .format(int(np.nansum(np.square(data - model) / data)),
                       (np.shape(data)[0] * np.shape(data)[1] - len(params))))
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual

def kepler_prf_1d(params,*args):
    """Residual between pixel data and 2D Kepler PRF model integrated over x and y"""

    data = args[0]
    prf = args[1]
    prfDelY = args[2]
    prfDelX = args[3]
    prfDimY = args[4]
    prfDimX = args[5]
    prfY0 = args[6]
    prfX0 = args[7]
    interpolation = args[8]
    verbose = args[9]
    fy, fx, y, x = params

    # integrate along X and Y
    dataY = data.sum(axis=1)
    dataX = data.sum(axis=0)
    prfY = prf.sum(axis=1)
    prfX = prf.sum(axis=0)

    # interpolate PRF centroid to new pixel position
    modelY = shift(prfY, [y], order=1, mode='constant')
    modelX = shift(prfX, [x], order=1, mode='constant')

    # extract the PRF model within the data limits
    modelY = modelY[prfY0:prfY0 + prfDimY]
    modelX = modelX[prfX0:prfX0 + prfDimX]

    # rebin the PRF image to the same size and dimension of the data image
    modelY = rebin2D(modelY, [np.shape(data)[0]], interpolation, True, False)
    modelY = modelY / prfDelY
    modelX = rebin2D(modelX, [np.shape(data)[1]], interpolation, True, False)
    modelX = modelX / prfDelX

    # calculate the sum squared difference between data and model
    residualY = np.nansum(np.square(dataY - modelY * fy))
    residualX = np.nansum(np.square(dataX - modelX * fx))

    return residualY + residualX

def BKJD2BJD(bkjd):
    """convert BKJD to BJD"""

    return bkjd + 2454833.0

def BJD2BKJD(bjd):
    """convert BJD to BKJD"""

    return bjd - 2454833.0


def inv_normal_cummulative_function(p):
    """inverse normal cummulative function

    Lower tail quantile for standard normal distribution function.

    This function returns an approximation of the inverse cumulative
    standard normal distribution function.  I.e., given P, it returns
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

    # Coefficients in rational approximations.
    a = [-3.969683028665376e+01, 2.209460984245205e+02,
          -2.759285104469687e+02, 1.383577518672690e+02,
          -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02,
          -1.556989798598866e+02, 6.680131188771972e+01,
          -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
          -2.400758277161838e+00, -2.549732539343734e+00,
          4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00]

    # Define break-points.
    plow  = 0.02425
    phigh = 1.0 - plow

    # Rational approximation for lower region.
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return ((((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1))

    # Rational approximation for upper region.
    if phigh < p:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -((((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                 ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1))

    # Rational approximation for central region.
    q = p - 0.5
    r = q * q
    return ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r+ b[3]) * r + b[4]) * r + 1))

def read_and_interpolate_prf(prfdir, module, output, column, row, xdim, ydim,
                             verbose=False, logfile='kepprf.log'):
    """
    Read PRF file and prepare the data to be used for evaluating PRF.

    Parameters
    ----------
    prfdir : str
        The full or relative directory path to a folder containing the Kepler
        PSF calibration. Calibration files can be downloaded from the Kepler
        focal plane characteristics page at the MAST.
    module : str
        The 'MODULE' keyword from TPF file.
    output : str
        The 'OUTPUT' keyword from TPF file.
    column : int
        The '1CRV5P' keyword from TPF[1] file.
    row : int
        The '2CRV5P' keyword from TPF[1] file.
    xdim : int
        The first part of the 'TDIM5' keyword from TPF[1] file.
    ydim : int
        The second part of the 'TDIM5' keyword from TPF[1] file.
    verbose : boolean
        Print informative messages and warnings to the shell and logfile?
    logfile : string
        Name of the logfile containing error and warning messages.

    Returns
    -------
    splineInterpolation
        You can get PRF at given position: 
        kepfunc.PRF2DET([flux], [x], [y], DATx, DATy, 1.0, 1.0, 0.0, splineInterpolation)
    DATx : numpy.array
        X-axis coordiantes of pixels for given TPF
    DATy : numpy.array
        Y-axis coordinates of pixels for given TPF
    prf : numpy.array
        PRF interpolated to given position on the camera
    PRFx : numpy.array
        X-axis coordinates of prf values
    PRFy : numpy.array
        Y-axis coordinates of prf values
    PRFx0 : int
    PRFy0 : int
    cdelt1p : numpy.array 
        CDELT1P values from 5 HDUs of PRF file.
    cdelt2p : numpy.array
        CDELT2P values from 5 HDUs of PRF file.
    prfDimX : int
        size of PRFx
    prfDimY : int
        size of PRFy

    """
    n_hdu = 5
    minimum_prf_weight = 1.e-6

    # determine suitable PRF calibration file
    if int(module) < 10:
        prefix = 'kplr0'
    else:
        prefix = 'kplr'
    prfglob = os.path.join(prfdir, prefix + module + '.' + output + '*_prf.fits')
    try:
        prffile = glob.glob(prfglob)[0]
    except:
        errmsg = "ERROR -- KEPPRF: No PRF file found in {0}".format(prfdir)
        kepmsg.err(logfile, errmsg, verbose)

    # read PRF images
    prfn = [0] * n_hdu
    crval1p = np.zeros(n_hdu, dtype='float32')
    crval2p = np.zeros(n_hdu, dtype='float32')
    cdelt1p = np.zeros(n_hdu, dtype='float32')
    cdelt2p = np.zeros(n_hdu, dtype='float32')
    for i in range(n_hdu):
        (prfn[i], _, _, crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i]) = \
            kepio.readPRFimage(prffile, i+1, logfile, verbose)
    prfn = np.array(prfn)
    PRFx = np.arange(0.5, np.shape(prfn[0])[1] + 0.5)
    PRFy = np.arange(0.5, np.shape(prfn[0])[0] + 0.5)
    PRFx = (PRFx - np.size(PRFx) / 2) * cdelt1p[0]
    PRFy = (PRFy - np.size(PRFy) / 2) * cdelt2p[0]

    # interpolate the calibrated PRF shape to the target position
    prf = np.zeros(np.shape(prfn[0]), dtype='float32')
    prfWeight = np.zeros(n_hdu, dtype='float32')
    ref_column = column + (xdim - 1.) / 2.
    ref_row = row + (ydim - 1.) / 2.
    for i in range(n_hdu):
        prfWeight[i] = math.sqrt(
            (ref_column - crval1p[i])**2 + (ref_row - crval2p[i])**2)
        if prfWeight[i] < minimum_prf_weight:
            prfWeight[i] = minimum_prf_weight
        prf += prfn[i] / prfWeight[i]
    prf /= (np.nansum(prf) * cdelt1p[0] * cdelt2p[0])

    # location of the data image centered on the PRF image (in PRF pixel units)
    prfDimY = int(ydim / cdelt1p[0])
    prfDimX = int(xdim / cdelt2p[0])
    PRFy0 = int(np.round((np.shape(prf)[0] - prfDimY) / 2))
    PRFx0 = int(np.round((np.shape(prf)[1] - prfDimX) / 2))
    DATx = np.arange(column, column + xdim)
    DATy = np.arange(row, row + ydim)

    # interpolation function over the PRF
    splineInterpolation = RectBivariateSpline(PRFx, PRFy, prf)

    return (splineInterpolation, DATx, DATy, prf, PRFx, PRFy, PRFx0, PRFy0,
            cdelt1p, cdelt2p, prfDimX, prfDimY)

