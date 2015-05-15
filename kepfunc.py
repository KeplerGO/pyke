import numpy, scipy, math, sys
import keparray
from math import modf, cos, sin, radians, exp
from scipy import ndimage, interpolate
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate
from scipy.interpolate import RectBivariateSpline, interp2d
from keparray import rebin2D
from numpy import square, nansum, shape, array, empty, zeros, absolute, size
from sys import stdout, exit
 
# -----------------------------------------------
# define functions

def poly0(): 
    return lambda p, x: p[0] + 0.0 * x

def poly1(): 
    return lambda p, x: p[0] + p[1] * x

def poly2(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x * x

def poly3(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3

def poly4(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4

def poly5(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5

def poly6(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5 + p[6] * x**6

def poly7(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5 + p[6] * x**6 + p[7] * x**7

def poly8(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5 + p[6] * x**6 + p[7] * x**7 + p[8] * x**8

def poly9(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5 + p[6] * x**6 + p[7] * x**7 + p[8] * x**8 + p[9] * x**9

def poly10(): 
    return lambda p, x: p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5 + p[6] * x**6 + p[7] * x**7 + p[8] * x**8 + p[9] * x**9 + p[10] * x**10

def poly1con(): 
    return lambda p, x: p[0] + x

def gauss(): 
    return lambda p, x: p[0] * scipy.exp(-(x - p[1])**2 / (2.0 * p[2]**2))

def gauss0(): 
    return lambda p, x: p[0] * scipy.exp(-x**2 / (2.0 * p[1]**2))

def congauss(): 
    return lambda p, x: p[0] + p[1] * scipy.exp(-(x - p[2])**2 / (2.0 * p[3]**2))

def moffat0(): 
    return lambda p, x: p[0] / (1.0 + (x / p[1])**2)**p[2]

def conmoffat(): 
    return lambda p, x: p[0] + p[1] / (1.0 + ((x - p[2]) / p[3])**2)**p[4]

def sine():
    return lambda p, x: p[0] * scipy.sin(2.0 * 3.14129 * x / p[1] - p[2])

def powerlaw():
    return lambda p, x: p[0] + p[1] * x

# -----------------------------------------------
# smooth the data using a window with requested size

def smooth(x,window_len=10,window='hanning'):

     """smooth the data using a window with requested size.
     
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
 
     t=linspace(-2,2,0.1)
     x=sin(t)+randn(len(t))*0.1
     y=smooth(x)
     
     see also: 
     
     numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
     scipy.signal.lfilter
  
     TODO: the window parameter could be the window itself if an array instead of a string   
     """
 
     if x.ndim != 1:
         raise ValueError, "smooth only accepts 1 dimension arrays."
 
     if x.size < window_len:
         raise ValueError, "Input vector needs to be bigger than window size."
 
 
     if window_len<3:
         return x
 
 
     if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
 
 
     s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
     if window == 'flat': #moving average
         w=numpy.ones(window_len,'d')
     else:
         w=eval('numpy.'+window+'(window_len)')
 
     y=numpy.convolve(w/w.sum(),s,mode='same')
     return y[window_len-1:-window_len+1]

#-------------------------------
def pei(law,wave,ebmv,rv,a_i,lambda_i,b_i,n_i):

# redden a spectrum using Pei Y.C., 1992 ApJ, 395, 130
# Rv = 3.08 : Milky Way (1)
# Rv = 3.16 : LMC       (2)
# Rv = 2.93 : SMC       (3)

# extinction at B (a_b)

    a_b = ebmv * (1. + rv)

# convert Angstroms to microns

    wave = wave / 1e4

# build function

    xi = 0.
    for i in range(6):
	term  = math.pow((wave / lambda_i[law,i]),n_i[law,i])
	term += math.pow((lambda_i[law,i] / wave),n_i[law,i])
	term += b_i[law,i]
	term  = a_i[law,i] / term
	xi   += term

# remove a_b normalization on the extinction curve

    a_lambda=a_b*xi
    if (wave < 0.08): a_lambda = 0.

# linearize extinction factor

    return 10.**(-a_lambda/2.512)

#-------------------------------
def pei_paramters():

# Data from Pei Y.C., 1992 ApJ, 395, 130 (Table 4).
# Rv = 3.08 : Milky Way (1)
# Rv = 3.16 : LMC       (2)
# Rv = 2.93 : SMC       (3)

    a_i = numpy.zeros([4,6])
    lambda_i = numpy.zeros([4,6])
    b_i = numpy.zeros([4,6])
    n_i = numpy.zeros([4,6])

# Milky Way Extinction Law

    a_i[1,0] = 165.        # BKG
    a_i[1,1] = 14.         # FUV
    a_i[1,2] = 0.045       # 2175 AA
    a_i[1,3] = 0.002       # 9.7 um
    a_i[1,4] = 0.002       # 18 um
    a_i[1,5] = 0.012       # FIR

    lambda_i[1,0] =  0.047 # BKG
    lambda_i[1,1] =  0.08  # FUV
    lambda_i[1,2] =  0.22  # 2175 AA
    lambda_i[1,3] =  9.7   # 9.7 um
    lambda_i[1,4] =  18.   # 18 um
    lambda_i[1,5] =  25.   # FIR

    b_i[1,0] = 90.         # BKG
    b_i[1,1] =  4.         # FUV
    b_i[1,2] = -1.95       # 2175 AA
    b_i[1,3] = -1.95       # 9.7 um
    b_i[1,4] = -1.8        # 18 um
    b_i[1,5] =  0.         # FIR

    n_i[1,0] = 2.          # BKG
    n_i[1,1] = 6.5         # FUV
    n_i[1,2] = 2.          # 2175 AA
    n_i[1,3] = 2.          # 9.7 um
    n_i[1,4] = 2.          # 18 um
    n_i[1,5] = 2.          # FIR

    return a_i, lambda_i, b_i, n_i


#------------------------------
# 1-d polynomial interpolation 

def polyval(x,c,tensor=True):

    c = array(c,ndmin=1,copy=0)
    if c.dtype.char in '?bBhHiIlLqQpP':
        c = c + 0.0
    if isinstance(x, (tuple, list)):
        x = asarray(x)
    if isinstance(x, numpy.ndarray) and tensor:
        c = c.reshape(c.shape + (1,)*x.ndim)

    c0 = c[-1] + x*0
    for i in range(2, len(c) + 1) :
        c0 = c[-i] + c0*x

    return c0


#------------------------------
# 2-d polynomial interpolation 

def polyval2d(x,y,c):

    try:
        x,y = array((x,y),copy=0)
    except:
        raise ValueError('x, y are incompatible')
    c = polyval(x,c)
    c = polyval(y,c,tensor=False)

    return c


#------------------------------
# 2-d Gaussian interpolation 

def PRFgauss2d(params,*args):

# notation from Vanderburg and Johnson (2014)
# pos: x, y position to extrapolate Gaussian to
# cen: x, y center of Gaussian
# A: amplitude of Gaussian
# sigma: x, y width of Gassian
# B: amplitude of rotation term
# D: background

# parameters

    cen = [params[0],params[1]]
    A = params[2]
    sigma = [params[3],params[4]]
    B = params[5]
    D = params[6]

# arguments

    posx = args[0]
    posy = args[1]
    flux = args[2]

    dx = posx - cen[0] 
    dy = posy - cen[1] 
    z = square(dx) / sigma[0]**2 + square(dy) / sigma[1]**2
    g = A * scipy.exp(-z - B * dx * dy) + D

    res = square(flux - g)

    return res

#------------------------------
# PRF interpolation function

def PRF2DET2(flux,OBJx,OBJy,DATx,DATy,splineInterpolation):

# where in the pixel is the source position?

    PRFfit = zeros((size(DATy),size(DATx)))
    for i in range(len(flux)):
        FRCx,INTx = modf(OBJx[i])
        FRCy,INTy = modf(OBJy[i])
        if FRCx > 0.5: 
            FRCx -= 1.0
            INTx += 1.0
        if FRCy > 0.5: 
            FRCy -= 1.0
            INTy += 1.0
        FRCx = -FRCx
        FRCy = -FRCy

# constuct model PRF in detector coordinates

        for (j,y) in enumerate(DATy):
            for (k,x) in enumerate(DATx):
                dy = y - INTy + FRCy
                dx = x - INTx + FRCx
                PRFfit[j,k] = PRFfit[j,k] + splineInterpolation(dy,dx) * flux[i]

    return PRFfit


#------------------------------
# PRF interpolation function

def PRF2DET(flux,OBJx,OBJy,DATx,DATy,wx,wy,a,splineInterpolation):

# trigonometry

    cosa = cos(radians(a))
    sina = sin(radians(a))

# where in the pixel is the source position?

    PRFfit = zeros((size(DATy),size(DATx)))
    for i in range(len(flux)):
        FRCx,INTx = modf(OBJx[i])
        FRCy,INTy = modf(OBJy[i])
        if FRCx > 0.5: 
            FRCx -= 1.0
            INTx += 1.0
        if FRCy > 0.5: 
            FRCy -= 1.0
            INTy += 1.0
        FRCx = -FRCx
        FRCy = -FRCy

# constuct model PRF in detector coordinates

        for (j,y) in enumerate(DATy):
            for (k,x) in enumerate(DATx):
                xx = x - INTx + FRCx
                yy = y - INTy + FRCy
                dx = xx * cosa - yy * sina
                dy = xx * sina + yy * cosa
                PRFfit[j,k] += PRFfit[j,k] + splineInterpolation(dy*wy,dx*wx) * flux[i]

    return PRFfit


#------------------------------
# PRF model 

def PRF(params,*args):

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
    
    f = empty((nsrc))
    x = empty((nsrc))
    y = empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]

# calculate PRF model binned to the detector pixel size

    PRFfit = PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)

# calculate the sum squared difference between data and model

#    PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))
    PRFres = nansum(square(DATimg - PRFfit))

# keep the fit centered

    if max(abs(col - x[0]),abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres


#------------------------------
# PRF model with variable background

def PRFwithBackground(params,*args):

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
    
    f = empty((nsrc))
    x = empty((nsrc))
    y = empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]
    b = array([params[nsrc*3:nsrc*3+bterms],params[nsrc*3+bterms:nsrc*3+bterms*2]]) 

# calculate PRF model binned to the detector pixel size

    PRFfit = PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)

# add background 

    if bterms == 1:
        PRFfit += params[nsrc*3]
    else:
        PRFfit += polyval2d(bx, by, b)

# calculate the sum squared difference between data and model

    PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))

# keep the fit centered

    if max(abs(col - x[0]),abs(row - y[0])) > 5.0:
        PRFres = 1.0e300

    return PRFres


#------------------------------
# PRF model with variable focus and background

def PRFwithFocusAndBackground(params,*args):

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
    
    f = empty((nsrc))
    x = empty((nsrc))
    y = empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]
    if bterms == 1:
        b = params[nsrc*3]
    else:
        b = array([params[nsrc*3:nsrc*3+bterms],params[nsrc*3+bterms:nsrc*3+bterms*2]]) 
    wx = params[-3]
    wy = params[-2]
    a = params[-1]

    try:
        PRFfit = PRF2DET(f,x,y,DATx,DATy,wx,wy,a,splineInterpolation)

# add background 

        if bterms == 1:
            PRFfit = PRFfit + b
        else:
            PRFfit = PRFfit + polyval2d(bx, by, b)

# calculate the sum squared difference between data and model

#        PRFres = nansum(square(DATimg - PRFfit))
#        PRFres = nansum(numpy.abs(square(DATimg - PRFfit) / PRFfit))
        PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))
    except:
        PRFres = 1.0e30
#    if wx > 1.15 or wx < 0.85:
#        PRFres = 1.0e30
#    if wy > 1.15 or wy < 0.85:
#        PRFres = 1.0e30

# keep the fit centered

    if max(abs(col - x[0]),abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres


#------------------------------
# PRF model with variable focus

def PRFwithFocus(params,*args):

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
    
    f = empty((nsrc))
    x = empty((nsrc))
    y = empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]
    wx = params[-3]
    wy = params[-2]
    a = params[-1]

# iterate over sources

    try:
        PRFfit = PRF2DET(f,x,y,DATx,DATy,wx,wy,0.0,splineInterpolation)

# calculate the sum squared difference between data and model

        PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))
    except:
        PRFres = 1.0e30
#    if wx > 1.15 or wx < 0.85:
#        PRFres = 1.0e30
#    if wy > 1.15 or wy < 0.85:
#        PRFres = 1.0e30

# keep the fit centered

    if max(abs(col - x[0]),abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres


#-----------------------------------------------------
# the residual between pixel data and 2D Kepler PRF model

def kepler_prf_2d(params,*args):

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
    f,y,x = params

# interpolate PRF centroid to new pixel position

    model = shift(prf,[y,x],order=1,mode='constant')

# extract the PRF model within the data limits

    model = model[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX]

# rebin the PRF image to the same size and dimension of the data image

    model = rebin2D(model,[numpy.shape(data)[0],numpy.shape(data)[1]],interpolation,True,False)
    model = model / prfDelY / prfDelX

# calculate the sum squared difference between data and model

    residual = nansum(square(data - model * f))

# write out parameters

    if verbose:
#        txt = '\rFlux = %.2f e-/s ' % f
#        txt += 'X = %.4f pix ' % (x * prfDelX)
#        txt += 'Y = %.4f pix ' % (y * prfDelY)
#        txt += ' ' * 5
#        sys.stdout.write(txt)
#        sys.stdout.flush()

        txt = '\rPearson\'s chisq = %d for %d dof' % \
            (int(nansum(square(data - model * f) / absolute(data))), (shape(data)[0] * shape(data)[1] - len(params)))
        txt += ' ' * 5 
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual


#-----------------------------------------------------------------
# the residual between pixel data and 2D Kepler multiple PRF model

def kepler_multi_prf_2d(params,*args):

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
    
    nsrc = len(params) / 3
    f = empty((nsrc))
    y = empty((nsrc))
    x = empty((nsrc))
    model = zeros((prfDimY+1,prfDimX+1))
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc+i]
        x[i] = params[nsrc*2+i]

# interpolate PRF centroid to new pixel position

        tmp = shift(prf,[y[i],x[i]],order=1,mode='constant')

# extract the PRF model within the data limits

        model = model + tmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX] * f[i]

# rebin the PRF image to the same size and dimension of the data image

    model = rebin2D(model,[shape(data)[0],shape(data)[1]],interpolation,True,False)
    model = model / prfDelY / prfDelX

# calculate the sum squared difference between data and model

    residual = nansum(square(data - model))

# write out parameters

    if verbose:
#        if nsrc == 1:
#            txt = '\rFlux = %.2f e-/s ' % f[1]
#            txt += 'X = %7.4f pix ' % (x[1] * prfDelX)
#            txt += 'Y = %7.4f pix ' % (y[1] * prfDelY)
#            txt += 'Pearson\'s chisq = %d for %d dof' % \
#                (int(nansum(square(data - model) / absolute(data))), (shape(data)[0] * shape(data)[1] - len(params)))
#        else:
        txt = '\rPearson\'s chisq = %d for %d dof' % \
            (int(nansum(square(data - model) / absolute(data))), (shape(data)[0] * shape(data)[1] - len(params)))
        txt += ' ' * 5 
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual


#---------------------------------------------------------------------------------
# the residual between pixel data and 2D Kepler multiple PRF model with background

def kepler_bkg_multi_prf_2d(params,*args):

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
    
    nsrc = (len(params) - 1) / 3
    f = empty((nsrc))
    y = empty((nsrc))
    x = empty((nsrc))
    b = params[nsrc*3]
    model = zeros((prfDimY+1,prfDimX+1))
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc+i]
        x[i] = params[nsrc*2+i]    

# interpolate PRF centroid to new pixel position

        tmp = shift(prf,[y[i],x[i]],order=1,mode='constant')

# extract the PRF model within the data limits

        model = model + tmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX] * f[i]

# rebin the PRF image to the same size and dimension of the data image

    model = rebin2D(model,[shape(data)[0],shape(data)[1]],interpolation,True,False)
    model = model / prfDelY / prfDelX
    model = model + b

# calculate the sum squared difference between data and model

    residual = nansum(square(data - model))

# write out parameters

    if verbose:
        txt = '\rPearson\'s chisq = %d for %d dof' % \
            (int(nansum(square(data - model) / data)), (shape(data)[0] * shape(data)[1] - len(params)))
        txt += ' ' * 5 
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual


#-------------------------------------------------------------------------------------------
# the residual between pixel data and 2D Kepler multiple PRF model with background and focus

def kepler_focus_multi_prf_2d(params,*args):

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
    
    nsrc = (len(params) - 2) / 3
    f = empty((nsrc))
    y = empty((nsrc))
    x = empty((nsrc))
    b = params[nsrc*3]
    w = params[nsrc*3+1]
    if w > 1.5:
        w = 1.5
    elif w < 1.0:
        w = 1.0
    for i in range(nsrc):
        f[i] = params[i]
        y[i] = params[nsrc+i]
        x[i] = params[nsrc*2+i]    

# dimensions of data image if it had PRF-sized pixels

    prfDimY = datDimY / prfDelY / w
    prfDimX = datDimX / prfDelX / w
    print w, prfDimY, prfDimX

# location of the data image centered on the PRF image (in PRF pixel units)

    prfY0 = (shape(prf)[0] - prfDimY) / 2
    prfX0 = (shape(prf)[1] - prfDimX) / 2

# iterate over each source identified in the mask, build model

    DY = 0.0; DX = 0.0
    if int(prfDimY) % 2 == 0: DY = 1.0
    if int(prfDimX) % 2 == 0: DX = 1.0
    model = zeros((prfDimY+DY,prfDimX+DX))
    for i in range(nsrc):

# interpolate PRF centroid to new pixel position

        tmp = shift(prf,[y[i]/w,x[i]/w],order=1,mode='constant')

# extract the PRF model within the data limits

        model = model + tmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX] * f[i]

# rebin the PRF image to the same size and dimension of the data image

    model = rebin2D(model,[shape(data)[0],shape(data)[1]],interpolation,True,False)
    model = model / prfDelY / prfDelX / w / w

# add background to model

    model = model + b

# calculate the sum squared difference between data and model

    residual = nansum(square(data - model))

# write out parameters

    if verbose:
        txt = '\rPearson\'s chisq = %d for %d dof' % \
            (int(nansum(square(data - model) / data)), (shape(data)[0] * shape(data)[1] - len(params)))
        txt += ' ' * 5 
        sys.stdout.write(txt)
        sys.stdout.flush()

    return residual


#--------------------------------------------------------------------------------
# the residual between pixel data and 2D Kepler PRF model integrated over x and y

def kepler_prf_1d(params,*args):

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
    fy,fx,y,x = params

# integrate along X and Y

    dataY = data.sum(axis=1)
    dataX = data.sum(axis=0)
    prfY = prf.sum(axis=1) 
    prfX = prf.sum(axis=0)

# interpolate PRF centroid to new pixel position

    modelY = shift(prfY,[y],order=1,mode='constant')
    modelX = shift(prfX,[x],order=1,mode='constant')

# extract the PRF model within the data limits

    modelY = modelY[prfY0:prfY0+prfDimY]
    modelX = modelX[prfX0:prfX0+prfDimX]

# rebin the PRF image to the same size and dimension of the data image

    modelY = rebin2D(modelY,[numpy.shape(data)[0]],interpolation,True,False)
    modelY = modelY / prfDelY
    modelX = rebin2D(modelX,[numpy.shape(data)[1]],interpolation,True,False)
    modelX = modelX / prfDelX

# calculate the sum squared difference between data and model

    residualY = nansum(square(dataY - modelY * fy))
    residualX = nansum(square(dataX - modelX * fx))
    return residualY + residualX

#--------------------------------------------------------------------------------
# convert BKJD to BJD

def BKJD2BJD(bkjd):

    return bkjd + 2454833.0

#--------------------------------------------------------------------------------
# convert BJD to BKJD

def BJD2BKJD(bjd):

    return bjd - 2454833.0

#------------------------------------------------------------------
# inverse normal cummulative function

def inv_normal_cummulative_function(p):

# Lower tail quantile for standard normal distribution function.
#
# This function returns an approximation of the inverse cumulative
# standard normal distribution function.  I.e., given P, it returns
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
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)

# Rational approximation for upper region.

    if phigh < p:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)

# Rational approximation for central region.

    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
        (((((b[0] * r + b[1]) * r + b[2]) * r+ b[3]) * r + b[4]) * r + 1)

