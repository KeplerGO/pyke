#!/usr/bin/env python

import kepmsg, kepstat, kepfunc, keparray
import pylab, scipy, numpy, math, random, sys
from pylab import polyval
from scipy import optimize, ndimage
from scipy.optimize import fmin_powell, fmin_tnc, fmin, leastsq
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift
from numpy import empty, zeros, shape, nansum, linspace
from keparray import rebin2D
from math import sqrt

# -----------------------------------------------------------
# linear least square polynomial fit using scipy

def leastsquare(functype,pinit,xdata,ydata,yerr,logfile,verbose):

    status = 0
    coeffs = []

# functional form

    if (functype == 'poly0'): fitfunc = kepfunc.poly0()
    if (functype == 'poly1'): fitfunc = kepfunc.poly1()
    if (functype == 'poly2'): fitfunc = kepfunc.poly2()
    if (functype == 'poly3'): fitfunc = kepfunc.poly3()
    if (functype == 'poly4'): fitfunc = kepfunc.poly4()
    if (functype == 'poly5'): fitfunc = kepfunc.poly5()
    if (functype == 'poly6'): fitfunc = kepfunc.poly6()
    if (functype == 'poly7'): fitfunc = kepfunc.poly7()
    if (functype == 'poly8'): fitfunc = kepfunc.poly8()
    if (functype == 'poly1con'): fitfunc = kepfunc.poly1con()
    if (functype == 'gauss'): fitfunc = kepfunc.gauss()
    if (functype == 'gauss0'): fitfunc = kepfunc.gauss0()
    if (functype == 'congauss'): fitfunc = kepfunc.congauss()
    if (functype == 'sine'): fitfunc = kepfunc.sine()
    if (functype == 'moffat0'): fitfunc = kepfunc.moffat0()
    if (functype == 'conmoffat'): fitfunc = kepfunc.conmoffat()

# define error coefficent calculation

    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

# if no data errors, substitude rms of fit

    if (yerr == None):
	yerr = []
	rerr = []
	for i in range(len(ydata)):
	    rerr.append(1.e10)
	try:
            out = optimize.leastsq(errfunc,pinit,args=(xdata,ydata,rerr),full_output=1)
	except:
	    message = 'ERROR -- KEPFIT.LEASTSQUARE: failed to fit data'
	    status = kepmsg.err(logfile,message,verbose)
            if functype == 'poly0':
                out = [numpy.mean(ydata),sqrt(numpy.mean(ydata))]
	if (functype == 'poly0' or functype == 'sineCompareBinPSF'):
	    coeffs.append(out[0])
	else:
	    coeffs = out[0]
	if (len(coeffs) > 1):
	    fit = fitfunc(coeffs,xdata)
	else:
	    fit = numpy.zeros(len(xdata))
	    for i in range(len(fit)):
		fit[i] = coeffs[0]
	sigma, status = kepstat.rms(ydata,fit,logfile,verbose)
	for i in range(len(ydata)):
	    yerr.append(sigma)

# fit data 

    try:
        out = optimize.leastsq(errfunc, pinit, args=(xdata, ydata, yerr), full_output=1)
    except:
	message = 'ERROR -- KEPFIT.LEASTSQUARE: failed to fit data'
	status = kepmsg.err(logfile,message,verbose)
        if functype == 'poly0':
            out = [numpy.mean(ydata),sqrt(numpy.mean(ydata))]

# define coefficients

    coeffs = []
    covar = []
    if (functype == 'poly0' or functype == 'poly1con' or 
	functype == 'sineCompareBinPSF'):
	coeffs.append(out[0])
	covar.append(out[1])
    else:
	coeffs = out[0]
	covar = out[1]

# calculate 1-sigma error on coefficients

    errors = []
    if (covar == None): 
	message = 'WARNING -- KEPFIT.leastsquare: NULL covariance matrix'
#	kepmsg.log(logfile,message,verbose)
    for i in range(len(coeffs)):
	if (covar != None and len(coeffs) > 1):
	    errors.append(sqrt(abs(covar[i][i])))
	else:
	    errors.append(coeffs[i])

# generate fit points for rms calculation

    if (len(coeffs) > 1):
	fit = fitfunc(coeffs,xdata)
    else:
	fit = numpy.zeros(len(xdata))
	for i in range(len(fit)):
	    fit[i] = coeffs[0]
    sigma, status = kepstat.rms(ydata,fit,logfile,verbose)

# generate fit points for plotting

    dx = xdata[len(xdata)-1] - xdata[0]
    plotx = linspace(xdata.min(),xdata.max(),10000)
    ploty = fitfunc(coeffs,plotx)
    if (len(coeffs) == 1):
	ploty = []
	for i in range(len(plotx)):
	    ploty.append(coeffs[0])
	ploty = numpy.array(ploty)

# reduced chi^2 calculation

    chi2 = 0
    dof = len(ydata) - len(coeffs)
    for i in range(len(ydata)):
	chi2 += (ydata[i] - fit[i])**2 / yerr[i]
    chi2 /= dof

    return coeffs, errors, covar, sigma, chi2, dof, fit, plotx, ploty, status

# -----------------------------------------------------------
# linear least square fit with sigma-clipping

def lsqclip(functype,pinit,x,y,yerr,rej_lo,rej_hi,niter,logfile,verbose):

# functype = unctional form
# pinit = initial guess for parameters
# x = list of x data
# y = list of y data
# yerr = list of 1-sigma y errors
# order = polynomial order
# rej_lo = lower rejection threshold (units=sigma)
# rej_hi = upper rejection threshold (units=sugma)
# niter = number of sigma-clipping iterations

    npts = []
    iiter = 0
    iterstatus = 1
    status = 0

# error catching

    if (len(x) == 0):
	status = kepmsg.exit('ERROR -- KEPFIT.LSQCLIP: x data array is empty')
    if (len(y) == 0):
	status = kepmsg.exit('ERROR -- KEPFIT.LSQCLIP: y data array is empty')
    if (len(x) < len(pinit)):
	kepmsg.warn(logfile,'WARNING -- KEPFIT.LSQCLIP: no degrees of freedom')
        

# sigma-clipping iterations

    while (iiter < niter and len(x) > len(pinit) and iterstatus > 0):
	iterstatus = 0
	tmpx = []
	tmpy = []
	tmpyerr = []
	npts.append(len(x))
	coeffs,errors,covar,sigma,chi2,dof,fit,plotx,ploty,status = \
	    leastsquare(functype,pinit,x,y,yerr,logfile,verbose)
	pinit = coeffs

# point-by-point sigma-clipping test

	for ix in range(npts[iiter]):
	    if (y[ix] - fit[ix] < rej_hi * sigma and
		fit[ix] - y[ix] < rej_lo * sigma):
		tmpx.append(x[ix])
		tmpy.append(y[ix])
		if (yerr != None): tmpyerr.append(yerr[ix])
	    else:
		iterstatus = 1
	x = scipy.array(tmpx)
	y = scipy.array(tmpy)
	if (yerr != None): yerr = scipy.array(tmpyerr)
	iiter += 1

# fudge scalar models

#    for i in range(len(plotx)):
	

# coeffs = best fit coefficients
# covar = covariance matrix
# iiter = number of sigma clipping iteration before convergence

    return coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status


# linear least square polynomial fit with sigma-clipping
# -----------------------------------------------------------

def poly(x,y,order,rej_lo,rej_hi,niter):

# x = list of x data
# y = list of y data
# order = polynomial order
# rej_lo = lower rejection threshold (units=sigma)
# rej_hi = upper rejection threshold (units=sugma)
# niter = number of sigma-clipping iterations

    import math
    from pylab import polyfit, polyval

    npts = []
    iiter = 0
    iterstatus = 1

# sigma-clipping iterations

    while (iiter < niter and iterstatus > 0):
	iterstatus = 0
	tmpx = []
	tmpy = []
	npts.append(len(x))
	coeffs = polyfit(x,y,order)
	fit = polyval(coeffs,x)

# calculate sigma of fit

	sig = 0
	for ix in range(npts[iiter]):
	    sig = sig + (y[ix] - fit[ix])**2
	sig = math.sqrt(sig / (npts[iiter] - 1))

# point-by-point sigma-clipping test

	for ix in range(npts[iiter]):
	    if (y[ix] - fit[ix] < rej_hi * sig and
		fit[ix] - y[ix] < rej_lo * sig):
		tmpx.append(x[ix])
		tmpy.append(y[ix])
	    else:
		iterstatus = 1
	x = tmpx
	y = tmpy
	iiter += 1

# coeffs = best fit coefficients
# iiter = number of sigma clipping iteration before convergence

    return coeffs, iiter


# Fit single PRF model to Kepler pixel mask data
# -----------------------------------------------------------

def fitPRF(flux,ydim,xdim,column,row,prfn,crval1p,crval2p,cdelt1p,cdelt2p,interpolation,
           tolerance,guess,type,verbose):

# construct input summed image

    status = 0
    if status == 0:
        imgflux = empty((ydim,xdim))
        n = 0
        for i in range(ydim):
            for j in range(xdim):
                imgflux[i,j] = flux[n]
                n += 1

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
            prf = prf / nansum(prf)

# dimensions of data image

    if status == 0:
        datDimY = shape(imgflux)[0]
        datDimX = shape(imgflux)[1]

# dimensions of data image if it had PRF-sized pixels

    if status == 0:
        prfDimY = datDimY / cdelt1p[0]
        prfDimX = datDimX / cdelt2p[0]

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfY0 = (shape(prf)[0] - prfDimY) / 2
        prfX0 = (shape(prf)[1] - prfDimX) / 2

# fit input image with model

    if status == 0:
        args = (imgflux,prf,cdelt1p[0],cdelt2p[0],prfDimY,prfDimX,prfY0,prfX0,
                interpolation,verbose)
        if type == '2D':
            [f,y,x] = fmin_powell(kepfunc.kepler_prf_2d,guess,args=args,xtol=tolerance,
                                  ftol=1.0,disp=False)
        elif type == '1D':
            guess.insert(0,guess[0])
            [fy,fx,y,x] = fmin_powell(kepfunc.kepler_prf_1d,guess,args=args,xtol=tolerance,
                                  ftol=1.0,disp=False)
            f = (fx + fy) / 2 

# calculate best-fit model

    if status == 0:
        prfMod = shift(prf,[y,x],order=1,mode='constant')
        prfMod = prfMod[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX]
        prfFit = rebin2D(prfMod,[numpy.shape(imgflux)[0],numpy.shape(imgflux)[1]],
                                  interpolation,True,False)
        prfFit = prfFit * f / cdelt1p[0] / cdelt2p[0]

# calculate residual between data and model

    if status == 0:
        prfRes = imgflux - prfFit

    return f, y * cdelt1p[0], x * cdelt2p[0], prfMod, prfFit, prfRes


# Fit multi- PRF model to Kepler pixel mask data
# -----------------------------------------------------------

def fitMultiPRF(flux,ydim,xdim,column,row,prfn,crval1p,crval2p,cdelt1p,cdelt2p,interpolation,
                tolerance,ftol,fluxes,columns,rows,type,verbose,logfile):

# caonstruct input summed image

    status = 0
    if status == 0:
        imgflux = empty((ydim,xdim))
        n = 0
        for i in xrange(ydim):
            for j in xrange(xdim):
                imgflux[i,j] = flux[n]
                n += 1

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
            prf = prf / nansum(prf)

# dimensions of data image

    if status == 0:
        datDimY = shape(imgflux)[0]
        datDimX = shape(imgflux)[1]

# dimensions of data image if it had PRF-sized pixels

    if status == 0:
        prfDimY = datDimY / cdelt1p[0]
        prfDimX = datDimX / cdelt2p[0]

# center of the data image (in CCD pixel units)

    if status == 0:
        datCenY = row + float(datDimY) / 2 - 0.5
        datCenX = column + float(datDimX) / 2 - 0.5

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfY0 = (shape(prf)[0] - prfDimY) / 2
        prfX0 = (shape(prf)[1] - prfDimX) / 2

# initial guess for fit parameters

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            y = rows.strip().split(',')
            x = columns.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i]) * numpy.nanmax(flux)
        except:
            f = fluxes
            y = rows
            x = columns
        for i in xrange(len(f)):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(y)):
                try:
                    guess.append((float(y[i]) - datCenY) / cdelt2p[0])
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(x)):
                try:
                    guess.append((float(x[i]) - datCenX) / cdelt1p[0])
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            if len(x) != len(y) or len(x) != len(f):
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        
# fit input image with model

    if status == 0:
        f = []
        x = []
        y = []
        nsrc = len(guess) / 3
        args = (imgflux,prf,cdelt1p[0],cdelt2p[0],prfDimY,prfDimX,prfY0,prfX0,
                interpolation,verbose)
        if type == '2D' and nsrc == 1:
            ans = fmin_powell(kepfunc.kepler_prf_2d,guess,args=args,xtol=tolerance,
                                  ftol=ftol,disp=False)
            f.append(ans[0])
            y.append(ans[1])
            x.append(ans[2]) 
        elif type == '1D' and nsrc == 1:
            guess.insert(0,guess[0])
            ans = fmin_powell(kepfunc.kepler_prf_1d,guess,args=args,xtol=tolerance,
                                  ftol=ftol,disp=False)
            f.append((ans[0] + ans[1]) / 2)
            y.append(ans[2])
            x.append(ans[3]) 
        else:
            ans = fmin_powell(kepfunc.kepler_multi_prf_2d,guess,args=args,xtol=tolerance,
                              ftol=ftol,disp=False)
            for i in xrange(nsrc):
                f.append(ans[i])
                y.append(ans[nsrc+i])
                x.append(ans[nsrc*2+i]) 

# calculate best-fit model

    if status == 0:
        prfMod = numpy.zeros((prfDimY+1,prfDimX+1))
        for i in xrange(nsrc):
            prfTmp = shift(prf,[y[i],x[i]],order=1,mode='constant')
            prfTmp = prfTmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX]
            prfMod = prfMod + prfTmp * f[i]
        prfFit = rebin2D(prfMod,[numpy.shape(imgflux)[0],numpy.shape(imgflux)[1]],
                         interpolation,True,False) / cdelt1p[0] / cdelt2p[0]

# calculate residual between data and model

    if status == 0:
        prfRes = imgflux - prfFit

# convert PRF pixels sizes to CCD pixel sizes

    if status == 0:
        for i in xrange(nsrc):
            y[i] = y[i] * cdelt1p[0] + datCenY
            x[i] = x[i] * cdelt2p[0] + datCenX

    return f, y, x, prfMod, prfFit, prfRes


# Fit multi- PRF model + constant background to Kepler pixel mask data
# -----------------------------------------------------------

def fitBackMultiPRF(flux,ydim,xdim,column,row,prfn,crval1p,crval2p,cdelt1p,cdelt2p,interpolation,
                tolerance,ftol,fluxes,columns,rows,bkg,type,verbose,logfile):

# caonstruct input summed image

    status = 0
    if status == 0:
        imgflux = empty((ydim,xdim))
        n = 0
        for i in xrange(ydim):
            for j in xrange(xdim):
                imgflux[i,j] = flux[n]
                n += 1

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
            prf = prf / nansum(prf)

# dimensions of data image

    if status == 0:
        datDimY = shape(imgflux)[0]
        datDimX = shape(imgflux)[1]

# dimensions of data image if it had PRF-sized pixels

    if status == 0:
        prfDimY = datDimY / cdelt1p[0]
        prfDimX = datDimX / cdelt2p[0]

# center of the data image (in CCD pixel units)

    if status == 0:
        datCenY = row + float(datDimY) / 2 - 0.5
        datCenX = column + float(datDimX) / 2 - 0.5

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfY0 = (shape(prf)[0] - prfDimY) / 2
        prfX0 = (shape(prf)[1] - prfDimX) / 2

# initial guess for fit parameters

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            y = rows.strip().split(',')
            x = columns.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i]) * numpy.nanmax(flux)
        except:
            f = fluxes
            y = rows
            x = columns
        b = bkg
        for i in xrange(len(f)):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(y)):
                try:
                    guess.append((float(y[i]) - datCenY) / cdelt2p[0])
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(x)):
                try:
                    guess.append((float(x[i]) - datCenX) / cdelt1p[0])
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            guess.append(b)
        if status == 0:
            if len(x) != len(y) or len(x) != len(f):
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        
# fit input image with model

    if status == 0:
        f = []
        x = []
        y = []
        nsrc = (len(guess) - 1) / 3
        args = (imgflux,prf,cdelt1p[0],cdelt2p[0],prfDimY,prfDimX,prfY0,prfX0,
                interpolation,verbose)
        ans = fmin_powell(kepfunc.kepler_bkg_multi_prf_2d,guess,args=args,xtol=tolerance,
                          ftol=ftol,disp=False)
        for i in xrange(nsrc):
            f.append(ans[i])
            y.append(ans[nsrc+i])
            x.append(ans[nsrc*2+i]) 
            b = ans[nsrc*3]

# calculate best-fit model

    if status == 0:
        prfMod = numpy.zeros((prfDimY+1,prfDimX+1))
        for i in xrange(nsrc):
            prfTmp = shift(prf,[y[i],x[i]],order=1,mode='constant')
            prfTmp = prfTmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX]
            prfMod = prfMod + prfTmp * f[i]
        prfFit = rebin2D(prfMod,[numpy.shape(imgflux)[0],numpy.shape(imgflux)[1]],
                         interpolation,True,False) / cdelt1p[0] / cdelt2p[0]
        prfFit = prfFit + b
        
# calculate residual between data and model

    if status == 0:
        prfRes = imgflux - prfFit

# convert PRF pixels sizes to CCD pixel sizes

    if status == 0:
        for i in xrange(nsrc):
            y[i] = y[i] * cdelt1p[0] + datCenY
            x[i] = x[i] * cdelt2p[0] + datCenX

    return f, y, x, b, prfMod, prfFit, prfRes


# Fit multi- PRF model + constant background with focus variations to Kepler pixel mask data
# ------------------------------------------------------------------------------------------

def fitFocusMultiPRF(flux,ydim,xdim,column,row,prfn,crval1p,crval2p,cdelt1p,cdelt2p,interpolation,
                tolerance,ftol,fluxes,columns,rows,bkg,wfac,type,verbose,logfile):

# caonstruct input summed image

    status = 0
    if status == 0:
        imgflux = empty((ydim,xdim))
        n = 0
        for i in xrange(ydim):
            for j in xrange(xdim):
                imgflux[i,j] = flux[n]
                n += 1

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
            prf = prf / nansum(prf)

# dimensions of data image

    if status == 0:
        datDimY = shape(imgflux)[0]
        datDimX = shape(imgflux)[1]

# dimensions of data image if it had PRF-sized pixels

#    if status == 0:
#        prfDimY = datDimY / cdelt1p[0]
#        prfDimX = datDimX / cdelt2p[0]

# center of the data image (in CCD pixel units)

    if status == 0:
        datCenY = row + float(datDimY) / 2 - 0.5
        datCenX = column + float(datDimX) / 2 - 0.5

# location of the data image centered on the PRF image (in PRF pixel units)

#    if status == 0:
#        prfY0 = (shape(prf)[0] - prfDimY) / 2
#        prfX0 = (shape(prf)[1] - prfDimX) / 2

# initial guess for fit parameters

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            y = rows.strip().split(',')
            x = columns.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i]) * numpy.nanmax(flux)
        except:
            f = fluxes
            y = rows
            x = columns
        b = bkg
        w = wfac
        for i in xrange(len(f)):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(y)):
                try:
                    guess.append((float(y[i]) - datCenY) / cdelt2p[0])
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(x)):
                try:
                    guess.append((float(x[i]) - datCenX) / cdelt1p[0])
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            guess.append(b)
            guess.append(w)
        if status == 0:
            if len(x) != len(y) or len(x) != len(f):
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        
# fit input image with model

    if status == 0:
        f = []
        x = []
        y = []
        nsrc = (len(guess) - 2) / 3
        args = (imgflux,prf,cdelt1p[0],cdelt2p[0],datDimY,datDimX,interpolation,verbose)
        ans = fmin_powell(kepfunc.kepler_focus_multi_prf_2d,guess,args=args,xtol=tolerance,
                          ftol=ftol,disp=False)
        for i in xrange(nsrc):
            f.append(ans[i])
            y.append(ans[nsrc+i])
            x.append(ans[nsrc*2+i]) 
        b = ans[nsrc*3]
        w = ans[nsrc*3+1]
        print ans
        print f,y,x,b,w

# calculate best-fit model

    if status == 0:
        prfDimY = datDimY / cdelt1p[0] / w
        prfDimX = datDimX / cdelt2p[0] / w
        prfY0 = (shape(prf)[0] - prfDimY) / 2
        prfX0 = (shape(prf)[1] - prfDimX) / 2
        DY = 0.0; DX = 0.0
        if int(prfDimY) % 2 == 0: DY = 1.0
        if int(prfDimX) % 2 == 0: DX = 1.0
        print w, prfDimY, prfDimX
        prfMod = zeros((prfDimY+DY,prfDimX+DX))
        for i in range(nsrc):
            prfTmp = shift(prf,[y[i]/w,x[i]/w],order=1,mode='constant')
            prfMod = prfMod + prfTmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX] * f[i]
        prfFit = rebin2D(prfMod,[shape(imgflux)[0],shape(imgflux)[1]],interpolation,True,False)
        prfFit = prfFit / cdelt1p[0] / cdelt2p[0] / w / w
        prfFit = prfFit + b
        
# calculate residual between data and model

    if status == 0:
        prfRes = imgflux - prfFit

# convert PRF pixels sizes to CCD pixel sizes

    if status == 0:
        for i in xrange(nsrc):
            y[i] = y[i] * cdelt1p[0] * w + datCenY
            x[i] = x[i] * cdelt2p[0] * w + datCenX

    return f, y, x, b, w, prfMod, prfFit, prfRes



# Fit multi-PRF model to Kepler pixel mask data
# -----------------------------------------------------------

def test(flux,ydim,xdim,column,row,prfn,crval1p,crval2p,cdelt1p,cdelt2p,interpolation,
                tolerance,ftol,fluxes,columns,rows,type,verbose,logfile):

# construct input summed image

    status = 0
    if status == 0:
        imgflux = empty((ydim,xdim))
        n = 0
        for i in xrange(ydim):
            for j in xrange(xdim):
                imgflux[i,j] = flux[n]
                n += 1

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
            prf = prf / nansum(prf)

# dimensions of data image

    if status == 0:
        datDimY = shape(imgflux)[0]
        datDimX = shape(imgflux)[1]

# dimensions of data image if it had PRF-sized pixels

    if status == 0:
        prfDimY = datDimY / cdelt1p[0]
        prfDimX = datDimX / cdelt2p[0]

# center of the data image (in CCD pixel units)

    if status == 0:
        datCenY = row + float(datDimY) / 2 - 0.5
        datCenX = column + float(datDimX) / 2 - 0.5

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfY0 = (shape(prf)[0] - prfDimY) / 2
        prfX0 = (shape(prf)[1] - prfDimX) / 2

# initial guess for fit parameters

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            y = rows.strip().split(',')
            x = columns.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i]) * numpy.nanmax(flux)
        except:
            f = fluxes
            y = rows
            x = columns
        for i in xrange(len(f)):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(y)):
                try:
                    guess.append((float(y[i]) - datCenY) / cdelt2p[0])
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(len(x)):
                try:
                    guess.append((float(x[i]) - datCenX) / cdelt1p[0])
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            if len(x) != len(y) or len(x) != len(f):
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        
# fit input image with model

    if status == 0:
        f = []
        x = []
        y = []
        nsrc = len(guess) / 3
        args = (imgflux,prf,cdelt1p[0],cdelt2p[0],prfDimY,prfDimX,prfY0,prfX0,
                interpolation,verbose)
        if type == '2D' and nsrc == 1:
            ans = fmin_powell(kepfunc.kepler_prf_2d,guess,args=args,xtol=tolerance,
                                  ftol=ftol,disp=False)
            f.append(ans[0])
            y.append(ans[1])
            x.append(ans[2]) 
        elif type == '1D' and nsrc == 1:
            guess.insert(0,guess[0])
            ans = fmin_powell(kepfunc.kepler_prf_1d,guess,args=args,xtol=tolerance,
                                  ftol=ftol,disp=False)
            f.append((ans[0] + ans[1]) / 2)
            y.append(ans[2])
            x.append(ans[3]) 
        else:
            ans = fmin_powell(kepfunc.kepler_multi_prf_2d,guess,args=args,xtol=tolerance,
                              ftol=ftol,disp=False)
            for i in xrange(nsrc):
                f.append(ans[i])
                y.append(ans[nsrc+i])
                x.append(ans[nsrc*2+i]) 

# calculate best-fit model

    if status == 0:
        prfMod = numpy.zeros((prfDimY+1,prfDimX+1))
        for i in xrange(nsrc):
            prfTmp = shift(prf,[y[i],x[i]],order=1,mode='constant')
            prfTmp = prfTmp[prfY0:prfY0+prfDimY,prfX0:prfX0+prfDimX]
            prfMod = prfMod + prfTmp * f[i]
        prfFit = rebin2D(prfMod,[numpy.shape(imgflux)[0],numpy.shape(imgflux)[1]],
                         interpolation,True,False) / cdelt1p[0] / cdelt2p[0]

# calculate residual between data and model

    if status == 0:
        prfRes = imgflux - prfFit

# convert PRF pixels sizes to CCD pixel sizes

    if status == 0:
        for i in xrange(nsrc):
            y[i] = y[i] * cdelt1p[0] + datCenY
            x[i] = x[i] * cdelt2p[0] + datCenX

    return f, y, x, prfMod, prfFit, prfRes


