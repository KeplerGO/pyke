from . import kepmsg, kepstat, kepfunc, keparray
import math
import numpy as np
from scipy import optimize, ndimage
from scipy.optimize import fmin_powell, fmin_tnc, fmin, leastsq
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift


def leastsquares(fitfunc, pinit, xdata, ydata, yerr, logfile, verbose):
    """linear least square polynomial fit using scipy"""
    errfunc = lambda p, x, y, err: np.sum(((y - fitfunc(p, x)) / err) ** 2)

    if yerr is None:
        yerr = np.ones(len(ydata))
    # fit data
    try:
        out = optimize.minimize(errfunc, pinit, args=(xdata, ydata, yerr))
    except:
        errmsg = 'ERROR -- KEPFIT.LEASTSQUARES: failed to fit data'
        kepmsg.err(logfile, errmsg, verbose)

    coeffs = out.x
    covar = out.hess_inv

    # calculate 1-sigma error on coefficients
    if covar is None:
        errmsg = 'ERROR -- KEPFIT.leastsquares: NULL covariance matrix'
        kepmsg.err(logfile, errmsg, verbose)
    else:
        if len(coeffs) > 1:
            errors = np.sqrt(np.diag(covar))
        else:
            errors = math.sqrt(covar)

    # generate fit points for rms calculation
    fit = fitfunc(coeffs, xdata)
    sigma = kepstat.rms(ydata, fit, logfile, verbose)

    # generate fit points for plotting
    plotx = np.linspace(xdata.min(), xdata.max(), 10000)
    ploty = fitfunc(coeffs, plotx)

    # reduced chi^2 calculation
    chi2 = 0
    dof = len(ydata) - len(coeffs)
    for i in range(len(ydata)):
        chi2 += (ydata[i] - fit[i]) ** 2 / yerr[i]
    chi2 /= dof

    return coeffs, errors, covar, sigma, chi2, dof, fit, plotx, ploty


def lsqclip(fit_func, pinit, x, y, yerr, rej_lo, rej_hi, niter, logfile,
            verbose):
    """linear least square fit with sigma-clipping"""

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

    # error catching
    if len(x) == 0:
        kepmsg.exit('ERROR -- KEPFIT.LSQCLIP: x data array is empty')
    if len(y) == 0:
        kepmsg.exit('ERROR -- KEPFIT.LSQCLIP: y data array is empty')
    if len(x) < len(pinit):
        kepmsg.warn(logfile, ("WARNING -- KEPFIT.LSQCLIP: no degrees of"
                              "freedom"), verbose)

    # sigma-clipping iterations
    while (iiter < niter and len(x) > len(pinit) and iterstatus > 0):
        iterstatus = 0
        tmpx = []
        tmpy = []
        tmpyerr = []
        npts.append(len(x))
        coeffs, errors, covar, sigma, chi2, dof, fit, plotx, ploty = \
            leastsquares(fit_func, pinit, x, y, yerr, logfile, verbose)
        pinit = coeffs

        # point-by-point sigma-clipping test
        for ix in range(npts[iiter]):
            if (y[ix] - fit[ix] < rej_hi * sigma and
                fit[ix] - y[ix] < rej_lo * sigma):
                tmpx.append(x[ix])
                tmpy.append(y[ix])
                if yerr is not None:
                    tmpyerr.append(yerr[ix])
            else:
                iterstatus = 1
        x = np.array(tmpx)
        y = np.array(tmpy)
        if yerr is not None:
            yerr = np.array(tmpyerr)
        iiter += 1

    # coeffs = best fit coefficients
    # covar = covariance matrix
    # iiter = number of sigma clipping iteration before convergence
    return coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty

def poly(x, y, order, rej_lo, rej_hi, niter):
    """linear least square polynomial fit with sigma-clipping"""

    # x = list of x data
    # y = list of y data
    # order = polynomial order
    # rej_lo = lower rejection threshold (units=sigma)
    # rej_hi = upper rejection threshold (units=sugma)
    # niter = number of sigma-clipping iterations

    npts = []
    iiter = 0
    iterstatus = 1

    # sigma-clipping iterations
    while (iiter < niter and iterstatus > 0):
        iterstatus = 0
        tmpx = []
        tmpy = []
        npts.append(len(x))
        coeffs = np.polyfit(x, y, order)
        fit = np.polyval(coeffs, x)

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

def fitPRF(flux, ydim, xdim, column, row, prfn, crval1p, crval2p, cdelt1p,
           cdelt2p, interpolation, tol, guess, mode, verbose):
    """Fit single PRF model to Kepler pixel mask data"""

    # construct input summed image
    imgflux = np.empty((ydim, xdim))
    n = 0
    for i in range(ydim):
        for j in range(xdim):
            imgflux[i, j] = flux[n]
            n += 1

    # interpolate the calibrated PRF shape to the target position
    prf = np.zeros(np.shape(prfn[0]), dtype='float32')
    prfWeight = np.zeros((5), dtype='float32')
    for i in range(5):
        prfWeight[i] = math.sqrt((column - crval1p[i]) ** 2
                                  + (row - crval2p[i]) ** 2)
        if prfWeight[i] == 0.0:
            prfWeight[i] = 1.0e6
        prf = prf + prfn[i] / prfWeight[i]
        prf = prf / np.nansum(prf)

    # dimensions of data image
    datDimY = np.shape(imgflux)[0]
    datDimX = np.shape(imgflux)[1]

    # dimensions of data image if it had PRF-sized pixels
    prfDimY = datDimY / cdelt1p[0]
    prfDimX = datDimX / cdelt2p[0]

    # location of the data image centered on the PRF image (in PRF pixel units)
    prfY0 = (np.shape(prf)[0] - prfDimY) / 2
    prfX0 = (np.shape(prf)[1] - prfDimX) / 2

    # fit input image with model
    args = (imgflux, prf, cdelt1p[0], cdelt2p[0], prfDimY, prfDimX, prfY0,
            prfX0, interpolation, verbose)
    if mode == '2D':
        [f, y, x] = fmin_powell(kepfunc.kepler_prf_2d, guess, args=args,
                                xtol=tol, ftol=1.0, disp=False)
    elif mode == '1D':
        guess.insert(0, guess[0])
        [fy, fx, y, x] = fmin_powell(kepfunc.kepler_prf_1d, guess, args=args,
                                     xtol=tol, ftol=1.0,disp=False)
        f = (fx + fy) / 2.0

    # calculate best-fit model
    prfMod = shift(prf, [y, x], order=1, mode='constant')
    prfMod = prfMod[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX]
    prfFit = keparray.rebin2D(prfMod, [np.shape(imgflux)[0],
                              np.shape(imgflux)[1]], interpolation, True,
                              False)
    prfFit = prfFit * f / cdelt1p[0] / cdelt2p[0]

    # calculate residual between data and model
    prfRes = imgflux - prfFit

    return f, y * cdelt1p[0], x * cdelt2p[0], prfMod, prfFit, prfRes

def fitMultiPRF(flux, ydim, xdim, column, row, prfn, crval1p, crval2p,
                cdelt1p, cdelt2p, interpolation, tol, ftol, fluxes,
                columns, rows, mode, verbose, logfile):
    """Fit multi-PRF model to Kepler pixel mask data"""

    # construct input summed image
    imgflux = np.empty((ydim,xdim))
    n = 0
    for i in range(ydim):
        for j in range(xdim):
            imgflux[i, j] = flux[n]
            n += 1

    # interpolate the calibrated PRF shape to the target position
    prf = np.zeros(np.shape(prfn[0]), dtype='float32')
    prfWeight = np.zeros((5), dtype='float32')
    for i in range(5):
        prfWeight[i] = math.sqrt((column - crval1p[i]) ** 2
                                 + (row - crval2p[i]) ** 2)
        if prfWeight[i] == 0.0:
            prfWeight[i] = 1.0e6
        prf = prf + prfn[i] / prfWeight[i]
        prf = prf / np.nansum(prf)

    # dimensions of data image
    datDimY = np.shape(imgflux)[0]
    datDimX = np.shape(imgflux)[1]

    # dimensions of data image if it had PRF-sized pixels
    prfDimY = datDimY / cdelt1p[0]
    prfDimX = datDimX / cdelt2p[0]

    # center of the data image (in CCD pixel units)
    datCenY = row + float(datDimY) / 2 - 0.5
    datCenX = column + float(datDimX) / 2 - 0.5

    # location of the data image centered on the PRF image (in PRF pixel units)
    prfY0 = (np.shape(prf)[0] - prfDimY) / 2
    prfX0 = (np.shape(prf)[1] - prfDimX) / 2

    # initial guess for fit parameters
    guess = []
    if len(x) != len(y) or len(x) != len(f):
        errmsg = ('ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                  'fluxes must have the same number of sources')
        kepmsg.err(logfile, message, verbose)
    else:
        for i in range(len(fluxes)):
            guess.append(float(fluxes[i]))
            guess.append((float(rows[i]) - datCenY) / cdelt2p[0])
            guess.append((float(columns[i]) - datCenX) / cdelt1p[0])

    # fit input image with model
    f, x, y = [], [], []
    nsrc = len(guess) // 3
    args = (imgflux, prf, cdelt1p[0], cdelt2p[0], prfDimY, prfDimX, prfY0,
            prfX0, interpolation, verbose)
    if mode == '2D' and nsrc == 1:
        ans = fmin_powell(kepfunc.kepler_prf_2d, guess, args=args, xtol=tol,
                          ftol=ftol, disp=False)
        f.append(ans[0])
        y.append(ans[1])
        x.append(ans[2])
    elif mode == '1D' and nsrc == 1:
        guess.insert(0,guess[0])
        ans = fmin_powell(kepfunc.kepler_prf_1d, guess, args=args, xtol=tol,
                          ftol=ftol, disp=False)
        f.append((ans[0] + ans[1]) / 2)
        y.append(ans[2])
        x.append(ans[3])
    else:
        ans = fmin_powell(kepfunc.kepler_multi_prf_2d, guess, args=args,
                          xtol=tol, ftol=ftol, disp=False)
        for i in range(nsrc):
            f.append(ans[i])
            y.append(ans[nsrc + i])
            x.append(ans[nsrc * 2 + i])

    # calculate best-fit model
    prfMod = np.zeros((prfDimY + 1, prfDimX + 1))
    for i in range(nsrc):
        prfTmp = shift(prf, [y[i], x[i]], order=1, mode='constant')
        prfTmp = prfTmp[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX]
        prfMod = prfMod + prfTmp * f[i]
    prfFit = keparray.rebin2D(prfMod, [np.shape(imgflux)[0],
                              np.shape(imgflux)[1]], interpolation, True,
                              False) / cdelt1p[0] / cdelt2p[0]

    prfRes = imgflux - prfFit

    # convert PRF pixels sizes to CCD pixel sizes
    for i in range(nsrc):
        y[i] = y[i] * cdelt1p[0] + datCenY
        x[i] = x[i] * cdelt2p[0] + datCenX

    return f, y, x, prfMod, prfFit, prfRes

def fitBackMultiPRF(flux, ydim, xdim, column, row, prfn, crval1p, crval2p,
                    cdelt1p, cdelt2p, interpolation, tol, ftol, fluxes,
                    columns, rows, bkg, mode, verbose, logfile):
    """Fit multi- PRF model + constant background to Kepler pixel mask data"""

    # construct input summed image
    imgflux = np.asarray(flux).reshape(ydim, xdim)

    # interpolate the calibrated PRF shape to the target position
    prf = np.zeros(np.shape(prfn[0]), dtype='float32')
    prfWeight = np.zeros((5), dtype='float32')
    for i in range(5):
        prfWeight[i] = math.sqrt((column - crval1p[i])**2
                                 + (row - crval2p[i])**2)
        if prfWeight[i] == 0.0:
            prfWeight[i] = 1.0e6
        prf = prf + prfn[i] / prfWeight[i]
        prf = prf / np.nansum(prf)

    # dimensions of data image
    datDimY = np.shape(imgflux)[0]
    datDimX = np.shape(imgflux)[1]

    # dimensions of data image if it had PRF-sized pixels
    prfDimY = datDimY / cdelt1p[0]
    prfDimX = datDimX / cdelt2p[0]

    # center of the data image (in CCD pixel units)
    datCenY = row + float(datDimY) / 2 - 0.5
    datCenX = column + float(datDimX) / 2 - 0.5

    # location of the data image centered on the PRF image (in PRF pixel units)
    prfY0 = (np.shape(prf)[0] - prfDimY) / 2
    prfX0 = (np.shape(prf)[1] - prfDimX) / 2

    # initial guess for fit parameters
    guess = []
    f, y, x = fluxes, rows, columns
    if len(x) != len(y) or len(x) != len(f):
        errmsg = ('ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                  'fluxes must have the same number of sources')
        kepmsg.err(logfile, message, verbose)
    else:
        for i in range(len(fluxes)):
            guess.append(float(fluxes[i]))
            guess.append((float(rows[i]) - datCenY) / cdelt2p[0])
            guess.append((float(columns[i]) - datCenX) / cdelt1p[0])
    guess.append(bkg)

    # fit input image with model
    f, x, y = [], [], []
    nsrc = (len(guess) - 1) // 3
    args = (imgflux, prf, cdelt1p[0], cdelt2p[0], prfDimY, prfDimX, prfY0,
            prfX0, interpolation, verbose)
    ans = fmin_powell(kepfunc.kepler_bkg_multi_prf_2d, guess, args=args,
                      xtol=tol, ftol=ftol, disp=False)
    for i in range(nsrc):
        f.append(ans[i])
        y.append(ans[nsrc + i])
        x.append(ans[nsrc * 2 + i])
        b = ans[nsrc * 3]

    # calculate best-fit model
    prfMod = np.zeros((prfDimY + 1, prfDimX + 1))
    for i in range(nsrc):
        prfTmp = shift(prf, [y[i], x[i]], order=1, mode='constant')
        prfTmp = prfTmp[prfY0:prfY0 + prfDimY, prfX0:prfX0 + prfDimX]
        prfMod = prfMod + prfTmp * f[i]
    prfFit = keparray.rebin2D(prfMod, [np.shape(imgflux)[0],
                              np.shape(imgflux)[1]], interpolation, True,
                              False) / cdelt1p[0] / cdelt2p[0]
    prfFit = prfFit + b

    # calculate residual between data and model
    prfRes = imgflux - prfFit

    # convert PRF pixels sizes to CCD pixel sizes
    for i in range(nsrc):
        y[i] = y[i] * cdelt1p[0] + datCenY
        x[i] = x[i] * cdelt2p[0] + datCenX

    return f, y, x, b, prfMod, prfFit, prfRes

def fitFocusMultiPRF(flux, ydim, xdim, column, row, prfn, crval1p, crval2p,
                     cdelt1p, cdelt2p, interpolation, tol, ftol, fluxes,
                     columns, rows, bkg, wfac, verbose, logfile):
    """Fit multi- PRF model + constant background with focus variations to
    Kepler pixel mask data
    """

    # caonstruct input summed image
    imgflux = np.asarray(flux).reshape(ydim, xdim)
    # interpolate the calibrated PRF shape to the target position

    prf = np.zeros(np.shape(prfn[0]), dtype='float32')
    prfWeight = np.zeros((5), dtype='float32')
    for i in range(5):
        prfWeight[i] = math.sqrt((column - crval1p[i]) ** 2
                                 + (row - crval2p[i]) ** 2)
        if prfWeight[i] == 0.0:
            prfWeight[i] = 1.0e6
        prf = prf + prfn[i] / prfWeight[i]
        prf = prf / np.nansum(prf)

    # dimensions of data image
    datDimY = np.shape(imgflux)[0]
    datDimX = np.shape(imgflux)[1]

    # center of the data image (in CCD pixel units)
    datCenY = row + float(datDimY) / 2 - 0.5
    datCenX = column + float(datDimX) / 2 - 0.5

    # initial guess for fit parameters
    guess = []
    f, y, x = fluxes, rows, columns
    b, w = bkg, wfac

    # initial guess for fit parameters
    guess = []
    if len(x) != len(y) or len(x) != len(f):
        errmsg = ('ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                  'fluxes must have the same number of sources')
        kepmsg.err(logfile, message, verbose)
    else:
        for i in range(len(fluxes)):
            guess.append(float(fluxes[i]))
            guess.append((float(rows[i]) - datCenY) / cdelt2p[0])
            guess.append((float(columns[i]) - datCenX) / cdelt1p[0])
    guess.append(b)
    guess.append(w)

    # fit input image with model
    f, x, y = [], [], []
    nsrc = (len(guess) - 2) / 3
    args = (imgflux, prf, cdelt1p[0], cdelt2p[0], datDimY, datDimX,
            interpolation, verbose)
    ans = fmin_powell(kepfunc.kepler_focus_multi_prf_2d, guess, args=args,
                      xtol=tol, ftol=ftol, disp=False)
    for i in range(nsrc):
        f.append(ans[i])
        y.append(ans[nsrc + i])
        x.append(ans[nsrc * 2 + i])
    b = ans[nsrc * 3]
    w = ans[nsrc * 3 + 1]
    print(ans)
    print(f, y, x, b, w)

    # calculate best-fit model
    prfDimY = datDimY / cdelt1p[0] / w
    prfDimX = datDimX / cdelt2p[0] / w
    prfY0 = (np.shape(prf)[0] - prfDimY) / 2
    prfX0 = (np.shape(prf)[1] - prfDimX) / 2
    DY, DX = 0.0, 0.0
    if int(prfDimY) % 2 == 0:
        DY = 1.0
    if int(prfDimX) % 2 == 0:
        DX = 1.0
    print(w, prfDimY, prfDimX)
    prfMod = np.zeros((prfDimY + DY, prfDimX + DX))
    for i in range(nsrc):
        prfTmp = shift(prf, [y[i] / w, x[i] / w], order=1, mode='constant')
        prfMod = (prfMod + prfTmp[prfY0:prfY0 + prfDimY,
                                  prfX0:prfX0 + prfDimX] * f[i])
    prfFit = keparray.rebin2D(prfMod, [np.shape(imgflux)[0],
                              np.shape(imgflux)[1]], interpolation, True,
                              False)
    prfFit = prfFit / cdelt1p[0] / cdelt2p[0] / w / w
    prfFit = prfFit + b

    # calculate residual between data and model
    prfRes = imgflux - prfFit

    # convert PRF pixels sizes to CCD pixel sizes
    for i in range(nsrc):
        y[i] = y[i] * cdelt1p[0] * w + datCenY
        x[i] = x[i] * cdelt2p[0] * w + datCenX

    return f, y, x, b, w, prfMod, prfFit, prfRes
