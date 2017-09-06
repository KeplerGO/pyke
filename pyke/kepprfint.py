import os
import glob
import math
import numpy as np
from scipy.interpolate import RectBivariateSpline

from . import kepio, kepmsg


all = ['read_and_interpolate_prf']


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
    for i in range(n_hdu):
        prfWeight[i] = math.sqrt(
            (column - crval1p[i])**2 + (row - crval2p[i])**2)
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

