from .utils import PyKEArgumentHelpFormatter
from . import kepmsg, kepio, kepkey
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cbook import is_numlike
from scipy.optimize import leastsq
from scipy.optimize import fmin
from scipy.interpolate import interp1d
from astropy.io import fits as pyfits
from tqdm import tqdm
import re

__all__ = ['kepcotrend']


def cut_bad_data(cad, date, flux, err):
    """
    this function finds cadences with good data and returns them
    """

    good_data_mask = np.logical_and(np.isfinite(date), np.isfinite(flux))
    date = date[good_data_mask]
    cad = cad[good_data_mask]
    flux = flux[good_data_mask]
    err = err[good_data_mask]

    return cad, date, flux, err, good_data_mask

def put_in_nans(good_data, flux):
    """
    Function finds the cadences where the data has been removed using
    cut_bad_data() and puts data back in. The flux data put back in is nan.
    This function is used when writing data to a FITS files.
    good_data == True means the datapoint is good!!
    """
    newflux = np.empty(len(good_data))
    newflux[:] = np.nan
    newflux[good_data] = flux

    return newflux

def get_pcomp_list_newformat(bvdat, pcomplist, newcad, short, scinterp):
    """
    Finds cotrending basis vectors which have been requested to be
    used by the user and adds them to an array.
    """
    pcomp = np.zeros((len(pcomplist), len(newcad)))
    for i in range(len(np.array(pcomplist))):
        j = int(np.array(pcomplist)[i])
        dat = bvdat.field('VECTOR_{}'.format(j))[~np.isnan(bvdat.field('CADENCENO'))]
        bvcadfull = bvdat.field('CADENCENO')[~np.isnan(bvdat.field('CADENCENO'))]
        #try:
        if short:
            #if the data is short cadence the interpolate the basis vectors
            bv_data = dat[np.in1d(bvdat.field('CADENCENO'), newcad)]
            bv_cad = bvcadfull[np.in1d(bvdat.field('CADENCENO'), newcad)]
            #funny things happen why I use interp1d for linear interpolation
            #so I have opted to use the numpy interp function for linear
            if scinterp == 'linear':
                intpl = np.interp(newcad, bv_cad, bv_data, left=bv_data[0],
                                  right=bv_data[-1])
                pcomp[i] = intpl
            else:
                intpl = interp1d(bv_cad, bv_data, kind=scinterp,
                                 bounds_error=False, fill_value=None)
                pcomp[i] = np.where(np.isnan(intpl(newcad)), 0, intpl(newcad))
                mid_pt = np.floor(np.median(np.arange(len(pcomp[i]))))
                p_len = len(pcomp[i])
                lower = [np.logical_and(np.arange(p_len) < mid_pt, pcomp[i] == 0)]
                upper = [np.logical_and(np.arange(p_len) > mid_pt, pcomp[i] == 0)]
                pcomp[i][lower] = bv_data[0]
                pcomp[i][upper] = bv_data[-1]
        else:
            pcomp[i] = dat[np.in1d(bvdat.field('CADENCENO'), newcad)]
    return pcomp

def make_sc_lc(obs_cad, bv_cad, flux):
    """
    make short cadence data look like long cadence data
    """
    newflux = np.zeros(len(bv_cad))
    for i in range(len(bv_cad)):
        mask = np.logical_and(obs_cad > bv_cad[i] - 15,
                              obs_cad < bv_cad[i] + 15)
        newflux[i] = obs_cad[mask]
    return newflux

def near_intpl(xout, xin, yin):
    """
    Interpolate the curve defined by (xin, yin) at points xout. The array
    xin must be monotonically increasing. The output has the same data type as
    the input yin.

    :param yin: y values of input curve
    :param xin: x values of input curve
    :param xout: x values of output interpolated curve
    :param method: interpolation method ('linear' | 'nearest')

    @:rtype: numpy array with interpolated curve
    """

    lenxin = len(xin)
    i1 = searchsorted(xin, xout)
    i1[i1 == 0] = 1
    i1[i1 == lenxin] = lenxin - 1
    x0 = xin[i1 - 1]
    x1 = xin[i1]
    y0 = yin[i1 - 1]
    y1 = yin[i1]

    return np.where(abs(xout - x0) < abs(xout - x1), y0, y1)

def get_pcomp_list(pcompdata, pcomplist, newcad):
    pcomp = np.zeros((len(pcomplist), len(newcad)))
    for i in range(len(np.array(pcomplist))):
        j = int(np.array(pcomplist)[i]) - 1
        dat = pcompdata[..., j + 2]
        pcomp[i] = dat[np.in1d(pcompdata[..., 1], newcad)]
    return pcomp

def do_lsq_uhat(pcomps, flux):
    """
    does a linear least squares fit of the basis vectors to the light curve
    using the 'matrix' method - U(transpose) * y = coeffs
    In my implimentation y is a horizontal 1D array and U is also a long thin
    array of length correpsonding to the number of basis vectors use. In
    effect what I have is my U is already transposed and y need to be
    transposed to used. First I convert to what is expected in leasts
    squares fitting
    """

    U_hat = np.matrix(pcomps).transpose()
    y_hat = np.matrix(flux).transpose()
    U_trans = U_hat.transpose()
    coeffs = - np.linalg.inv(U_trans * U_hat) * U_trans * y_hat

    return coeffs

def do_lsq_nlin(pcomps, flux):
    """
    does a linear least squares fit of the basis vectors to the light curve
    using the 'lst_sq' method - this performs a Levenberg-Marquart least
    squares fit. The initial guess is an array of zeros
    """

    guess = np.append(np.array([1.]), np.zeros(len(pcomps) - 1))
    t = leastsq(fitfunct, guess, args=(pcomps, flux), full_output=0)
    return - np.array(t[0])

def do_lsq_fmin(pcomps, flux):
    """
    performs a simplex fit of the basis vectors to the light curve.
    Initial guess is an array with 1. as the first element and zero as the
    value of all other elements in the array
    """

    guess = np.append(np.array([1.]), np.zeros(len(pcomps) - 1))
    t = fmin(fitfunct_fmin, guess, args=(pcomps, flux))
    return -np.array(t)

def do_lsq_fmin_pow(pcomps, flux, order):
    """
    performs a simplex fit of the basis vectors to the light curve.
    Initial guess is an array with 1. as the first element and zero as the
    value of all other elements in the array
    """

    guess = np.array([1, 0])
    initial = fmin(fitfunct_fmin_pow, guess, args=(pcomps[0:2], flux, order))
    guess = np.append(initial, np.zeros(len(pcomps) - 2))
    t = fmin(fitfunct_fmin_pow, guess, args=(pcomps, flux, order))
    return - np.array(t)

def fitfunct_fmin(scale, pcomp, zeroflux):
    outflux = fitfunct(scale, pcomp, zeroflux)
    sumsq = np.sum(np.abs(outflux))
    return sumsq

def fitfunct_fmin_pow(scale, pcomp, zeroflux, order):
    outflux = fitfunct(scale, pcomp, zeroflux)
    sumsq = np.sum(np.power(np.abs(outflux), order))
    return sumsq

def fitfunct(scale, pcomp, zeroflux):
    outflux = np.copy(zeroflux)
    outflux -= np.dot(scale, pcomp)
    return outflux

def chi2_gtf(obs, expect, err, dof):
    """
    calculates a chi squared of the model fit to the data
    """

    chisqu = 0.
    obs = obs
    expect  = expect
    err = err
    for i in range(len(obs)):
        chisqu += ((obs[i] - expect[i]) / err[i]) ** 2
    chisqu = chisqu / float(dof)
    return chisqu

def rms(model, data):
    """
    calculates a root mean square of the model fit to the data
    """

    rms = math.sqrt(np.sum((model - data) ** 2) / len(model))
    return rms

def do_lst_iter(bvs, cad, flux, nsigma, niter, method, order):
    """
    performs an iterative fit of the basis vectors to the light curve
    after each fit outliers further than nsigma from the fit are removed and the fit recalculated.

    The sigma actually means a median absolute deviation from the median.
    """
    iiter = 1
    fluxnew = np.copy(flux)
    lcnew = np.copy(cad)
    bvsnew = np.copy(bvs)
    if method == 'lst_sq':
        t = do_lsq_nlin(bvsnew, fluxnew)
    elif method == 'simplex':
        t = do_lsq_fmin_pow(bvsnew, fluxnew, order)
    elif method == 'llsq':
        t = do_lsq_uhat(bvsnew, fluxnew)

    bvsum = np.dot(t.T, bvsnew).reshape(-1)
    while (iiter < niter):
        iiter += 1
        matchrange = 1.4826 * nsigma * MAD_model(np.subtract(fluxnew, bvsum))
        mask = np.asarray(abs(fluxnew - bvsum) < matchrange)
        mask = mask.flatten()
        fluxnew = fluxnew[mask]
        lcnew = lcnew[mask]
        try:
            bvsnew = np.copy(bvsnew2)
        except:
            pass
        bvsnew2 = newpcompsarray(bvsnew, mask)
        for i in range(np.shape(bvsnew)[0]):
            bvsnew2[i] = bvsnew[i][mask]
        if method == 'llsq':
            t = do_lsq_uhat(bvsnew2, fluxnew)
        elif method == 'lst_sq':
            t = do_lsq_nlin(bvsnew2, fluxnew)
        elif method == 'simplex':
            t = do_lsq_fmin_pow(bvsnew2, fluxnew, order)
        bvsum = np.dot(t.T, bvsnew2).reshape(-1)

    return t, mask

def newpcompsarray(pcomp, mask):
    pcompnew = np.zeros((np.shape(pcomp)[0], len(mask[mask])))
    return pcompnew


def MAD_model(xx, minSd=1E-16):
    """Median Absolute Deviation"""
    absdev = abs(xx)
    mad = np.median(absdev, 0)
    mad = np.maximum(mad, np.multiply(np.ones(mad.shape, np.float32), (minSd / 1.48)))
    mad = np.asarray(mad)
    return mad


def make_outfile(fitsfile, outfile, flux_new, bvsum, version):
    """
    creates a fits file identical to the input fits file save from
    containing two extra columns - CBVSAP_MODL and CBVSAP_FLUX which are the
    sum of basis vectors fit to the data and the resulting corrected flux
    after the basis vector fit has been subtracted
    """

    if version == 1:
        unit = 'e-/cadence'
        flux_new = flux_new * 1625.3514 #convert to e-/cadence
    elif version == 2:
        unit = 'e-/s'
    col1 = pyfits.Column(name='CBVSAP_MODL', format='E13.7   ', unit=unit,
                         array=bvsum)
    col2 = pyfits.Column(name='CBVSAP_FLUX', format='E13.7   ', unit=unit,
                         array=flux_new)
    cols = fitsfile[1].columns + col1 + col2
    fitsfile[1] = pyfits.BinTableHDU.from_columns(cols,
                                                  header=fitsfile[1].header)
    fitsfile.writeto(outfile)

def do_plot(date, flux_old, flux_new, bvsum, cad, good_data, cad_nans, version,
    maskdata, outfile, noninteractive):

    plt.figure(figsize=[15, 8])
    plt.clf()

    if version == 1:
        barytime0 = float(int(date[0] / 100) * 100.0)
        date_sub = date - barytime0
        xlab = r'BJD $-$ {}'.format(barytime0+2400000.)
    elif version == 2:
        barytime0 = float(int((date[0] + 54833.) / 100) * 100.0)
        date_sub = date + 54833. - barytime0
        xlab = r'BJD $-$ {}'.format(barytime0+2400000.)

    try:
        nrm1 = len(str(int(flux_old.max()))) - 1
    except:
        nrm1 = 0
    flux_old_sub = flux_old / 10 ** nrm1
    bvsum_sub = bvsum / 10 ** nrm1
    ylab1 = r'10$^%d$ e$^-$ s$^{-1}$' % nrm1

    try:
        nrm2 = len(str(int(flux_new.max()))) - 1
    except:
        nrm2 = 0
    flux_new_sub = flux_new / 10 ** nrm2
    ylab2 = r'10$^%d$ e$^-$ s$^{-1}$' % nrm2

    xmin = min(date_sub)
    xmax = max(date_sub)
    ymin1 = min(min(flux_old_sub), min(bvsum_sub))
    ymax1 = max(max(flux_old_sub), max(bvsum_sub))
    ymin2 = min(flux_new_sub)
    ymax2 = max(flux_new_sub)
    xr = xmax - xmin
    yr1 = ymax1 - ymin1
    yr2 = ymax2 - ymin2

    ax1 = plt.subplot(211)

    blocks = split_on_nans(good_data,cad_nans)
    for i in range(len(blocks)):
        if i == 0:
            block = [blocks[0], blocks[i]]
        else:
            block = [blocks[i - 1], blocks[i]]
        mask = np.logical_and(cad >= block[0], cad <= block[1])
        plot_x = date_sub[mask]
        plot_y = flux_old_sub[mask]
        if np.nan in plot_y:
            break
        plt.scatter(plot_x, plot_y, color='#363636', linestyle='-', linewidth=1.0, marker='.', s=5)
        plot_y = bvsum_sub[mask]
        plt.plot(plot_x, plot_y, color='#c0392b', linestyle='-', linewidth=2.0)
    date2 = np.insert(date_sub, [0], [date_sub[0]])
    date2 = np.append(date2, [date_sub[-1]])
    flux2 = np.insert(flux_old_sub,[0], [0.0])
    flux2 = np.append(flux2, [0.0])
    plt.fill(date2, flux2, color='#a8a7a7', linewidth=0.0, alpha=0.2, label='Data')
    if maskdata is not None:
        for m in maskdata:
            pos = np.where((barytime0 + 2400000 + date2 > m[0]) & (barytime0 + 2400000 + date2 <= m[1]))[0]
            plt.fill_between(date2[pos], flux2[pos].min(), flux2[pos].max(), color='#c0392b', linewidth=0.0, alpha=0.3, label='Masked')

    plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
    if ymin1 - yr1 * 0.01 <= 0.0:
        plt.ylim(1.0e-10, ymax1 + yr1 * 0.01)
    else:
        plt.ylim(ymin1 - yr1 * 0.01, ymax1 + yr1 * 0.01)
    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab1, {'color' : 'k'})
    plt.grid(ls='--', alpha=0.3)
    plt.legend()

    ax2 = plt.subplot(212, sharex=ax1)
    for i in range(len(blocks)):
        if i == 0:
            block = [blocks[0], blocks[i]]
        else:
            block = [blocks[i - 1], blocks[i]]
        mask = np.logical_and(cad >= block[0], cad <= block[1])
        plot_x = date_sub[mask]
        plot_y = flux_new_sub[mask]
        if np.nan in plot_y:
            break
        plt.scatter(plot_x, plot_y, color='#363636', linestyle='-', linewidth=1.0, marker='.', s=5)
        plot_y = bvsum_sub[mask]

    date2 = np.insert(date_sub, [0], [date_sub[0]])
    date2 = np.append(date2, [date_sub[-1]])
    flux2 = np.insert(flux_new_sub, [0], [0.0])
    flux2 = np.append(flux2, [0.0])
    plt.fill(date2, flux2, fc='#a8a7a7', alpha=0.2,  linewidth=0.0)
    if maskdata is not None:
        for m in maskdata:
            pos = np.where((barytime0 + 2400000 + date2 >m[0]) & (barytime0 + 2400000 + date2 <=m [1]))[0]
            plt.fill_between(date2[pos], flux2[pos].min(), flux2[pos].max(), color='#c0392b', linewidth=0.0, alpha=0.3)
    plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)

    if ymin2-yr2*0.01 <= 0.0:
        plt.ylim(1.0e-10, ymax2 + yr2 * 0.01)
    else:
        plt.ylim(ymin2 - yr2 * 0.01, ymax2 + yr2 * 0.01)

    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab2, {'color' : 'k'})
    plt.grid(ls='--',alpha=0.3)
    plt.subplots_adjust(0.1, 0.1, 0.94, 0.94, 0.0, 0.0)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

    # render plot
    plt.savefig(re.sub('.fits', '.png', outfile), bbox_inches='tight')
    if not noninteractive:
        plt.show()

def split_on_nans(good_data, cad):
    blocks = []
    time_of_nans = cad[~good_data]
    if good_data[0]:
        blocks.append(cad[0])
    for i in range(1, len(time_of_nans)):
        if time_of_nans[i] - time_of_nans[i - 1] > 1:
            blocks.append(time_of_nans[i])
        if good_data[-1]:
            blocks.append(cad[-1])
    return blocks

def kepcotrend(infile, bvfile, listbv, outfile=None, fitmethod='llsq',
               fitpower=1, iterate=False, sigma=None, maskfile='',
               scinterp='linear', plot=False, noninteractive=False,
               overwrite=False, verbose=False, logfile='kepcotrend.log'):
    """
    kepcotrend -- Remove systematic trends Kepler light curves using
    cotrending basis vectors. The cotrending basis vectors files can be found
    here: http://archive.stsci.edu/kepler/cbv.html

    Simple Aperture Photometry (SAP) data often contain systematic trends
    associated with the spacecraft, detector and environment rather than the
    target. See the the Kepler data release notes for descriptions of
    systematics and the cadences that they affect. Within the Kepler pipeline
    these contaminants are treated during Pre-search Data Conditioning (PDC)
    and cleaned data are provided in the light curve files archived at MAST
    within the column PDCSAP_FLUX. The Kepler pipeline attempts to remove
    systematics with a combination of data detrending and cotrending against
    engineering telemetry from the spacecraft such as detector temperatures.
    These processes are imperfect but tackled in the spirit of correcting as
    many targets as possible with enough accuracy for the mission to meet
    exoplanet detection specifications.

    The imperfections in the method are most apparent in variable stars, those
    stars that are of most interest for stellar astrophysics. The PDC
    correction can occasionally hamper data analysis or, at worst, destroy
    astrophysical signal from the target. While data filtering (``kepoutlier``,
    ``kepfilter``) and data detrending with analytical functions
    (``kepdetrend``) often provide some mitigation for data artifacts, these
    methods require assumptions and often result in lossy data. An alternative
    viable approach is to identify the photometric variability common to all of
    the stars neighboring the target and subtract those trends from the target.
    In principle, the correct choice, weighting and subtraction of these common
    trends will leave behind a corrected flux time series which better
    represents statistically the true signal from the target.

    While GOs, KASC members and archive users wait for the Kepler project to
    release quarters of data, they do not have access to all the light curve
    data neighboring their targets and so cannot take the ensemble approach
    themselves without help. To mitigate this problem the Kepler Science Office
    have made available ancillary data which describes the systematic trends
    present in the ensemble flux data for each CCD channel. These data are
    known as the Cotrending Basis Vectors (CBVs). More details on the method
    used to generate these basis vectors will be provided in the Kepler Data
    Processing Handbook soon, but until that time, a summary of the method is
    given here. To create the initial basis set, that is the flux time series'
    that are used to make the cotrending basis vectors:

    The time series photometry of each star on a specific detector channel is
    normalized by its own median flux. One (unity) is subtracted from each time
    series so that the median value of the light curve is zero.
    The time series is divided by the root-mean square of the photometry.
    The correlation between each time series on the CCD channel is calculated
    using the median and root-mean square normalized flux.
    The median absolute correlation is then calculated for each star.
    All stars on the channel are sorted into ascending order of correlation.
    The 50 percent most correlated stars are selected.
    The median normalized fluxes only (as opposed to the root-mean square
    normalized fluxes) are now used for the rest of the process Singular Value
    Decomposition is applied to the matrix of correlated sources to create
    orthonormal basis vectors from the U matrix, sorted by their singular
    values.

    The archived cotrending basis vectors are a reduced-rank representation of
    the full set of basis vectors and consist of the 16 leading columns.

    To correct a SAP light curve, :math:`Fsap`, for systematic features,
    ``kepcotrend`` employs the cotrending basis vectors :math:`CBVi`. The task
    finds the coefficients :math:`Ai` which minimize

    .. math::

        Fcbv = Fsap - \sum_{i} Ai \cdot CBV_i

    The corrected light curve, Fcbv, can be tailored to the needs of the user
    and their scientific objective. The user decides which combination of basis
    vectors best removes systematics from their specific Kepler SAP light
    curve. In principle the user can choose any combination of cotrending basis
    vectors to fit to the data. However, experience suggests that most choices
    will be to decide how many sequential basis vectors to include in the fit,
    starting with first vector. For example a user is much more likely to
    choose a vector combination 1, 2, 3, 4, 5, 6 etc over e.g. a combination 1,
    2, 5, 7, 8, 10, 12. The user should always include at least the first two
    basis vectors. The number of basis vectors used is directly related to the
    scientific aims of the user and the light curve being analyzed and
    experimental iteration towards a target-specific optimal basis set is
    recommended. Occasionally kepcotrend over-fits the data and removes real
    astrophysical signal. This is particularly prevalent if too many basis
    vectors are used. A good rule of thumb is to start with two basis vectors
    and increase the number until there is no improvement, or signals which are
    thought to be astrophysical start to become distorted.

    The user is given a choice of fitting algorithm to use. For most purposes
    the linear least squares method is both the fastest and the most accurate
    because it gives the exact solution to the least squares problem. However
    we have found a few situations where the best solution, scientifically,
    comes from using the simplex fitting algorithm which performs something
    other than a least squares fit. Performing a least absolute residuals fit
    (fitpower=1.0), for example, is more robust to outliers.

    There are instances when the fit performs sub-optimally due to the presence
    of certain events in the light curve. For this reason we have included two
    options which can be used individually or simultaneously to improve the fit
    - iterative fitting and data masking. Iterative fitting performs the fit
    and rejects data points that are greater than a specified distance from the
    optimal fit before re-fitting. The lower threshold for data clipping is
    provided by the user as the number of sigma from the best fit. The clipping
    threshold is more accurately defined as the number of Median Absolute
    Deviations (MADs) multiplied by 1.4826. The distribution of MAD will be
    identical to the distribution of standard deviation if the distribution is
    Gaussian. We use MAD because in highly non-Gaussian distributions MAD is
    more robust to outliers than standard deviation.

    The code will print out the coefficients fit to each basis vector, the
    root-mean square of the fit and the chi-squared value of the fit. The rms
    and the chi-squared value include only the data points included in the fit
    so if an iterative fit is performed these clipped values are not included
    in this calculation.

    Parameters
    ----------
    infile : str
        the input file in the FITS format obtained from MAST
    outfile : str
        the output will be a fits file in the same style as the input file but
        with two additional columns: CBVSAP_MODL and CBVSAP_FLUX. The first of
        these is the best fitting linear combination of basis vectors.
        The second is the new flux with the basis vector sum subtracted. This
        is the new flux value.
    bvfile : str
        the name of the FITS file containing the basis vectors
    listbv : list of integers
        the basis vectors to fit to the data
    fitmethod : str
        fit using either the 'llsq' or the 'simplex' method. 'llsq' is usually
        the correct one to use because as the basis vectors are orthogonal.
        Simplex gives you option of using a different merit function - ie. you
        can minimise the least absolute residual instead of the least squares
        which weights outliers less
    fitpower : float
        if using a simplex you can chose your own power in the metir function
        - i.e. the merit function minimises :math:`abs(Obs - Mod)^P`.
        :math:`P = 2` is least squares, :math:`P = 1` minimises least absolutes
    iterate : bool
        should the program fit the basis vectors to the light curve data then
        remove data points further than 'sigma' from the fit and then refit
    maskfile : str
        this is the name of a mask file which can be used to define regions of
        the flux time series to exclude from the fit. The easiest way to create
        this is by using ``keprange`` from the PyKE set of tools. You can also
        make this yourself with two BJDs on each line in the file specifying
        the beginning and ending date of the region to exclude.
    scinterp : str
        the basis vectors are only calculated for long cadence data, therefore if
        you want to use short cadence data you have to interpolate the basis
        vectors. There are several methods to do this, the best of these probably
        being nearest which picks the value of the nearest long cadence data
        point.
        The options available are:

        * linear
        * nearest
        * zero
        * slinear
        * quadratic
        * cubic
    plot : bool
        Plot the data and result?
    non-interactive : bool
        If True, prevents the matplotlib window to pop up.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepcotrend kplr005110407-2009350155506_llc.fits ~/cbv/kplr2009350155506-q03-d25_lcbv.fits
        '1 2 3' --plot --verbose
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPCOTREND -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' bvfile={}'.format(bvfile)
            + ' listbv={} '.format(listbv)
            + ' fitmethod={}'.format(fitmethod)
            + ' fitpower={}'.format(fitpower)
            + ' iterate={}'.format(iterate)
            + ' sigma_clip={}'.format(sigma)
            + ' mask_file={}'.format(maskfile)
            + ' scinterp={}'.format(scinterp)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPCOTREND started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPCOTREND: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    if not kepio.fileexists(bvfile):
        message = 'ERROR -- KEPCOTREND: ' + bvfile + ' does not exist.'
        kepmsg.err(logfile, message, verbose)
    #lsq_sq - nonlinear least squares fitting and simplex_abs have been
    #removed from the options in PyRAF but they are still in the code!
    if fitmethod not in ['llsq','matrix','lst_sq','simplex_abs','simplex']:
        errmsg = 'Fit method must either: llsq, matrix, lst_sq or simplex'
        kepmsg.err(logfile, errmsg, verbose)

    if not is_numlike(fitpower) and fitpower is not None:
        errmsg = 'Fit power must be an real number or None'
        kepmsg.err(logfile, errmsg, verbose)

    if fitpower is None:
        fitpower = 1.

    # input data
    short = False
    try:
        test = str(instr[0].header['FILEVER'])
        version = 2
    except KeyError:
        version = 1

    table = instr[1].data
    if version == 1:
        if str(instr[1].header['DATATYPE']) == 'long cadence':
            quarter = str(instr[1].header['QUARTER'])
            module = str(instr[1].header['MODULE'])
            output = str(instr[1].header['OUTPUT'])
            channel = str(instr[1].header['CHANNEL'])
            lc_cad_o = table.field('cadence_number')
            lc_date_o = table.field('barytime')
            lc_flux_o = table.field('ap_raw_flux') / 1625.3468 #convert to e-/s
            lc_err_o = table.field('ap_raw_err') / 1625.3468 #convert to e-/s
        elif str(instr[1].header['DATATYPE']) == 'short cadence':
            short = True
            quarter = str(instr[1].header['QUARTER'])
            module = str(instr[1].header['MODULE'])
            output = str(instr[1].header['OUTPUT'])
            channel = str(instr[1].header['CHANNEL'])
            lc_cad_o = table.field('cadence_number')
            lc_date_o = table.field('barytime')
            lc_flux_o = table.field('ap_raw_flux') / 54.178 #convert to e-/s
            lc_err_o = table.field('ap_raw_err') / 54.178 #convert to e-/s

    elif version >= 2:
        if str(instr[0].header['OBSMODE']) == 'long cadence':
            quarter = str(instr[0].header['QUARTER'])
            module = str(instr[0].header['MODULE'])
            output = str(instr[0].header['OUTPUT'])
            channel = str(instr[0].header['CHANNEL'])
            lc_cad_o = table.field('CADENCENO')
            lc_date_o = table.field('TIME')
            lc_flux_o = table.field('SAP_FLUX')
            lc_err_o = table.field('SAP_FLUX_ERR')
        elif str(instr[0].header['OBSMODE']) == 'short cadence':
            short = True
            quarter = str(instr[0].header['QUARTER'])
            module = str(instr[0].header['MODULE'])
            output = str(instr[0].header['OUTPUT'])
            channel = str(instr[0].header['CHANNEL'])
            lc_cad_o = table.field('CADENCENO')
            lc_date_o = table.field('TIME')
            lc_flux_o = table.field('SAP_FLUX')
            lc_err_o = table.field('SAP_FLUX_ERR')

    if str(quarter) == str(4) and version == 1:
        lc_cad_o = lc_cad_o[lc_cad_o >= 11914]
        lc_date_o = lc_date_o[lc_cad_o >= 11914]
        lc_flux_o = lc_flux_o[lc_cad_o >= 11914]
        lc_err_o = lc_err_o[lc_cad_o >= 11914]

    if short and scinterp == None:
        errmsg = ('You cannot select None as the interpolation method '
                  'because you are using short cadence data and '
                  'therefore must use some form of interpolation. I '
                  'reccommend nearest if you are unsure.')
        kepmsg.err(logfile, errmsg, verbose)

    bvfiledata = pyfits.open(bvfile)
    bvdata = bvfiledata['MODOUT_{0}_{1}'.format(module, output)].data

    if int(bvfiledata[0].header['QUARTER']) != int(quarter):
        errmsg = ('CBV file and light curve file are from different '
                  'quarters. CBV file is from Q{0} and light curve is '
                  'from Q{1}'.format(int(bvfiledata[0].header['QUARTER']),
                                     int(quarter)))
        kepmsg.err(logfile, errmsg, verbose)

    if int(quarter) == 4 and int(module) == 3:
        errmsg = ('Approximately twenty days into Q4 Module 3 failed. '
                  'As a result, Q4 light curves contain these 20 day '
                  'of data. However, we do not calculate CBVs for '
                  'this section of data.')
        kepmsg.err(logfile, errmsg, verbose)

    #cut out infinites and zero flux columns
    lc_cad, lc_date, lc_flux, lc_err, good_data = cut_bad_data(lc_cad_o,
                                      lc_date_o, lc_flux_o, lc_err_o)
    #get a list of basis vectors to use from the list given
    #accept different seperators
    if len(listbv) == 1:
        bvlist = [listbv]
    else:
        listbv = listbv.strip()
        if listbv[1] in [' ', ',', ':', ';', '|', ', ']:
            separator = str(listbv)[1]
        else:
            message = ('You must separate your basis vector numbers to use '
                       'with \' \' \',\' \':\' \';\' or \'|\' and the '
                       'first basis vector to use must be between 1 and 9')
            kepmsg.err(logfile, message, verbose)

        bvlist = np.fromstring(listbv, dtype=int, sep=separator)

    if bvlist[0] == 0:
        errmsg = 'Must use at least one basis vector'
        kepmsg.err(logfile, errmsg, verbose)
    if short:
        bvdata.field('CADENCENO')[:] = ((((bvdata.field('CADENCENO')[:] +
                                        (7.5 / 15.) ) * 30.) - 11540.).round())
    bvectors = get_pcomp_list_newformat(bvdata, bvlist, lc_cad, short, scinterp)
    medflux = np.median(lc_flux)
    n_flux = (lc_flux / medflux) - 1
    n_err = np.sqrt(lc_err * lc_err / (medflux * medflux))

    if maskfile != '':
        domasking = True
        if not kepio.fileexists(maskfile):
            errmsg = 'Maskfile {} does not exist'.format(maskfile)
            kepmsg.err(logfile, errmsg, verbose)
    else:
        domasking = False

    if domasking:
        lc_date_masked = np.copy(lc_date)
        n_flux_masked = np.copy(n_flux)
        lc_cad_masked = np.copy(lc_cad)
        n_err_masked = np.copy(n_err)
        maskdata = np.atleast_2d(np.genfromtxt(maskfile, delimiter=','))
        mask = np.ones(len(lc_date_masked), dtype=bool)
        for maskrange in maskdata:
            if version == 1:
                start = maskrange[0] - 2400000.0
                end = maskrange[1] - 2400000.0
            elif version == 2:
                start = maskrange[0] - 2454833.
                end = maskrange[1] - 2454833.
            masknew = np.logical_xor(lc_date < start, lc_date > end)
            mask = np.logical_and(mask,masknew)
        lc_date_masked = lc_date_masked[mask]
        n_flux_masked = n_flux_masked[mask]
        lc_cad_masked = lc_cad_masked[mask]
        n_err_masked = n_err_masked[mask]
    else:
        lc_date_masked = np.copy(lc_date)
        n_flux_masked = np.copy(n_flux)
        lc_cad_masked = np.copy(lc_cad)
        n_err_masked = np.copy(n_err)

    bvectors_masked = get_pcomp_list_newformat(bvdata, bvlist, lc_cad_masked,
                                               short, scinterp)

    if iterate and sigma is None:
        errmsg = 'If fitting iteratively you must specify a clipping range'
        kepmsg.err(logfile, errmsg, verbose)

    #uses Pvals = yhat * U_transpose
    if iterate:
        coeffs, fittedmask = do_lst_iter(bvectors_masked, lc_cad_masked,
                                         n_flux_masked, sigma, 50., fitmethod,
                                         fitpower)
    else:
        if fitmethod == 'lst_sq':
            coeffs = do_lsq_nlin(bvectors_masked, n_flux_masked)
        elif fitmethod == 'simplex':
            coeffs = do_lsq_fmin_pow(bvectors_masked, n_flux_masked, fitpower)
        else:
            coeffs = do_lsq_uhat(bvectors_masked, n_flux_masked)

    coeffs = np.asarray(coeffs)
    flux_after = medflux * (n_flux + np.dot(coeffs.T, bvectors) + 1).reshape(-1)
    flux_after_masked = medflux * (n_flux_masked + np.dot(coeffs.T, bvectors_masked) + 1).reshape(-1)
    bvsum = np.dot(coeffs.T, bvectors).reshape(-1)
    bvsum_masked = np.dot(coeffs.T, bvectors_masked).reshape(-1)
    bvsum_nans = put_in_nans(good_data, bvsum)
    flux_after_nans = put_in_nans(good_data, flux_after)

    if plot:
        if not domasking:
            maskdata = None
        newmedflux = np.median(flux_after + 1)
        bvsum_un_norm = newmedflux * (1 - bvsum)
        do_plot(lc_date, lc_flux, flux_after, bvsum_un_norm, lc_cad,
                good_data, lc_cad_o, version, maskdata, outfile, noninteractive)

    print("Writing output file {}...".format(outfile))
    make_outfile(instr, outfile, flux_after_nans, bvsum_nans, version)
    # close input file
    instr.close()
    #print some results to screen:
    print('      -----      ')
    if iterate:
        flux_fit = n_flux_masked[fittedmask]
        sum_fit = bvsum_masked[fittedmask]
        err_fit = n_err_masked[fittedmask]
    else:
        flux_fit = n_flux_masked
        sum_fit = bvsum_masked
        err_fit = n_err_masked

    print('reduced chi2: {}'.format(chi2_gtf(flux_fit, sum_fit, err_fit,
                                             len(flux_fit) - len(coeffs))))
    print('rms: {}'.format(medflux * rms(flux_fit, sum_fit)))

    for i in range(len(coeffs)):
        print('Coefficient of CBV #{0}: {1}'.format(i + 1, coeffs[i]))
    print('      -----      ')

    # end time
    kepmsg.clock('KEPCOTREND completed at', logfile, verbose)

def kepcotrend_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Remove systematic trends in photometry using'
                          ' cotrending basis vectors (CBV)'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('cbvfile', help='Name of file containing the CBVs',
                        type=str)
    parser.add_argument('listbv', help='The CBVs to use', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepcotrend.'),
                        default=None)
    parser.add_argument('--method', '-m', help='Fitting method',
                        default='llsq', dest='fitmethod', type=str,
                        choices=['llsq', 'simplex', 'lst_sq'])
    parser.add_argument('--fitpower', '-f',
                        help='The index of the merit function (simplex only)',
                        default=1, type=float)
    parser.add_argument('--iterate', action='store_true',
                        help='Fit iteratively ', dest='iterate')
    parser.add_argument('--sigmaclip', type=float,
                        help='Sigma clip value when iteratively fitting',
                        default=None, dest='sigma')
    parser.add_argument('--maskfile', '-q',
                        help='Name of file containing a mask', default='',
                        dest='maskfile', type=str)
    parser.add_argument('--scinterp', type=str,
                        help='Short cadence interpolation method',
                        default='linear',
                        choices=['linear', 'nearest', 'slinear', 'quadratic',
                                 'cubic'])
    parser.add_argument('--plot', '-p', action='store_true',
                        help='Plot result?')
    parser.add_argument('--non-interactive', action='store_true',
                        help='Pop up matplotlib plot window?',
                        dest='noninteractive')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepcotrend.log', type=str)
    args = parser.parse_args()
    kepcotrend(args.infile, args.cbvfile, args.listbv, args.outfile,
               args.fitmethod, args.fitpower, args.iterate, args.sigma,
               args.maskfile, args.scinterp, args.plot, args.noninteractive,
               args.overwrite, args.verbose, args.logfile)
