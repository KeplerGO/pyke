from .utils import PyKEArgumentHelpFormatter
import re
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from . import kepmsg, kepio, kepkey, kepplot, kepfit, kepfunc


__all__ = ['kepsff']


def kepsff(infile, outfile=None, datacol='DETSAP_FLUX', cenmethod='moments',
           stepsize=5., npoly_cxcy=1, sigma_cxcy=10.0, npoly_ardx=6,
           npoly_dsdt=2, sigma_dsdt=3.0, npoly_arfl=3, sigma_arfl=3.0,
           plot=False, overwrite=False, verbose=False, logfile='kepsff.log'):
    """
    kepsff -- remove motion-correlated noise from aperture light curve data

    The systematic noise within aperture-constructed K2 photometry is dominated
    by the effects of boresight roll under the force of solar pressure. Roll is
    minimized by orienting the spacecraft in pointing directions that minimize
    spacecraft structure asymmetry normal to the Sun-vector and compensated for
    over campaign durations by regular thruster firings, typically every 6
    hours. kepsff provides an open source python implementation of the
    photometric correction for K2 motion systematics provided by Vanderburg and
    Johnson (2014).The method will also work upon archived Kepler data and
    provides a relatively-simple and CPU-friendly alternative to cotrending
    against basis-vector derived ensemble photometry (``kepcotrend``),
    pixel-level principal component analysis (``keppca``), and PSF fitting
    (``kepprfphot``). As well, as computational speed, this method has the
    advantage of requiring no external data in order to perform corrections.
    All required data is extracted from individual Kepler light curves stored
    at the archive, or K2 light curves constructed from archived target pixel
    files using the kepextract task. In the example figure above, the 12th
    magnitude target, EPIC 202093219, has a median 6.5-hour standard deviation
    of 60 parts-per-million (according to the tool kepstddev. After motion
    correction using ``kepsff``, this quantity is reduced to 36 parts-per-million.

    The name ``kepsff`` is derived from "Self-Flat-Fielding" (SFF) and is
    propagated from the Vanderburg and Johnson (2014) paper. However the name
    is somewhat misleading because the effects corrected for in the K2 case are
    dominated by aperture losses and source crowding rather than flat-field
    variations. In essence, ``kepsff`` corrects motion-induced
    aperture-photometry artifacts after constructing a position-flux relation
    from a detrended (astrophysics-removed) version of the light curve. This
    tool will find most-employment upon K2 data, especially in the early
    mission phases when the data archives will contain calibrated target pixel
    files but no light curves. However, in line with the philosophy of the PyKE
    tools, ``kepsff`` provides the user with a degree of customization to meet
    the demands of a broad range of target types and science goals. The
    long-term goal is to provide Kepler and K2 users with a range of reduction
    options and this tool is part of that software kit. Also, in line with the
    PyKE philosophy, these algorithms, as coded, are not perfect. We recommend
    and welcome tinkering and further development by the community in order to
    make these tools better. This case is also precisely what the K2 community
    requires in the sense that the archive users, here Andrew Vanderburg and
    John Johnson, are seeking solutions to their own data reduction needs,
    becoming self-sufficient, and reporting their findings and methods to us
    all.

    Preparation work is required before kepsff can provide a good correction.
    First, a light curve needs to be extracted a from Target Pixel File
    (``kepextract``) using a user-customized pixel aperture defined within
    ``kepmask``. The resulting light curve also needs to be detrended, i.e. the
    astrophysics needs to be removed as much as possible (kepflatten) in order
    for the correction to be as precise as possible. Note that the astrophysics
    is not removed permanently, but care has to be taken to prepare the data
    this way before confidence in the correction can be built for individual
    light curves on a case-by-case basis.

    Parameters
    ----------
    nfile : str
        The name of an input light curve FITS file. This file needs to contain
        a minimum set of information: Mid-exposure times, cadence numbers,
        simple-aperture flux, moment-derived and/or PSF-derived target
        centroids, quality flags and Kepler standard keywords. All of these
        properties are provided in archived light curves from the Kepler
        mission, but need to be calculated for the K2 mission. The PyKE tool
        kepextract calculates all these properties and provides an input file
        suitable for ``kepsff``. In addition, the input file also requires a
        column of detrended data, where long-term astrophysics and
        orbit-related systematics have been removed. The PyKE tool kepflatten
        performs this pre-step, producing a FITS column of detrended data (at
        low-frequencies) called DETSAP_FLUX.
    outfile : str
        The name of the output light curve FITS file. This file has the same
        structure as the input file with three content changes. The aperture
        photometry, DETSAP_FLUX, and detrended photometry, DETSAP_FLUX are
        corrected for motion-correlated noise. The quality flag column,
        SAP_QUALITY, includes a new bit-flag of value 131072 that identifies
        exposures collected during on-board thruster firings, typically
        occurring on a 6-hr cadence.
    datacol : str
        The name of the FITS data column containing detrended flux data (i.e.
        low-frequency structure removed). If the input file is a product of
        kepflatten then this column is named DETSAP_FLUX.
    cenmethod : str
        kepsff requires target position data on the CCD detector in order to
        correlate spacecraft boresight motion with systematic signal within the
        flux time-series. The typical Kepler/K2 light curve files have
        place-holders for two different measures of source centroid:

        * ``moments`` -- center-of-light within the pixel aperture
          (alternatively known as the zeroth-moment of the flux distribution
          over the pixels). These two arrays, the X and Y centroid over time,
          are stored in FITS columns MOM_CENTR1 and MOM_CENTR2.
        * ``psf`` -- the PSF-fitting method, the X and Y locations of best-fit
          PSF models to the aperture pixels. These are rarely available within
          archived Kepler data but can be stored in PSF_CENTR1 and PSF_CENTR2.
          Both sets of centroids are provided in the output from kepextract,
          where the moments data is generally preferred over the PSF data which
          derives from the simplifying assumptions that sources are not
          confused and are well-characterized by symmetric Gaussian profiles.
    stepsize : float [days]
        Ultimately, kepsff will construct a time-series of flux-correction
        factors based upon a correlation curve between spacecraft boresight
        motion and aperture-flux variation. Perhaps due to differential
        velocity aberration, this calibration is not stable over time. The
        somewhat clunky solution to this in the current version of the tool is
        re-calibrate the correlation and apply the subsequent correction over a
        series of discrete time windows. stepsize identifies the length of
        these windows in units of days. Some trial and error is required to
        optimize light curve correction. Window sizes of 4-6 days often provide
        the optimal correction, but astrophysical structure on these timescales
        can lead the user to adopt longer or shorter step sizes in order to
        reduce the unwanted impact of astrophysics on the correction.
    npoly_cxcy : int
        There is a sequence of four polynomial fits required to characterize
        and filter the input data in order to produce a flux-motion correlation
        curve. The following sequence of fit parameters is therefore quite
        dis-orienting until some familiarity is gained with hands-on
        experience. The diagnostic plots, provided when plot=True provide
        significant help here with the four fits delivered in the four
        left-most panels (figure 1). The first fit is to the centroid
        coordinates selected by cenmethod. These data points are provided in
        the top-left panel of the figure and generally trace very-close to a
        linear path across the detector. Fine-point motion is dominated by roll
        around the boresight and therefore small-angle, near-linear behavior is
        unsurprising. The two black lines represent the eigenvectors of the
        covariance matrix of the centroids. These eigenvectors are used to
        rotate the centroid distribution so that the direction of motion (x')
        is aligned to the x-axis of the rotated system. The rotation
        functionality requires user-optimization in one aspect. Obvious
        outliers need to be removed from the sample before the eigenvectors are
        calculated. We achieve this iterative sigma-clipping. The user defines
        the order of a polynomial fit (1st order is recommended) and a
        threshold distance from the best fit where data points outside that
        threshold are removed before the remaining data are re-fit. Iterations
        continue until no additional data points are rejected. The rejection
        threshold is provided by the argument sigma_cxcy. Points lying below
        the rejection threshold are plotted as green symbols whereas data
        points rejected by this process will be plotted red. In the example
        above, no data points are rejected. The most likely cause of outliers
        is data obtained while the spacecraft is not guiding on reference
        stars. We predict such occasions to be rare during optimal operations,
        but such observations were common during the engineering tests,
        including the first half of campaign 0.
    sigma_cxcy : float [sigma]
        Threshold for rejecting centroid data. The threshold is the number of
        1-sigma standard deviation distances from the best-fit. Points further
        from the best-fit than the threshold provided will be considered
        unrepresentative of the flux-motion correlation that the tool is
        attempting to characterize. These points are not thrown away, they are
        corrected in the final calculation with the rest of the data, but they
        play no further part in the calibration of the correlation.
    npoly_ardx : int
        Step 2 (top-middle panel) is the determination of the relatively small
        deviation of centroid motion from linear. Here, in the rotated frame,
        we fit a high-order polynomial to the centroid data filtered and
        rotated in the first step. The green curve is the integrated path
        length along that curve (s) minus the linear term. The calculation is
        provided by equation 4 of Vanderburg and Johnson (2014). The red curve
        is a polynomial fit to the integral which subsequently allows us to
        cheaply calculate the absolute motion of the centroids relative to
        linear location along the dominant eigenvector in the unrotated frame.
        Producing a 'fit of a fit' is somewhat inelegant and open source
        developers could be challenged here to add some additional
        mathematical rigor to this step. However this section of the
        calculation has, so far, not proved to be critical - the path length
        integral is very close to linear, so cleaning this calculation up is
        low-priority currently.
    npoly_dsdt : int
        The third step (lower-left panel) is a prescription to filter out
        exposures collected during thruster firings. Firing occur typically on
        a 6-hr cadence and occur when the spacecraft has rolled far enough to
        trigger a pointing adjustment. These data cannot be corrected using the
        current method. They are flagged here so that they do not bias the
        motion-flux correlation calculation and so that they can be documented
        within the quality flags of the output file. A low-order polynomial is
        fit to the time derivative along the path length (ds/dt) using 3-sigma
        iterative clipping. ``npoly_dsdt`` provides the user with flexibility
        to choose the polynomial order. In the plot, the best fit is shown as a
        red, solid line.
    sigma_dsdt : float [sigma]
        This is a threshold limit in units of 1-sigma standard deviation from
        the best-fit. Data points falling outside of the threshold are
        more-likely-than-not collected during thruster firings. These points
        are flagged in the SAP_QUALITY column of the output file with the bit
        131072 and is not employed to calibrate the flux-motion relation. The
        threshold curve is plot as a red, dashed line and rejected points are
        colored red.
    npoly_arfl : int
        The fourth panel (lower-middle) provides the subsequent target motion
        vs aperture flux distribution. The red line is polynomial fit, the
        order of which is chosen by the user here. Iterative
        :math:`\sigma`-clipping is again required - the many outliers below
        the best-fit are astrophysical in nature (binary star eclipses) and
        would bias the fit if they were unrecognized or simply allowed to.
        The user has flexibility to change the order and clipping threshold
        in order to provide the most robust fit and there is a degree of
        objectivity in this process. A low polynomial is generally adequate,
        this example employs a 3rd order function.
    sigma_arfl : float [sigma]
        This is the threshold limit in units of 1-sigma standard deviation from
        the best-fit polynomial to the motion-flux relation. Data points
        falling outside of the threshold are in effect rejected from the
        calibration. Such points potentially contain astrophysical signal that
        would both bias the calculation and damp the astrophysical signal
        within the data after the correction has been made. The best-fit red
        line provides the correction factors within the data window. For each
        point along the aperture-derived flux curve (upper-right), the position
        along the arc can be calculated and the multiplicative correction
        factor to the flux can determined directly from the best-fit
        motion-flux relation. The subsequent corrected light curve is provided
        lower-right on the same plotting scale as above. The red points are
        those flagged as potential thruster firing and their 6-hr cadence
        suggests these events have generally been flagged well.
    plot : bool
        If true, diagnostic plots identical to figure 1 above will be rendered
        and also saved as PNG files. There will be a different plot for every
        time window defined by the stepsize parameters. If the output FITS file
        is named filename.fits then each PNG file will be named
        filename_nn.png, where nn is a sequential number beginning with 1.
    overwrite : bool
        Overwrite the output FITS file? if **overwrite** is **False** and an
        existing file has the same name as outfile then the task will stop with
        an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = ('KEPSFF -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' cenmethod={}'.format(cenmethod)
            + ' stepsize={}'.format(stepsize)
            + ' npoly_cxcy={}'.format(npoly_cxcy)
            + ' sigma_cxcy={}'.format(sigma_cxcy)
            + ' npoly_ardx={}'.format(npoly_ardx)
            + ' npoly_dsdt={}'.format(npoly_dsdt)
            + ' sigma_dsdt={}'.format(sigma_dsdt)
            + ' npoly_arfl={}'.format(npoly_arfl)
            + ' sigma_arfl={}'.format(sigma_arfl)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPSFF started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPSFF: {} exists. Use overwrite=True'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)
    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence
    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # determine sequence of windows in time
    frametim = instr[1].header['FRAMETIM']
    num_frm = instr[1].header['NUM_FRM']
    exptime = frametim * num_frm / 86400
    tstart = table.field('TIME')[0]
    tstop = table.field('TIME')[-1]
    winedge = np.arange(tstart, tstop, stepsize)
    if tstop > winedge[-1] + stepsize / 2:
        winedge = np.append(winedge, tstop)
    else:
        winedge[-1] = tstop
    winedge = (winedge - tstart) / exptime
    winedge = winedge.astype(int)
    if len(table.field('TIME')) > winedge[-1] + 1:
        winedge = np.append(winedge, len(table.field('TIME')))
    elif len(table.field('TIME')) < winedge[-1]:
        winedge[-1] = len(table.field('TIME'))
    # step through the time windows
    for iw in tqdm(range(1, len(winedge))):
        t1 = winedge[iw - 1]
        t2 = winedge[iw]
        # filter input data table
        work1 = np.array([table.field('TIME')[t1:t2],
                          table.field('CADENCENO')[t1:t2],
                          table.field(datacol)[t1:t2],
                          table.field('MOM_CENTR1')[t1:t2],
                          table.field('MOM_CENTR2')[t1:t2],
                          table.field('PSF_CENTR1')[t1:t2],
                          table.field('PSF_CENTR2')[t1:t2],
                          table.field('SAP_QUALITY')[t1:t2]],'float64')
        work1 = np.rot90(work1, 3)
        work2 = work1[(work1[:, 0] == 0.0) | (work1[:, 0] > 1e5)]


        # assign table columns
        intime = work2[:, 7] + bjdref
        cadenceno = work2[:, 6].astype(int)
        indata = work2[:, 5]
        mom_centr1 = work2[:, 4]
        mom_centr2 = work2[:, 3]
        psf_centr1 = work2[:, 2]
        psf_centr2 = work2[:, 1]
        sap_quality = work2[:, 0]
        if cenmethod == 'moments':
            centr1 = np.copy(mom_centr1)
            centr2 = np.copy(mom_centr2)
        else:
            centr1 = np.copy(psf_centr1)
            centr2 = np.copy(psf_centr2)

        # fit centroid data with low-order polynomial
        cfit = np.zeros((len(centr2)))
        csig = np.zeros((len(centr2)))
        functype = getattr(kepfunc, 'poly' + str(npoly_cxcy))
        pinit = np.array([np.nanmean(centr2)])
        if npoly_cxcy > 0:
            for j in range(npoly_cxcy):
                pinit = np.append(pinit,0.0)
        try:
            coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit, centr1, centr2, None,
                               sigma_cxcy, sigma_cxcy, 10, logfile, verbose)
            for j in range(len(coeffs)):
                cfit += coeffs[j] * np.power(centr1, j)
                csig[:] = sigma
        except:
            warnmsg = ('WARNING -- KEPSFF: could not fit centroid data with'
                       ' polynomial. There are no data points within the'
                       ' range of input rows {} - {}. Either increase the'
                       ' stepsize (with an appreciation of the effects on'
                       ' light curve quality this will have!), or better yet'
                       ' - cut the timeseries up to remove large gaps in the'
                       ' input light curve using kepclip.'.format(t1, t2))
            kepmsg.warn(logfile, warnmsg, verbose)
            continue

        # reject outliers
        time_good = np.array([], 'float64')
        centr1_good = np.array([], 'float32')
        centr2_good = np.array([], 'float32')
        flux_good = np.array([], 'float32')
        cad_good = np.array([], 'int')
        for i in range(len(cfit)):
            if abs(centr2[i] - cfit[i]) < sigma_cxcy * csig[i]:
                time_good = np.append(time_good, intime[i])
                centr1_good = np.append(centr1_good, centr1[i])
                centr2_good = np.append(centr2_good, centr2[i])
                flux_good = np.append(flux_good, indata[i])
                cad_good = np.append(cad_good, cadenceno[i])

        # covariance matrix for centroid time series
        centr = np.concatenate([[centr1_good] - np.nanmean(centr1_good),
                                [centr2_good] - np.nanmean(centr2_good)])
        covar = np.cov(centr)
        # eigenvector eigenvalues of covariance matrix
        [_, evec] = np.linalg.eigh(covar)
        ex = np.arange(-10.0, 10.0, 0.1)
        epar = evec[1, 1] / evec[0, 1] * ex
        enor = evec[1, 0] / evec[0, 0] * ex
        ex = ex + np.nanmean(centr1)
        epar = epar + np.nanmean(centr2_good)
        enor = enor + np.nanmean(centr2_good)
        # rotate centroid data
        centr_rot = np.dot(evec.T, centr)
        # fit polynomial to rotated centroids
        rfit = np.zeros((len(centr2)))
        rsig = np.zeros((len(centr2)))
        functype = getattr(kepfunc, 'poly' + str(npoly_ardx))
        pinit = np.array([np.nanmean(centr_rot[0,:])])
        pinit = np.array([1.0])
        if npoly_ardx > 0:
            for j in range(npoly_ardx):
                pinit = np.append(pinit,0.0)
        try:
            coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit, centr_rot[1, :],
                               centr_rot[0, :], None, 100.0, 100.0, 1, logfile,
                               verbose)
        except:
            warnmsg = ('WARNING -- KEPSFF: could not fit rotated centroid data'
                      ' with polynomial')
            kepmsg.warn(logfile, warnmsg, verbose)
            continue
        rx = np.linspace(np.nanmin(centr_rot[1, :]),
                         np.nanmax(centr_rot[1, :]), 100)
        ry = np.zeros((len(rx)))
        for i in range(len(coeffs)):
            ry = ry + coeffs[i] * np.power(rx,i)

        # calculate arclength of centroids
        s = np.zeros((len(rx)))
        for i in range(1, len(s)):
            work3 = ((ry[i] - ry[i - 1]) / (rx[i] - rx[i-1])) ** 2
            s[i] = s[i - 1] + np.sqrt(1.0 + work3) * (rx[i] - rx[i-1])

        # fit arclength as a function of strongest eigenvector
        sfit = np.zeros((len(centr2)))
        ssig = np.zeros((len(centr2)))
        functype = getattr(kepfunc, 'poly' + str(npoly_ardx))
        pinit = np.array([np.nanmean(s)])
        if npoly_ardx > 0:
            for j in range(npoly_ardx):
                pinit = np.append(pinit, 0.0)
        try:
            acoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit, rx, s, None, 100.0, 100.0, 100,
                               logfile, verbose)
        except:
            warnmsg = ('WARNING -- KEPSFF: could not fit arclength data'
                       ' with polynomial')
            kepmsg.warn(logfile, warnmsg, verbose)
            continue

        # correlate arclength with detrended flux
        t = np.copy(time_good)
        c = np.copy(cad_good)
        y = np.copy(flux_good)
        z = centr_rot[1, :]
        x = np.zeros((len(z)))
        for i in range(len(acoeffs)):
            x = x + acoeffs[i] * np.power(z, i)

        # calculate time derivative of arclength s
        dx = np.zeros((len(x)))
        for i in range(1, len(x)):
            dx[i] = (x[i] - x[i - 1]) / (t[i] - t[i - 1])
        dx[0] = dx[1]

        # fit polynomial to derivative and flag outliers (thruster firings)
        dfit = np.zeros((len(dx)))
        dsig = np.zeros((len(dx)))
        functype = getattr(kepfunc, 'poly' + str(npoly_dsdt))
        pinit = np.array([np.nanmean(dx)])
        if npoly_dsdt > 0:
            for j in range(npoly_dsdt):
                pinit = np.append(pinit, 0.0)
        try:
            dcoeffs, errors, covar, iiter, dsigma, chi2, dof, fit, dumx, dumy = \
                kepfit.lsqclip(functype, pinit, t, dx, None, 3.0, 3.0, 10,
                               logfile, verbose)
        except:
            warnmsg = ('WARNING -- KEPSFF: could not fit arclength derivative'
                       ' with polynomial.')
            kepmsg.warn(logfile, warnmsg, verbose)
            continue
        for i in range(len(dcoeffs)):
            dfit = dfit + dcoeffs[i] * np.power(t, i)
        centr1_pnt = np.array([], 'float32')
        centr2_pnt = np.array([], 'float32')
        time_pnt = np.array([], 'float64')
        flux_pnt = np.array([], 'float32')
        dx_pnt = np.array([], 'float32')
        s_pnt = np.array([], 'float32')
        time_thr = np.array([], 'float64')
        flux_thr = np.array([], 'float32')
        dx_thr = np.array([], 'float32')
        thr_cadence = []
        for i in range(len(t)):
            if (dx[i] < dfit[i] + sigma_dsdt * dsigma
                and dx[i] > dfit[i] - sigma_dsdt * dsigma):
                time_pnt = np.append(time_pnt, time_good[i])
                flux_pnt = np.append(flux_pnt, flux_good[i])
                dx_pnt = np.append(dx_pnt, dx[i])
                s_pnt = np.append(s_pnt, x[i])
                centr1_pnt = np.append(centr1_pnt, centr1_good[i])
                centr2_pnt = np.append(centr2_pnt, centr2_good[i])
            else:
                time_thr = np.append(time_thr, time_good[i])
                flux_thr = np.append(flux_thr, flux_good[i])
                dx_thr = np.append(dx_thr, dx[i])
                thr_cadence.append(cad_good[i])

        # fit arclength-flux correlation
        cfit = np.zeros((len(time_pnt)))
        csig = np.zeros((len(time_pnt)))
        functype = getattr(kepfunc, 'poly' + str(npoly_arfl))
        pinit = np.array([np.nanmean(flux_pnt)])
        if npoly_arfl > 0:
            for j in range(npoly_arfl):
                pinit = np.append(pinit, 0.0)
        try:
            ccoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plx, ply = \
                kepfit.lsqclip(functype, pinit, s_pnt, flux_pnt, None,
                               sigma_arfl, sigma_arfl, 100, logfile, verbose)
        except:
            warnmsg = ('WARNING -- KEPSFF: could not fit arclength-flux'
                       ' correlation with polynomial')
            kepmsg.warn(logfile, warnmsg, verbose)
            continue

        # correction factors for unfiltered data
        centr = np.concatenate([[centr1] - np.nanmean(centr1_good),
                                [centr2] - np.nanmean(centr2_good)])
        centr_rot = np.dot(evec.T, centr)
        yy = np.copy(indata)
        zz = centr_rot[1, :]
        xx = np.zeros((len(zz)))
        cfac = np.zeros((len(zz)))
        for i in range(len(acoeffs)):
            xx = xx + acoeffs[i] * np.power(zz, i)
        for i in range(len(ccoeffs)):
            cfac = cfac + ccoeffs[i] * np.power(xx, i)

        # apply correction to flux time-series
        out_detsap = indata / cfac
        # split time-series data for plotting
        tim_gd = np.array([], 'float32')
        flx_gd = np.array([], 'float32')
        tim_bd = np.array([], 'float32')
        flx_bd = np.array([], 'float32')
        for i in range(len(indata)):
            if intime[i] in time_pnt:
                tim_gd = np.append(tim_gd, intime[i])
                flx_gd = np.append(flx_gd, out_detsap[i])
            else:
                tim_bd = np.append(tim_bd, intime[i])
                flx_bd = np.append(flx_bd, out_detsap[i])
        if plot:
            # plot style and size
            #kepplot.define(16, 14, logfile, verbose)
            plt.figure(figsize=[20, 8])
            plt.clf()

            # plot x-centroid vs y-centroid
            ax = kepplot.location([0.04, 0.57, 0.16, 0.41])
            px = np.copy(centr1)
            py = np.copy(centr2)
            pxmin = px.min()
            pxmax = px.max()
            pymin = py.min()
            pymax = py.max()
            pxr = pxmax - pxmin
            pyr = pymax - pymin
            pad = 0.05
            if pxr > pyr:
                dely = (pxr - pyr) / 2
                plt.xlim(pxmin - pxr * pad, pxmax + pxr * pad)
                plt.ylim(pymin - dely - pyr * pad, pymax + dely + pyr * pad)
            else:
                delx = (pyr - pxr) / 2
                plt.ylim(pymin - pyr * pad, pymax + pyr * pad)
                plt.xlim(pxmin - delx - pxr * pad, pxmax + delx + pxr * pad)
            plt.plot(px, py, color='#980000', markersize=5, marker='D', ls='')
            plt.plot(centr1_good, centr2_good,color='#009900', markersize=5,
                     marker='D', ls='')
            plt.plot(ex,epar,color='k',ls='-')
            plt.plot(ex,enor,color='k',ls='-')
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            kepplot.labels('CCD Column', 'CCD Row', 'k', 16)
            plt.grid()

            # plot arclength fits vs drift along strongest eigenvector
            ax = kepplot.location([0.24, 0.57, 0.16, 0.41])
            px = rx - rx[0]
            py = s - rx - (s[0] - rx[0])
            py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
            kepplot.RangeOfPlot(px, py, 0.05, False)
            plt.plot(px, py, color='#009900', markersize=5, marker='D', ls='')
            px = plotx - rx[0]
            py = ploty - plotx - (s[0] - rx[0])
            py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
            plt.plot(px, py, color='r', ls='-', lw=3)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            ylab = re.sub(' e\S+', ' pixels)', ylab)
            ylab = re.sub(' s\S+', '', ylab)
            ylab = re.sub('Flux', 's $-$ x\'', ylab)
            kepplot.labels('Linear Drift [x\'] (pixels)', ylab, 'k', 16)
            plt.grid()

            # plot time derivative of arclength s
            ax = kepplot.location([0.04,0.08,0.16,0.41])
            px = np.copy(time_pnt)
            py = np.copy(dx_pnt)
            px, xlab = kepplot.cleanx(px,logfile,verbose)
            kepplot.RangeOfPlot(px, dx, 0.05, False)
            plt.plot(px, py, color='#009900', markersize=5, marker='D', ls='')
            try:
                px = np.copy(time_thr)
                py = np.copy(dx_thr)
                px, xlab = kepplot.cleanx(px, logfile, verbose)
                plt.plot(px, py, color='#980000', markersize=5, marker='D', ls='')
            except:
                pass
            px = np.copy(t)
            py = np.copy(dfit)
            px, xlab = kepplot.cleanx(px, logfile, verbose)
            plt.plot(px, py, color='r', ls='-', lw=3)
            py = np.copy(dfit + sigma_dsdt * dsigma)
            plt.plot(px, py, color='r', ls='--', lw=3)
            py = np.copy(dfit-sigma_dsdt*dsigma)
            plt.plot(px,py,color='r',ls='--',lw=3)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            kepplot.labels(xlab, 'ds/dt (pixels day$^{-1}$)', 'k', 16)
            plt.grid()
            # plot relation of arclength vs detrended flux
            ax = kepplot.location([0.24, 0.08, 0.16, 0.41])
            px = np.copy(s_pnt)
            py = np.copy(flux_pnt)
            py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
            kepplot.RangeOfPlot(px, py, 0.05, False)
            plt.plot(px, py, color='#009900', markersize=5, marker='D', ls='')
            plt.plot(plx, ply, color='r', ls='-', lw=3)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            kepplot.labels('Arclength [s] (pixels)', ylab, 'k', 16)
            plt.grid()

            # plot aperture photometry
            kepplot.location([0.44, 0.53, 0.55, 0.45])
            px, xlab = kepplot.cleanx(intime, logfile, verbose)
            py, ylab = kepplot.cleany(indata, 1.0, logfile, verbose)
            kepplot.RangeOfPlot(px, py, 0.01, True)
            kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2, True)
            kepplot.labels(' ',ylab,'k',16)
            plt.setp(plt.gca(),xticklabels=[])
            kepplot.labels(xlab,re.sub('Flux','Aperture Flux',ylab),'k',16)
            plt.grid()
            kepplot.location([0.44, 0.08, 0.55, 0.45])
            kepplot.RangeOfPlot(px, py, 0.01, True)
            px, xlab = kepplot.cleanx(tim_gd, logfile, verbose)
            py, ylab = kepplot.cleany(flx_gd, 1.0, logfile, verbose)
            kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2, True)
            try:
                px, xlab = kepplot.cleanx(tim_bd,logfile,verbose)
                py = np.copy(flx_bd)
                plt.plot(px, py,color='#980000',markersize=5,marker='D',ls='')
            except:
                pass
            kepplot.labels(xlab,re.sub('Flux', 'Corrected Flux', ylab), 'k', 16)
            plt.grid()
            # render plot
            plt.show()
            plt.savefig(re.sub('.fits','_%d.png' % (iw + 1), outfile))

        # correct fluxes within the output file
        intime = work1[:, 7] + bjdref
        cadenceno = work1[:, 6].astype(int)
        indata = work1[:, 5]
        mom_centr1 = work1[:, 4]
        mom_centr2 = work1[:, 3]
        psf_centr1 = work1[:, 2]
        psf_centr2 = work1[:, 1]
        centr1 = np.copy(mom_centr1)
        centr2 = np.copy(mom_centr2)
        centr = np.concatenate([[centr1] - np.nanmean(centr1_good),
                               [centr2] - np.nanmean(centr2_good)])
        centr_rot = np.dot(evec.T, centr)
        yy = np.copy(indata)
        zz = centr_rot[1,:]
        xx = np.zeros((len(zz)))
        cfac = np.zeros((len(zz)))
        for i in range(len(acoeffs)):
            xx = xx + acoeffs[i] * np.power(zz, i)
        for i in range(len(ccoeffs)):
            cfac = cfac + ccoeffs[i] * np.power(xx, i)
        out_detsap = yy / cfac
        instr[1].data.field('SAP_FLUX')[t1:t2] /= cfac
        instr[1].data.field('PDCSAP_FLUX')[t1:t2] /= cfac
        try:
            instr[1].data.field('DETSAP_FLUX')[t1:t2] /= cfac
        except:
            pass

        # add quality flag to output file for thruster firings
        for i in range(len(intime)):
            if cadenceno[i] in thr_cadence:
                instr[1].data.field('SAP_QUALITY')[t1 + i] += 131072
    # write output file
    kepmsg.log(logfile, "Writing output file {}...".format(outfile), True)
    instr.writeto(outfile)
    # close input file
    instr.close()
    # end time
    kepmsg.clock('KEPSFF completed at', logfile, verbose)

def kepsff_main():
    import argparse
    parser = argparse.ArgumentParser(
             description='Correct aperture photmetry using target motion',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input FITS file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepsff.'),
                        default=None)
    parser.add_argument('--datacol', default='DETSAP_FLUX',
                        help='Name of data column', type=str)
    parser.add_argument('--cenmethod', default='moments',
                        help='Centroid method', type=str,
                        choices=['moments','psf'])
    parser.add_argument('--stepsize', default=5.,
                        help='Stepsize over which to calibrate data [days]',
                        type=float)
    parser.add_argument('--npoly_cxcy', default=1,
                        help='Order of ploynomial fit to target centroids',
                        type=int)
    parser.add_argument('--sigma_cxcy', default=10.0,
                        help=('Sigma-clipping threshold for fit to target'
                              ' centroids [sigma]'),
                        type=float)
    parser.add_argument('--npoly_ardx', default=6,
                        help=('Order of ploynomial fit for thruster firing'
                              ' detection'), type=int)
    parser.add_argument('--npoly_dsdt', default=2,
                        help=('Order of ploynomial fit for thruster firing'
                              ' detection'), type=int)
    parser.add_argument('--sigma_dsdt', default=3.0,
                        help=('Sigma-clipping threshold for thruster firing'
                              ' detection [sigma]'), type=float)
    parser.add_argument('--npoly_arfl', default=3,
                        help=('Order of ploynomial for arclength-flux'
                              ' calibration'), type=int)
    parser.add_argument('--sigma_arfl', default=3.0,
                        help=('Sigma-clipping threshold for arclength-flux'
                              ' calibration [sigma]'),
                        type=float)
    parser.add_argument('--plot', action='store_true',
                        help='Save hardcopies of the plots?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepsff.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepsff(args.infile, args.outfile, args.datacol, args.cenmethod,
           args.stepsize, args.npoly_cxcy, args.sigma_cxcy, args.npoly_ardx,
           args.npoly_dsdt, args.sigma_dsdt, args.npoly_arfl,
           args.sigma_arfl, args.plot, args.overwrite, args.verbose,
           args.logfile)
