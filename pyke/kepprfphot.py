import math
import multiprocessing
import itertools
import glob
import sys
import time
import re
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline
from . import kepio, kepmsg, kepkey, kepplot, kepfit, kepfunc
from .utils import PyKEArgumentHelpFormatter


__all__ = ['kepprfphot']


def kepprfphot(infile, prfdir, columns, rows, fluxes, border=0,
               background=False, focus=False, ranges='0,0', xtol=1e-4,
               ftol=1e-2, qualflags=False, outfile=None, plot=False, overwrite=False,
               verbose=False, logfile='kepprfphot.log'):
    """
    kepprfphot -- Fit a PSF model to time series observations within a Target
    Pixel File

    Parameters
    ----------
    nfile : str
        The name of a MAST standard format FITS file containing Kepler Target
        Pixel data within the first data extension.
    columns : str or list
        A starting guess for the CCD column position(s) of the source(s) that
        are to be fit. The model is unlikely to converge if the guess is too
        far away from the correct location. A rule of thumb is to provide a
        guess within 1 CCD pixel of the true position. If more than one source
        is being modeled then the column positions of each are separated by a
        comma. The same number of sources in the columns, rows and fluxes field
        is a requirement of this task.
    rows : str or list
        A starting guess for the CCD row position(s) of the source(s) that are
        to be fit. The model is unlikely to converge if the guess is too far
        away from the correct location. A rule of thumb is to provide a guess
        within 1 CCD pixel of the true position. If more than one source is
        being modeled then the row positions of each are separated by a comma.
        The same number of sources in the columns, rows and fluxes field is a
        requirement of this task.
    fluxes : str or list
        A starting guess for the flux(es) of the source(s) that are to be fit.
        Fit convergence is not particularly reliant on the accuracy of these
        guesses, but the fit will converge faster the more accurate the guess.
        If more than one source is being modeled then the row positions of
        each are separated by a comma. The same number of sources in the
        columns, rows and fluxes field is a requirement of this task.
    prfdir : str
        The full or relative directory path to a folder containing the Kepler
        PSF calibration. Calibration files can be downloaded from the Kepler
        focal plane characteristics page at the MAST here:
        http://archive.stsci.edu/missions/kepler/fpc/prf/.
    border : int
        If a background is included in the fit then it is modeled as a
        two-dimensional polynomial. This parameter is the polynomial order.
        A zero-order polynomial is generally recommended.
    background : bool
        Whether to include a background component in the model. If ``True``
        the background will be represented by a two-dimensional polynomial of
        order border. This functionality is somewhat experimental, with one eye
        upon potential background gradients across large masks or on those
        detectors more prone to pattern noise. Generally it is recommended to
        set background as ``False``.
    focus : bool
        Whether to include pixel scale and focus rotation with the fit
        parameters of the model. This is also an experimental function. This
        approach does not attempt to deal with inter- or intra-pixel
        variations. The recommended use is currently to set focus as ``False``.
    ranges : str
        The user can choose specific time ranges of data on which to work. This
        could, for example, avoid removing known stellar flares from a dataset
        Time ranges are supplied as comma-separated pairs of Barycentric Julian
        Dates (BJDs). Multiple ranges are separated by a semi-colon.
        An example containing two time ranges is::

            '2455012.48517,2455014.50072;2455022.63487,2455025.08231'

        If the user wants to correct the entire time series then providing
        ranges = '0,0' will tell the task to operate on the whole time series.
    xtol : float
        The dimensionless, relative model parameter convergence criterion for
        the fit algorithm.
    ftol : float
        The dimensionless, relative model residual convergence criterion for
        the fit algorithm.
    qualflags : bool
        If qualflags is ``False``, archived observations flagged with any
        quality issue will not be fit.
    outfile : str
        kepprfphot creates two types of output file containing fit results and
        diagnostics. ``outfile.png`` contains a time series plot of fit
        parameters, residuals and chi-squared. ``outfile.fits`` contains a
        table of the same properties, consistent with Kepler archive light
        curve files. The FITS column PSF_FLUX contains the flux time-series in
        units of e-/s derived by integrating under the best-fit PRF model.
        PSF_BKG provides the best-fit background (if calculated) averaged over
        all mask pixels in units of e-/s/pixel. PSF_CENTR1 provides the
        best-fit PSF centroid position in the CCD column direction, in CCD
        pixel units. Similarly, PSF_CENTR2 provides the best-fit PSF centroid
        position in the CCD row direction, in CCD pixel units. If calculated,
        PSF_FOCUS1 and PSF_FOCUS2 provide scale factors in the column and row
        dimensions by which the CCD pixel scale is adjusted to approximate
        focus variation. PSF_ROTATION provides the angle by which the scaled
        PSF model was rotated on the focal plane in order to yield a best fit.
        The table column PSF_RESIDUAL provides the sum of all mask pixels
        after the best-fit model has been subtracted from the data. PSF_CHI2
        delivers the best-fit chi-squred statistic for each observation.
    plot : bool
        Plot fit results to the screen?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepprfphot kplr012557548-2012004120508_lpd-targ.fits.gz --columns 95
          --rows 1020 --fluxes 1.0 --border 0 --prfdir ../kplr2011265_prf --xtol 1e-7 --ftol 1e-7
          --plot --verbose

          --------------------------------------------------------------
          KEPPRFPHOT --  infile=kplr012557548-2012004120508_lpd-targ.fits.gz
          columns=95 rows=1020 fluxes=1.0 border=0 background=False
          focus=False prfdir=../kplr2011265_prf ranges=0,0 xtol=1e-07 ftol=1e-07
          qualflags=False plot=True overwrite=True verbose=True logfile=kepprfphot.log

          KEPPRFPHOT started at: Wed Jun 14 15:33:30 2017

                KepID: 12557548
           RA (J2000): 290.96622
          Dec (J2000): 51.50472
               KepMag: 15.692
             SkyGroup: 4
               Season: 1
              Channel: 32
               Module: 10
               Output: 4

           19% nrow = 740 t = 0.1 sec

    .. image:: ../_static/images/api/kepprfphot.png
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPPRFPHOT -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' columns={}'.format(columns)
            + ' rows={}'.format(rows)
            + ' fluxes={}'.format(fluxes)
            + ' border={}'.format(border)
            + ' background={}'.format(background)
            + ' focus={}'.format(focus)
            + ' prfdir={}'.format(prfdir)
            + ' ranges={}'.format(ranges)
            + ' xtol={}'.format(xtol)
            + ' ftol={}'.format(ftol)
            + ' qualflags={}'.format(qualflags)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPPRFPHOT started at', logfile, verbose)

    f = fluxes
    x = columns
    y = rows
    nsrc = len(f)

    if len(x) != nsrc or len(y) != nsrc:
        errmsg = ("ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and "
                  "fluxes must have the same number of sources")
        kepmsg.err(logfile, errmsg, verbose)

    guess = list(f) + list(x) + list(y)

    if background:
        if border == 0:
            guess.append(0.0)
        else:
            for i in range((border + 1) * 2):
                guess.append(0.0)
    if focus:
        guess = guess + [1.0, 1.0, 0.0]

    # overwrite output file
    for i in range(nsrc):
        outfile = '{0}_{1}.fits'.format(outfile, i)
        if overwrite:
            kepio.overwrite(outfile, logfile, verbose)
        if kepio.fileexists(outfile):
            errmsg = 'ERROR -- KEPPRFPHOT: {} exists. Use --overwrite'.format(outfile)
            kepmsg.err(logfile, errmsg, verbose)

    # open TPF FITS file
    try:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, barytime = \
            kepio.readTPF(infile, 'TIME', logfile, verbose)
    except:
        message = 'ERROR -- KEPPRFPHOT: is %s a Target Pixel File? ' % infile
        kepmsg.err(logfile,message,verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, tcorr = \
        kepio.readTPF(infile,'TIMECORR', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, cadno = \
        kepio.readTPF(infile,'CADENCENO',logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, fluxpixels = \
        kepio.readTPF(infile,'FLUX', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, errpixels = \
        kepio.readTPF(infile,'FLUX_ERR', logfile, verbose)
    try:
        kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, poscorr1 = \
            kepio.readTPF(infile, 'POS_CORR1', logfile, verbose)
    except:
        poscorr1 = np.zeros((len(barytime)), dtype='float32')
        poscorr1[:] = np.nan
    try:
        kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, poscorr2 = \
            kepio.readTPF(infile, 'POS_CORR2', logfile, verbose)
    except:
        poscorr2 = np.zeros((len(barytime)), dtype='float32')
        poscorr2[:] = np.nan
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, qual = \
        kepio.readTPF(infile,'QUALITY',logfile,verbose)
    struct = pyfits.open(infile)
    tstart, tstop, bjdref, cadence = kepio.timekeys(struct, infile, logfile, verbose)

    # input file keywords and mask map
    cards0 = struct[0].header.cards
    cards1 = struct[1].header.cards
    cards2 = struct[2].header.cards
    maskmap = np.copy(struct[2].data)
    npix = np.size(np.nonzero(maskmap)[0])

    # print target data
    if verbose:
        print('')
        print('      KepID: {}'.format(kepid))
        print(' RA (J2000): {}'.format(ra))
        print('Dec (J2000): {}'.format(dec))
        print('     KepMag: {}'.format(kepmag))
        print('   SkyGroup: {}'.format(skygroup))
        print('     Season: {}'.format(season))
        print('    Channel: {}'.format(channel))
        print('     Module: {}'.format(module))
        print('     Output: {}'.format(output))
        print('')

    # read PRF file and interpolate
    result = kepfunc.read_and_interpolate_prf(prfdir=prfdir, module=module,
                                              output=output, column=column,
                                              row=row, xdim=xdim, ydim=ydim,
                                              verbose=verbose, logfile=logfile)
    splineInterpolation = result[0]
    DATx = result[1]
    DATy = result[2]
    PRFx = result[4]
    PRFy = result[5]

    # construct mesh for background model
    bx = np.arange(1., float(xdim + 1))
    by = np.arange(1., float(ydim + 1))
    xx, yy = np.meshgrid(np.linspace(bx.min(), bx.max(), xdim),
                         np.linspace(by.min(), by.max(), ydim))
    # Get time ranges for new photometry, flag good data
    barytime += bjdref
    tstart, tstop = kepio.timeranges(ranges, logfile, verbose)
    incl = np.zeros((len(barytime)), dtype='int')
    for rownum in range(len(barytime)):
        for winnum in range(len(tstart)):
            if (barytime[rownum] >= tstart[winnum]
                and barytime[rownum] <= tstop[winnum]
                and (qual[rownum] == 0 or qualflags)
                and np.isfinite(barytime[rownum])
                and np.isfinite(np.nansum(fluxpixels[rownum, :]))):
                incl[rownum] = 1
    if not np.in1d(1,incl):
        message = ('ERROR -- KEPPRFPHOT: No legal data within the'
                   ' range {}'.format(ranges))
        kepmsg.err(logfile, message, verbose)
    # filter out bad data
    n = 0
    nincl = (incl == 1).sum()
    tim = np.zeros((nincl), 'float64')
    tco = np.zeros((nincl), 'float32')
    cad = np.zeros((nincl), 'float32')
    flu = np.zeros((nincl, len(fluxpixels[0])), 'float32')
    fer = np.zeros((nincl, len(fluxpixels[0])), 'float32')
    pc1 = np.zeros((nincl), 'float32')
    pc2 = np.zeros((nincl), 'float32')
    qua = np.zeros((nincl), 'float32')
    for rownum in range(len(barytime)):
        if incl[rownum] == 1:
            tim[n] = barytime[rownum]
            tco[n] = tcorr[rownum]
            cad[n] = cadno[rownum]
            flu[n,:] = fluxpixels[rownum]
            fer[n,:] = errpixels[rownum]
            pc1[n] = poscorr1[rownum]
            pc2[n] = poscorr2[rownum]
            qua[n] = qual[rownum]
            n += 1
    barytime = tim * 1.0
    tcorr = tco * 1.0
    cadno = cad * 1.0
    fluxpixels = flu * 1.0
    errpixels = fer * 1.0
    poscorr1 = pc1 * 1.0
    poscorr2 = pc2 * 1.0
    qual = qua * 1.0

    # initialize plot arrays
    t = np.array([], dtype='float64')
    fl, dx, dy, bg, fx, fy, fa, rs, ch = [], [], [], [], [], [], [], [], []
    for i in range(nsrc):
        fl.append(np.array([], dtype='float32'))
        dx.append(np.array([], dtype='float32'))
        dy.append(np.array([], dtype='float32'))
    # Preparing fit data message
    progress = np.arange(nincl)
    if verbose:
        txt  = 'Preparing...'
        sys.stdout.write(txt)
        sys.stdout.flush()
    # single processor version
    oldtime = 0.0
    for rownum in range(np.min([80, len(barytime)])):
        try:
            if barytime[rownum] - oldtime > 0.5:
                ftol = 1.0e-10; xtol = 1.0e-10
        except:
            pass
        args = (fluxpixels[rownum, :], errpixels[rownum, :], DATx, DATy, nsrc,
                border, xx, yy, PRFx, PRFy, splineInterpolation, guess, ftol,
                xtol, focus, background, rownum, 80, float(x[i]),
                float(y[i]), False)
        guess = PRFfits(args)
        ftol = ftol
        xtol = xtol
        oldtime = barytime[rownum]
    # Fit the time series: multi-processing
    anslist = []
    cad1 = 0
    cad2 = 50
    for i in range(int(nincl/50) + 1):
        try:
            fluxp = fluxpixels[cad1:cad2, :]
            errp = errpixels[cad1:cad2, :]
            progress = np.arange(cad1, cad2)
        except:
            fluxp = fluxpixels[cad1:nincl, :]
            errp = errpixels[cad1:nincl, :]
            progress = np.arange(cad1, nincl)
        try:
            args = itertools.izip(fluxp, errp, itertools.repeat(DATx),
                                  itertools.repeat(DATy),
                                  itertools.repeat(nsrc),
                                  itertools.repeat(border),
                                  itertools.repeat(xx),
                                  itertools.repeat(yy),
                                  itertools.repeat(PRFx),
                                  itertools.repeat(PRFy),
                                  itertools.repeat(splineInterpolation),
                                  itertools.repeat(guess),
                                  itertools.repeat(ftol),
                                  itertools.repeat(xtol),
                                  itertools.repeat(focus),
                                  itertools.repeat(background), progress,
                                  itertools.repeat(np.arange(cad1,nincl)[-1]),
                                  itertools.repeat(float(x[0])),
                                  itertools.repeat(float(y[0])),
                                  itertools.repeat(True))
            p = multiprocessing.Pool()
            model = [0.0]
            model = p.imap(PRFfits, args, chunksize=1)
            p.close()
            p.join()
            cad1 += 50; cad2 += 50
            ans = np.array([np.array(item) for item in zip(*model)])
            try:
                anslist = np.concatenate((anslist, ans.transpose()), axis=0)
            except:
                anslist = ans.transpose()
            guess = anslist[-1]
            ans = anslist.transpose()
        except:
            pass

    # single processor version
    oldtime = 0.0; ans = []
    for rownum in range(nincl):
        proctime = time.time()
        try:
            if barytime[rownum] - oldtime > 0.5:
                ftol = 1.0e-10; xtol = 1.0e-10
        except:
            pass
        args = (fluxpixels[rownum, :], errpixels[rownum, :], DATx, DATy, nsrc,
                border, xx, yy, PRFx, PRFy, splineInterpolation, guess, ftol,
                xtol, focus, background, rownum, nincl, float(x[0]),
                float(y[0]), True)
        guess = PRFfits(args)
        ans.append(guess)
        ftol = ftol; xtol = xtol; oldtime = barytime[rownum]
    ans = np.array(ans).transpose()

    # unpack the best fit parameters
    flux, OBJx, OBJy = [], [], []
    na = np.shape(ans)[1]
    for i in range(nsrc):
        flux.append(ans[i, :])
        OBJx.append(ans[nsrc + i, :])
        OBJy.append(ans[nsrc * 2 + i, :])
    try:
        bterms = border + 1
        if bterms == 1:
            b = ans[nsrc * 3, :]
        else:
            b = np.array([])
            bkg = []
            for i in range(na):
                bcoeff = np.array([ans[nsrc * 3:nsrc * 3 + bterms, i],
                                   ans[nsrc * 3 + bterms:nsrc * 3 + bterms * 2, i]])
                bkg.append(kepfunc.polyval2d(xx, yy, bcoeff))
                b = np.append(b, np.nanmean(bkg[-1].reshape(bkg[-1].size)))
    except:
        b = np.zeros(na)
    if focus:
        wx = ans[-3, :]
        wy = ans[-2, :]
        angle = ans[-1, :]
    else:
        wx = np.ones(na)
        wy = np.ones(na)
        angle = np.zeros(na)

    # constuct model PRF in detector coordinates
    residual, chi2 = [], []
    for i in range(na):
        f = np.empty(nsrc)
        x = np.empty(nsrc)
        y = np.empty(nsrc)
        for j in range(nsrc):
            f[j] = flux[j][i]
            x[j] = OBJx[j][i]
            y[j] = OBJy[j][i]
        PRFfit = kepfunc.PRF2DET(f, x, y, DATx, DATy, wx[i], wy[i], angle[i],
                                 splineInterpolation)
        if background and bterms == 1:
            PRFfit = PRFfit + b[i]
        if background and bterms > 1:
            PRFfit = PRFfit + bkg[i]

        # calculate residual of DATA - FIT
        xdim = np.shape(xx)[1]
        ydim = np.shape(yy)[0]
        DATimg = np.empty((ydim, xdim))
        n = 0
        for k in range(ydim):
            for j in range(xdim):
                DATimg[k,j] = fluxpixels[i, n]
                n += 1
        PRFres = DATimg - PRFfit
        residual.append(np.nansum(PRFres) / npix)
        # calculate the sum squared difference between data and model
        chi2.append(abs(np.nansum(np.square(DATimg - PRFfit) / PRFfit)))
    # load the output arrays
    otime = barytime - bjdref
    otimecorr = tcorr
    ocadenceno = cadno
    opos_corr1 = poscorr1
    opos_corr2 = poscorr2
    oquality = qual
    opsf_bkg = b
    opsf_focus1 = wx
    opsf_focus2 = wy
    opsf_rotation = angle
    opsf_residual = residual
    opsf_chi2 = chi2
    opsf_flux_err = np.empty((na))
    opsf_flux_err.fill(np.nan)
    opsf_centr1_err = np.empty((na))
    opsf_centr1_err.fill(np.nan)
    opsf_centr2_err = np.empty((na))
    opsf_centr2_err.fill(np.nan)
    opsf_bkg_err = np.empty((na))
    opsf_bkg_err.fill(np.nan)
    opsf_flux, opsf_centr1, opsf_centr2 = [], [], []
    for i in range(nsrc):
        opsf_flux.append(flux[i])
        opsf_centr1.append(OBJx[i])
        opsf_centr2.append(OBJy[i])

    # load the plot arrays
    t = barytime
    for i in range(nsrc):
        fl[i] = flux[i]
        dx[i] = OBJx[i]
        dy[i] = OBJy[i]
    bg = b
    fx = wx
    fy = wy
    fa = angle
    rs = residual
    ch = chi2

    # construct output primary extension
    for j in range(nsrc):
        hdu0 = pyfits.PrimaryHDU()
        for i in range(len(cards0)):
            if cards0[i].keyword not in hdu0.header.keys():
                hdu0.header[cards0[i].keyword] = (cards0[i].value,
                                                  cards0[i].comment)
            else:
                hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
        kepkey.history(call, hdu0, outfile, logfile, verbose)
        outstr = pyfits.HDUList(hdu0)

        # construct output light curve extension
        col1 = pyfits.Column(name='TIME', format='D', unit='BJD - 2454833',
                             array=otime)
        col2 = pyfits.Column(name='TIMECORR', format='E', unit='d',
                             array=otimecorr)
        col3 = pyfits.Column(name='CADENCENO', format='J', array=ocadenceno)
        col4 = pyfits.Column(name='PSF_FLUX', format='E', unit='e-/s',
                             array=opsf_flux[j])
        col5 = pyfits.Column(name='PSF_FLUX_ERR', format='E', unit='e-/s',
                             array=opsf_flux_err)
        col6 = pyfits.Column(name='PSF_BKG', format='E', unit='e-/s/pix',
                             array=opsf_bkg)
        col7 = pyfits.Column(name='PSF_BKG_ERR', format='E', unit='e-/s',
                             array=opsf_bkg_err)
        col8 = pyfits.Column(name='PSF_CENTR1', format='E', unit='pixel',
                             array=opsf_centr1[j])
        col9 = pyfits.Column(name='PSF_CENTR1_ERR', format='E', unit='pixel',
                             array=opsf_centr1_err)
        col10 = pyfits.Column(name='PSF_CENTR2', format='E', unit='pixel',
                              array=opsf_centr2[j])
        col11 = pyfits.Column(name='PSF_CENTR2_ERR', format='E', unit='pixel',
                              array=opsf_centr2_err)
        col12 = pyfits.Column(name='PSF_FOCUS1', format='E', array=opsf_focus1)
        col13 = pyfits.Column(name='PSF_FOCUS2', format='E', array=opsf_focus2)
        col14 = pyfits.Column(name='PSF_ROTATION', format='E', unit='deg',
                              array=opsf_rotation)
        col15 = pyfits.Column(name='PSF_RESIDUAL', format='E', unit='e-/s',
                              array=opsf_residual)
        col16 = pyfits.Column(name='PSF_CHI2', format='E', array=opsf_chi2)
        col17 = pyfits.Column(name='POS_CORR1', format='E', unit='pixel',
                              array=opos_corr1)
        col18 = pyfits.Column(name='POS_CORR2', format='E', unit='pixel',
                              array=opos_corr2)
        col19 = pyfits.Column(name='SAP_QUALITY', format='J', array=oquality)
        cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8,
                               col9, col10, col11, col12, col13, col14, col15,
                               col16, col17, col18, col19])
        hdu1 = pyfits.BinTableHDU.from_columns(cols)
        for i in range(len(cards1)):
            if (cards1[i].keyword not in hdu1.header.keys()
                and cards1[i].keyword[:4] not in ['TTYP', 'TFOR', 'TUNI',
                                                  'TDIS', 'TDIM', 'WCAX',
                                                  '1CTY', '2CTY', '1CRP',
                                                  '2CRP', '1CRV', '2CRV',
                                                  '1CUN', '2CUN', '1CDE',
                                                  '2CDE', '1CTY', '2CTY',
                                                  '1CDL', '2CDL', '11PC',
                                                  '12PC', '21PC', '22PC']):
                hdu1.header[cards1[i].keyword] = (cards1[i].value,
                                                  cards1[i].comment)
        outstr.append(hdu1)

        # construct output mask bitmap extension
        hdu2 = pyfits.ImageHDU(maskmap)
        for i in range(len(cards2)):
            if cards2[i].keyword not in hdu2.header.keys():
                hdu2.header[cards2[i].keyword] = (cards2[i].value,
                                                  cards2[i].comment)
            else:
                hdu2.header.cards[cards2[i].keyword].comment = cards2[i].comment
        outstr.append(hdu2)
        # write output file
        print("Writing output file {}...\n".format(outfile + '_' + str(j) + '.fits'))
        outstr.writeto(outfile + '_' + str(j) + '.fits', checksum=True)
        # close input structure
        struct.close()

    # clean up x-axis unit
    barytime0 = float(int(t[0] / 100) * 100.0)
    t -= barytime0
    t = np.insert(t,[0],[t[0]])
    t = np.append(t,[t[-1]])
    xlab = 'BJD $-$ %d' % barytime0

    # plot the light curves
    bg = np.insert(bg, [0], [-1.0e10])
    bg = np.append(bg, -1.0e10)
    fx = np.insert(fx, [0], [fx[0]])
    fx = np.append(fx, fx[-1])
    fy = np.insert(fy, [0], [fy[0]])
    fy = np.append(fy, fy[-1])
    fa = np.insert(fa, [0], [fa[0]])
    fa = np.append(fa, fa[-1])
    rs = np.insert(rs, [0], [-1.0e10])
    rs = np.append(rs, -1.0e10)
    ch = np.insert(ch, [0], [-1.0e10])
    ch = np.append(ch, -1.0e10)
    for i in range(nsrc):
        # clean up y-axis units
        nrm = math.ceil(math.log10(np.nanmax(fl[i]))) - 1.0
        fl[i] /= 10 ** nrm
        if nrm == 0:
            ylab1 = 'e$^-$ s$^{-1}$'
        else:
            ylab1 = '10$^{%d}$ e$^-$ s$^{-1}$' % nrm
        xx = np.copy(dx[i])
        yy = np.copy(dy[i])
        ylab2 = 'offset (pixels)'
        # data limits
        xmin = np.nanmin(t)
        xmax = np.nanmax(t)
        ymin1 = np.nanmin(fl[i])
        ymax1 = np.nanmax(fl[i])
        ymin2 = np.nanmin(xx)
        ymax2 = np.nanmax(xx)
        ymin3 = np.nanmin(yy)
        ymax3 = np.nanmax(yy)
        ymin4 = np.nanmin(bg[1:-1])
        ymax4 = np.nanmax(bg[1:-1])
        ymin5 = np.nanmin([np.nanmin(fx), np.nanmin(fy)])
        ymax5 = np.nanmax([np.nanmax(fx), np.nanmax(fy)])
        ymin6 = np.nanmin(fa[1:-1])
        ymax6 = np.nanmax(fa[1:-1])
        ymin7 = np.nanmin(rs[1:-1])
        ymax7 = np.nanmax(rs[1:-1])
        ymin8 = np.nanmin(ch[1:-1])
        ymax8 = np.nanmax(ch[1:-1])
        xr = xmax - xmin
        yr1 = ymax1 - ymin1
        yr2 = ymax2 - ymin2
        yr3 = ymax3 - ymin3
        yr4 = ymax4 - ymin4
        yr5 = ymax5 - ymin5
        yr6 = ymax6 - ymin6
        yr7 = ymax7 - ymin7
        yr8 = ymax8 - ymin8
        fl[i] = np.insert(fl[i], [0], [0.0])
        fl[i] = np.append(fl[i], 0.0)
        # define size of plot on monitor screen
        plt.figure(str(i + 1) + ' ' + str(time.asctime(time.localtime())),
                   figsize=[12,16])
        # delete any fossil plots in the matplotlib window
        plt.clf()
        # position first axes inside the plotting window
        ax = plt.axes([0.11, 0.523, 0.78, 0.45])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        # nox-label
        plt.setp(plt.gca(), xticklabels=[])
        # plot flux vs time
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for j in range(1, len(t)-1):
            dt = t[j] - t[j-1]
            if dt < work1:
                ltime = np.append(ltime, t[j])
                ldata = np.append(ldata, fl[i][j])
            else:
                plt.plot(ltime, ldata, color='#0000ff', linestyle='-',
                         linewidth=1.0)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        plt.plot(ltime, ldata, color='#0000ff', linestyle='-', linewidth=1.0)
        # plot the fill color below data time series, with no data gaps
        plt.fill(t,fl[i],fc='#ffff00',linewidth=0.0,alpha=0.2)
        # define plot x and y limits
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin1 - yr1 * 0.01 <= 0.0:
            plt.ylim(1.0e-10, ymax1 + yr1 * 0.01)
        else:
            plt.ylim(ymin1 - yr1 * 0.01, ymax1 + yr1 * 0.01)

        plt.ylabel('Source (' + ylab1 + ')', {'color' : 'k'})

        # make grid on plot
        plt.grid()

        # plot centroid tracks - position second axes inside the plotting window
        if focus and background:
            axs = [0.11, 0.433, 0.78, 0.09]
        elif background or focus:
            axs = [0.11, 0.388, 0.78, 0.135]
        else:
            axs = [0.11, 0.253, 0.78, 0.27]
        ax1 = plt.axes(axs)

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.setp(plt.gca(),xticklabels=[])

        # plot dx vs time
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for j in range(1, len(t)-1):
            dt = t[j] - t[j-1]
            if dt < work1:
                ltime = np.append(ltime, t[j])
                ldata = np.append(ldata, xx[j-1])
            else:
                ax1.plot(ltime, ldata, color='r', linestyle='-', linewidth=1.0)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        ax1.plot(ltime, ldata, color='r', linestyle='-', linewidth=1.0)

        # define plot x and y limits
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        plt.ylim(ymin2 - yr2 * 0.03, ymax2 + yr2 * 0.03)

        # plot labels
        ax1.set_ylabel('X-' + ylab2, color='k', fontsize=11)

        # position second axes inside the plotting window
        ax2 = ax1.twinx()

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.setp(plt.gca(), xticklabels=[])

        # plot dy vs time
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for j in range(1, len(t)-1):
            dt = t[j] - t[j-1]
            if dt < work1:
                ltime = np.append(ltime, t[j])
                ldata = np.append(ldata, yy[j-1])
            else:
                ax2.plot(ltime, ldata, color='g', linestyle='-', linewidth=1.0)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        ax2.plot(ltime, ldata, color='g', linestyle='-', linewidth=1.0)

        # define plot y limits
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        plt.ylim(ymin3 - yr3 * 0.03, ymax3 + yr3 * 0.03)

        # plot labels
        ax2.set_ylabel('Y-' + ylab2, color='k',fontsize=11)

        # background - position third axes inside the plotting window
        if background and focus:
            axs = [0.11, 0.343, 0.78, 0.09]
        if background and not focus:
            axs = [0.11, 0.253, 0.78, 0.135]
        if background:
            ax1 = plt.axes(axs)

            # force tick labels to be absolute rather than relative
            plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.setp(plt.gca(), xticklabels=[])

            # plot background vs time
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1, len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = np.append(ltime, t[j])
                    ldata = np.append(ldata, bg[j])
                else:
                    ax1.plot(ltime, ldata, color='#0000ff', linestyle='-',
                             linewidth=1.0)
                    ltime = np.array([], dtype='float64')
                    ldata = np.array([], dtype='float32')
            ax1.plot(ltime, ldata, color='#0000ff', linestyle='-',
                     linewidth=1.0)

            # plot the fill color below data time series, with no data gaps
            plt.fill(t, bg, fc='#ffff00', linewidth=0.0, alpha=0.2)

            # define plot x and y limits
            plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            plt.ylim(ymin4 - yr4 * 0.03, ymax4 + yr4 * 0.03)

            # plot labels
            ax1.set_ylabel('Background \n(e$^-$ s$^{-1}$ pix$^{-1}$)',
                           multialignment='center', color='k',fontsize=11)
            plt.grid()

        # position focus axes inside the plotting window
        if focus and background:
            axs = [0.11, 0.253, 0.78, 0.09]
        if focus and not background:
            axs = [0.11, 0.253, 0.78, 0.135]
        if focus:
            ax1 = plt.axes(axs)

            # force tick labels to be absolute rather than relative
            plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.setp(plt.gca(), xticklabels=[])

            # plot x-axis PSF width vs time
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = np.append(ltime, t[j])
                    ldata = np.append(ldata, fx[j])
                else:
                    ax1.plot(ltime, ldata, color='r', linestyle='-',
                             linewidth=1.0)
                    ltime = np.array([], dtype='float64')
                    ldata = np.array([], dtype='float32')
            ax1.plot(ltime, ldata, color='r', linestyle='-', linewidth=1.0)

            # plot y-axis PSF width vs time
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = np.append(ltime, t[j])
                    ldata = np.append(ldata, fy[j])
                else:
                    ax1.plot(ltime, ldata, color='g', linestyle='-',
                             linewidth=1.0)
                    ltime = np.array([], dtype='float64')
                    ldata = np.array([], dtype='float32')
            ax1.plot(ltime, ldata, color='g', linestyle='-', linewidth=1.0)

            # define plot x and y limits
            plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            plt.ylim(ymin5 - yr5 * 0.03, ymax5 + yr5 * 0.03)

            # plot labels
            ax1.set_ylabel('Pixel Scale\nFactor',
                           multialignment='center', color='k',fontsize=11)

            # Focus rotation - position second axes inside the plotting window
            ax2 = ax1.twinx()

            # force tick labels to be absolute rather than relative
            plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.setp(plt.gca(), xticklabels=[])

            # plot dy vs time
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = np.append(ltime, t[j])
                    ldata = np.append(ldata, fa[j])
                else:
                    ax2.plot(ltime, ldata, color='#000080', linestyle='-',
                             linewidth=1.0)
                    ltime = np.array([], dtype='float64')
                    ldata = np.array([], dtype='float32')
            ax2.plot(ltime, ldata, color='#000080', linestyle='-',
                     linewidth=1.0)

            # define plot y limits
            plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            plt.ylim(ymin6 - yr6 * 0.03, ymax6 + yr6 * 0.03)

            # plot labels
            ax2.set_ylabel('Rotation (deg)', color='k',fontsize=11)

        # fit residuals - position fifth axes inside the plotting window
        axs = [0.11, 0.163, 0.78, 0.09]
        ax1 = plt.axes(axs)

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.setp(plt.gca(), xticklabels=[])

        # plot residual vs time
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for j in range(1, len(t)-1):
            dt = t[j] - t[j-1]
            if dt < work1:
                ltime = np.append(ltime, t[j])
                ldata = np.append(ldata, rs[j])
            else:
                ax1.plot(ltime, ldata, color='b', linestyle='-', linewidth=1.0)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        ax1.plot(ltime, ldata, color='b', linestyle='-', linewidth=1.0)

        # plot the fill color below data time series, with no data gaps
        plt.fill(t, rs, fc='#ffff00', linewidth=0.0, alpha=0.2)

        # define plot x and y limits
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        plt.ylim(ymin7 - yr7 * 0.03, ymax7 + yr7 * 0.03)

        # plot labels
        ax1.set_ylabel('Residual \n(e$^-$ s$^{-1}$)',
                       multialignment='center', color='k', fontsize=11)

        plt.grid()

        # fit chi square - position sixth axes inside the plotting window
        axs = [0.11, 0.073, 0.78, 0.09]
        ax1 = plt.axes(axs)

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # plot background vs time
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for j in range(1,len(t)-1):
            dt = t[j] - t[j-1]
            if dt < work1:
                ltime = np.append(ltime, t[j])
                ldata = np.append(ldata, ch[j])
            else:
                ax1.plot(ltime, ldata, color='b', linestyle='-', linewidth=1.0)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        ax1.plot(ltime, ldata, color='b', linestyle='-', linewidth=1.0)

        # plot the fill color below data time series, with no data gaps
        plt.fill(t, ch, fc='#ffff00', linewidth=0.0, alpha=0.2)

        # define plot x and y limits
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        plt.ylim(ymin8 - yr8 * 0.03, ymax8 + yr8 * 0.03)

        # plot labels
        ax1.set_ylabel('$\chi^2$ (%d dof)' % (npix - len(guess) - 1),
                       color='k', fontsize=11)
        plt.xlabel(xlab, {'color' : 'k'})

        # make grid on plot
        plt.grid()

        # render plot
        plt.savefig(outfile + '_' + str(i) + '.png')
        plt.show()

    # stop time
    kepmsg.clock('\n\nKEPPRFPHOT ended at',logfile,verbose)

def PRFfits(args):

    # start time
    proctime = time.time()

    # extract image from the time series
    xdim = np.shape(args[6])[1]
    ydim = np.shape(args[6])[0]
    DATimg = np.empty((ydim,xdim))
    DATerr = np.empty((ydim,xdim))
    n = 0
    for i in range(ydim):
        for j in range(xdim):
            DATimg[i,j] = args[0][n]
            DATerr[i,j] = args[1][n]
            n += 1

    # minimize data and model
    if args[14] and args[15]:
        argm = (args[2], args[3], DATimg, DATerr, args[4], args[5], args[6],
                args[7], args[10], args[18], args[19])
        ans = fmin_powell(kepfunc.PRFwithFocusAndBackground, args[11],
                          args=argm, xtol=args[12], ftol=args[13], disp=False)
    elif args[14] and not args[15]:
        argm = (args[2], args[3], DATimg, DATerr, args[4], args[10], args[18],
                args[19])
        ans = fmin_powell(kepfunc.PRFwithFocus, args[11], args=argm,
                          xtol=args[12], ftol=args[13], disp=False)
    elif args[15] and not args[14]:
        argm = (args[2], args[3], DATimg, DATerr, args[4], args[5], args[6],
                args[7], args[10], args[18], args[19])
        ans = fmin_powell(kepfunc.PRFwithBackground, args[11], args=argm,
                          xtol=args[12], ftol=args[13], disp=False)
    else:
        argm = (args[2], args[3], DATimg, DATerr, args[4], args[10], args[18],
                args[19])
        ans = fmin_powell(kepfunc.PRF, args[11], args=argm, xtol=args[12],
                          ftol=args[13], disp=False)

    # print progress
    if args[20]:
        txt  = '\r%3d%% ' % ((float(args[16]) + 1.0) / float(args[17]) * 100.0)
        txt += 'nrow = %d ' % (args[16]+1)
        txt += 't = %.1f sec' % (time.time() - proctime)
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return ans

def kepprfphot_main():

    import argparse

    parser = argparse.ArgumentParser(
             description='Fitting PRF model to Target Pixel time series',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input target pixel file',
                        type=str)
    parser.add_argument('--prfdir',
                        help='Folder containing PRF files',
                        type=str)
    parser.add_argument('--columns',
                        help='Column number of each source to be fit',
                        nargs='+', type=float)
    parser.add_argument('--rows', help='Row number of each source to be fit',
                        nargs='+', type=float)
    parser.add_argument('--fluxes',
                        help='Relative flux of each source to be fit',
                        nargs='+', type=float)
    parser.add_argument('--border',
                        help='Order of background polynmial fit', default=0,
                        type=int)
    parser.add_argument('--background', action='store_true',
                        help='Fit background?')
    parser.add_argument('--focus', action='store_true',
                        help='Fit focus changes?')
    parser.add_argument('--ranges', default='0,0', help='Time ranges to fit',
                        type=str)
    parser.add_argument('--xtol', default=1.0e-4, help='Fit parameter xtol',
                        type=float)
    parser.add_argument('--ftol', default=1.0e-2,
                        help='Fit minimization tolerance', type=float)
    parser.add_argument('--qualflags', action='store_true',
                        help='Fit data that have quality flags?')
    parser.add_argument('--outfile',
                        help=('Root name of output light curve FITS files.'
                              ' If None, root name is infile-kepprfphot.'),
                        default=None)
    parser.add_argument('--plot', action='store_true',
                        help='Plot fit results?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepprfphot.log',
                        help='Name of ascii log file', type=str)
    args = parser.parse_args()
    kepprfphot(args.infile, args.prfdir, args.columns, args.rows, args.fluxes,
               args.border, args.background, args.focus, args.ranges,
               args.xtol, args.ftol, args.qualflags, args.outfile, args.plot,
               args.overwrite, args.verbose, args.logfile)
