import numpy as np
import time
import math
import glob
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy.io import fits as pyfits
from scipy import interpolate, optimize, ndimage, stats
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import interpolation
from . import kepio, kepmsg, kepplot, kepfunc, kepstat
from .utils import PyKEArgumentHelpFormatter


__all__ = ['kepprf']


def kepprf(infile, prfdir, frameno, columns, rows, fluxes, background=False,
           border=1, focus=False, xtol=1e-4, ftol=1., outfile=None, plot=False,
           imscale='linear', cmap='YlOrBr', apercol='#ffffff', verbose=False,
           logfile='kepprf.log'):
    """
    kepprf -- Fit a PSF model to a specific image within a Target Pixel File

    Fit a PSF model, combined with spacecraft jitter and pixel scale
    drift (the Pixel Response Function; PRF) to a single observation of Kepler
    target pixels.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler Target
        Pixel data within the first data extension.
    prfdir : str
        The full or relative directory path to a folder containing the Kepler
        PSF calibration. Calibration files can be downloaded from the Kepler
        focal plane characteristics page at the MAST.
    frameno : int
        The cadence number in the input file data containing the pixels to
        plot. If the chosen observation has a non-zero quality flag set or the
        pixel set contains only NULLs then the task will halt with an error
        message.
    columns : list
        A starting guess for the CCD column position(s) of the source(s) that
        are to be fit. The model is unlikely to converge if the guess is too
        far away from the correct location. A rule of thumb is to provide a
        guess within 1 CCD pixel of the true position. If more than one source
        is being modeled then the column positions of each are separated by a
        comma. The same number of sources in the columns, rows and fluxes
        field is a requirement of this task.
    rows : list
        A starting guess for the CCD row position(s) of the source(s) that are
        to be fit. The model is unlikely to converge if the guess is too far
        away from the correct location. A rule of thumb is to provide a guess
        within 1 CCD pixel of the true position. If more than one source is
        being modeled then the row positions of each are separated by a comma.
        The same number of sources in the columns, rows and fluxes field is a
        requirement of this task.
    fluxes : list
        A starting guess for the flux(es) of the source(s) that are to be fit.
        Fit convergence is not particularly reliant on the accuracy of these
        guesses, but the fit will converge faster the more accurate the guess.
        If more than one source is being modeled then the row positions of
        each are separated by a comma. The same number of sources in the
        columns, rows and fluxes field is a requirement of this task.
    background : boolean
        Whether to include a background component in the model. If `True` then
        the background will be represented by a two-dimensional polynomial of
        order `border`. This functionality is somewhat experimental, with one
        eye upon potential background gradients across large masks or on those
        detectors more prone to pattern noise. Generally it is recommended to
        set background as `False`.
    border : int
        If a background is included in the fit then it is modeled as a
        two-dimensional polynomial. This parameter is the polynomial order. A
        zero-order polynomial is generally recommended.
    focus : boolean
        Whether to incude pixel scale and focus rotation with the fit
        parameters of the model. This is also an experimental function. This
        approach does not attempt to deal with inter- or intra-pixel
        variations. The recommended use is currently to set focus as `False`.
    xtol : float
        The dimensionless, relative model parameter convergence criterion for
        the fit algorithm.
    ftol : float
        The dimensionless, relative model residual convergence criterion for
        the fit algorithm.
    outfile : str
        Name of an optional output plot file containing the results of kepprf.
        An example is provided in Figure 1. Typically this is a PNG format
        file. If no file is required, outfile can be 'None' or blank, in
        which case the plot will be generated but the plot will not be saved
        to a file. Any existing file with this name will be automatically
        overwritten.
    plot : boolean
        Plot fit results to the screen?
    imscale : str
        kepprf can plot images with three choices of image scales. The choice
        is made using this argument.
        The options are:
        * linear
        * logarithmic
        * squareroot
    cmap : str
        matplotlib's color map
    verbose : boolean
        Print informative messages and warnings to the shell and logfile?
    logfile : string
        Name of the logfile containing error and warning messages.

    Examples
    --------
    Using the command line tool ``kepprf``, one can fit multiple PRFs to a given
    frame in a target pixel file as follows

    .. code-block:: bash

        $ kepprf kplr008256049-2010174085026_lpd-targ.fits --prfdir ~/kplr2011265_prf/
        --frameno 1000 --columns 830 831 --rows 242 241 --fluxes 1.0 0.1 --plot --verbose

              KepID: 8256049
                BJD: 2455296.903574196
         RA (J2000): 298.67861
        Dec (J2000): 44.1755
             KepMag: 15.654
           SkyGroup: 53
             Season: 3
            Channel: 81
             Module: 24
             Output: 1

        Convergence time = 0.15390515327453613s

        Flux = 3978.040625752744 e-/s X = 829.8259431097927 pix Y = 242.3810334478628 pix
        Flux = 4734.069273790539 e-/s X = 830.990805551025 pix Y = 240.97340366638306 pix

                        Total flux in mask = 10747.293440638587 e-/s
                       Target flux in mask = 3793.041929468528 e-/s
                    Total flux in aperture = 6365.551487630484 e-/s
                   Target flux in aperture = 3110.924803570053 e-/s
          Target flux fraction in aperture = 78.2024392468689 %
        Contamination fraction in aperture = 51.12874650978488 %

               Residual flux = -0.5748827605745994 e-/s
        Pearsons chi^2 test = 296.12077907844986 for 13 dof
                 Chi^2 test = 19803.55879917441 for 13 dof

    .. image:: ../_static/images/api/kepprf.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.png".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPPRF -- '
            + ' infile={}'.format(infile)
            + ' frameno={}'.format(frameno)
            + ' columns={}'.format(columns)
            + ' rows={}'.format(rows)
            + ' fluxes={}'.format(fluxes)
            + ' prfdir={}'.format(prfdir)
            + ' background={}'.format(background)
            + ' border={}'.format(border)
            + ' focus={}'.format(focus)
            + ' xtol={}'.format(xtol)
            + ' ftol={}'.format(xtol)
            + ' outfile={}'.format(outfile)
            + ' plot={}'.format(plot)
            + ' imscale={}'.format(imscale)
            + ' cmap={}'.format(cmap)
            + ' apercol={}'.format(apercol)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call + '\n', verbose)

    # start time
    kepmsg.clock('KEPPRF started at',logfile,verbose)

    # construct inital guess vector for fit
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

    try:
        kepid, channel, skygroup, module, output, quarter, season, ra, \
        dec, column, row, kepmag, xdim, ydim, barytime = \
        kepio.readTPF(infile, 'TIME', logfile, verbose)
    except:
        errmsg = "ERROR -- KEPPRF: is {} a Target Pixel File? ".format(infile)
        kepmsg.err(logfile, errmsg, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, tcorr = \
        kepio.readTPF(infile, 'TIMECORR', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, cadno = \
        kepio.readTPF(infile, 'CADENCENO', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, fluxpixels = \
        kepio.readTPF(infile,'FLUX',logfile,verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, errpixels = \
        kepio.readTPF(infile, 'FLUX_ERR', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, qual = \
        kepio.readTPF(infile, 'QUALITY', logfile, verbose)

    # read mask defintion data from TPF file

    maskimg, pixcoord1, pixcoord2 = kepio.readMaskDefinition(infile, logfile, verbose)
    npix = np.size(np.nonzero(maskimg)[0])

    # print target data
    if verbose:
        print('')
        print('      KepID: {}'.format(kepid))
        print('        BJD: {}'.format(barytime[frameno-1] + 2454833.0))
        print(' RA (J2000): {}'.format(ra))
        print('Dec (J2000): {}'.format(dec))
        print('     KepMag: {}'.format(kepmag))
        print('   SkyGroup: {}'.format(skygroup))
        print('     Season: {}'.format(str(season)))
        print('    Channel: {}'.format(channel))
        print('     Module: {}'.format(module))
        print('     Output: {}'.format(output))
        print('')

    # is this a good row with finite timestamp and pixels?
    if (not np.isfinite(barytime[frameno-1])
        or np.nansum(fluxpixels[frameno-1,:]) == np.nan):
        errmsg = ("ERROR -- KEPFIELD: Row {0} is a bad quality timestamp"
                  .format(frameno))
        kepmsg.err(logfile, errmsg, verbose)

    # construct input pixel image
    flux = fluxpixels[frameno-1,:]
    ferr = errpixels[frameno-1,:]

    # image scale and intensity limits of pixel data
    n = 0
    DATimg = np.empty((ydim, xdim))
    ERRimg = np.empty((ydim, xdim))
    for i in range(ydim):
        for j in range(xdim):
            DATimg[i, j] = flux[n]
            ERRimg[i, j] = ferr[n]
            n += 1

    # read and interpolate PRF
    (splineInterpolation, DATx, DATy, prf, _, _, PRFx0, PRFy0, cdelt1p, 
        cdelt2p, prfDimX, prfDimY) = kepfunc.read_and_interpolate_prf(
        prfdir=prfdir, module=module, output=output, column=column, row=row, 
        xdim=xdim, ydim=ydim, verbose=verbose, logfile=logfile)

    # construct mesh for background model
    if background:
        bx = np.arange(1.,float(xdim+1))
        by = np.arange(1.,float(ydim+1))
        xx, yy = np.meshgrid(np.linspace(bx.min(), bx.max(), xdim),
                             np.linspace(by.min(), by.max(), ydim))

    # fit PRF model to pixel data
    start = time.time()
    if focus and background:
        args = (DATx, DATy, DATimg, ERRimg, nsrc, border, xx, yy,
                splineInterpolation, float(x[0]), float(y[0]))
        ans = fmin_powell(kepfunc.PRFwithFocusAndBackground, guess, args=args,
                          xtol=xtol, ftol=ftol, disp=False)
    elif focus and not background:
        args = (DATx, DATy, DATimg, ERRimg, nsrc, splineInterpolation,
                float(x[0]), float(y[0]))
        ans = fmin_powell(kepfunc.PRFwithFocus, guess, args=args, xtol=xtol,
                          ftol=ftol, disp=False)
    elif background and not focus:
        args = (DATx, DATy, DATimg, ERRimg, nsrc, border, xx, yy,
                splineInterpolation, float(x[0]), float(y[0]))
        ans = fmin_powell(kepfunc.PRFwithBackground, guess, args=args,
                          xtol=xtol, ftol=ftol, disp=False)
    else:
        args = (DATx, DATy, DATimg, ERRimg, nsrc, splineInterpolation,
                float(x[0]), float(y[0]))
        ans = fmin_powell(kepfunc.PRF, guess, args=args, xtol=xtol,
                          ftol=ftol, disp=False)
    print("Convergence time = {}s\n".format(time.time() - start))

    # pad the PRF data if the PRF array is smaller than the data array
    flux = []
    OBJx = []
    OBJy = []
    PRFmod = np.zeros((prfDimY, prfDimX))
    if PRFy0 < 0 or PRFx0 < 0.0:
        PRFmod = np.zeros((prfDimY, prfDimX))
        superPRF = np.zeros((prfDimY + 1, prfDimX + 1))
        superPRF[abs(PRFy0):abs(PRFy0) + np.shape(prf)[0],
                 abs(PRFx0):abs(PRFx0) + np.shape(prf)[1]] = prf
        prf = superPRF * 1.0
        PRFy0 = 0
        PRFx0 = 0
        # rotate the PRF model around its center
    if focus:
        angle = ans[-1]
        prf = interpolation.rotate(prf, -angle, reshape=False,
                                   mode='nearest')
    # iterate through the sources in the best fit PSF model
    for i in range(nsrc):
        flux.append(ans[i])
        OBJx.append(ans[nsrc + i])
        OBJy.append(ans[nsrc * 2 + i])
        # calculate best-fit model
        y = (OBJy[i] - np.mean(DATy)) / cdelt1p[0]
        x = (OBJx[i] - np.mean(DATx)) / cdelt2p[0]
        prfTmp = interpolation.shift(prf, [y, x], order=3, mode='constant')
        prfTmp = prfTmp[PRFy0:PRFy0 + prfDimY, PRFx0:PRFx0 + prfDimX]
        PRFmod = PRFmod + prfTmp * flux[i]
        wx = 1.0
        wy = 1.0
        angle = 0
        b = 0.0
        # write out best fit parameters
        if verbose:
            txt = ("Flux = {0} e-/s X = {1} pix Y = {2} pix"
                   .format(flux[i], OBJx[i], OBJy[i]))
            kepmsg.log(logfile, txt, True)
    if background:
        bterms = border + 1
        if bterms == 1:
            b = ans[nsrc * 3]
        else:
            bcoeff = np.array([ans[nsrc*3:nsrc*3+bterms],
                               ans[nsrc*3+bterms:nsrc*3+bterms*2]])
            bkg = kepfunc.polyval2d(xx, yy, bcoeff)
            b = np.nanmean(bkg.reshape(bkg.size))
        txt = "\n   Mean background = {0} e-/s".format(b)
        kepmsg.log(logfile, txt, True)
    if focus:
        wx = ans[-3]
        wy = ans[-2]
        angle = ans[-1]
    if focus:
        if not background:
            kepmsg.log(logfile, '', True)
        kepmsg.log(logfile, " X/Y focus factors = {0}/{1}".format(wx, wy),
                   True)
        kepmsg.log(logfile, "PRF rotation angle = {0} deg".format(angle),
                   True)
    # measure flux fraction and contamination
    PRFall = kepfunc.PRF2DET(flux, OBJx, OBJy, DATx, DATy, wx, wy, angle,
                             splineInterpolation)
    PRFone = kepfunc.PRF2DET([flux[0]], [OBJx[0]], [OBJy[0]], DATx, DATy,
                             wx, wy, angle, splineInterpolation)
    FluxInMaskAll = np.nansum(PRFall)
    FluxInMaskOne = np.nansum(PRFone)
    FluxInAperAll = 0.0
    FluxInAperOne = 0.0
    for i in range(1, ydim):
        for j in range(1, xdim):
            if kepstat.bitInBitmap(maskimg[i, j], 2):
                FluxInAperAll += PRFall[i, j]
                FluxInAperOne += PRFone[i, j]
    FluxFraction = FluxInAperOne / flux[0]
    try:
        Contamination = (FluxInAperAll - FluxInAperOne) / FluxInAperAll
    except:
        Contamination = 0.0

    kepmsg.log(logfile, "\n                Total flux in mask = {0} e-/s"
            .format(FluxInMaskAll), True)
    kepmsg.log(logfile, "               Target flux in mask = {0} e-/s"
            .format(FluxInMaskOne), True)
    kepmsg.log(logfile, "            Total flux in aperture = {0} e-/s"
            .format(FluxInAperAll), True)
    kepmsg.log(logfile, "           Target flux in aperture = {0} e-/s"
            .format(FluxInAperOne), True)
    kepmsg.log(logfile, "  Target flux fraction in aperture = {0} %"
            .format(FluxFraction * 100.0), True)
    kepmsg.log(logfile, "Contamination fraction in aperture = {0} %"
            .format(Contamination * 100.0), True)

    # construct model PRF in detector coordinates
    PRFfit = PRFall + 0.0
    if background and bterms == 1:
        PRFfit = PRFall + b
    if background and bterms > 1:
        PRFfit = PRFall + bkg

    # calculate residual of DATA - FIT
    PRFres = DATimg - PRFfit
    FLUXres = np.nansum(PRFres) / npix

    # calculate the sum squared difference between data and model
    Pearson = abs(np.nansum(np.square(DATimg - PRFfit) / PRFfit))
    Chi2 = np.nansum(np.square(DATimg - PRFfit) / np.square(ERRimg))
    DegOfFreedom = npix - len(guess) - 1
    try:
        kepmsg.log(logfile, "\n       Residual flux = {0} e-/s"
                .format(FLUXres), True)
        kepmsg.log(logfile, "Pearson\'s chi^2 test = {0} for {1} dof"
                .format(Pearson, DegOfFreedom), True)
    except:
        pass
    kepmsg.log(logfile, "          Chi^2 test = {0} for {1} dof"
            .format(Chi2, DegOfFreedom), True)

    # image scale and intensity limits for plotting images
    imgdat_pl, zminfl, zmaxfl = kepplot.intScale2D(DATimg, imscale)
    imgprf_pl, zminpr, zmaxpr = kepplot.intScale2D(PRFmod, imscale)
    imgfit_pl, zminfi, zmaxfi = kepplot.intScale2D(PRFfit, imscale)
    imgres_pl, zminre, zmaxre = kepplot.intScale2D(PRFres, 'linear')
    if imscale == 'linear':
        zmaxpr *= 0.9
    elif imscale == 'logarithmic':
        zmaxpr = np.max(zmaxpr)
        zminpr = zmaxpr / 2

    plt.figure(figsize=[12, 10])
    plt.clf()
    plotimage(imgdat_pl, zminfl, zmaxfl, 1, row, column, xdim, ydim, 0.07,
              0.53, 'observation', cmap)
    plotimage(imgprf_pl, zminpr, zmaxpr, 2, row, column, xdim, ydim, 0.44,
              0.53, 'model', cmap)
    kepplot.borders(maskimg, xdim, ydim, pixcoord1, pixcoord2, 1, apercol,
                    '--', 0.5)
    kepplot.borders(maskimg, xdim, ydim, pixcoord1, pixcoord2, 2, apercol,
                    '-', 3.0)
    plotimage(imgfit_pl, zminfl, zmaxfl, 3, row, column, xdim, ydim, 0.07,
              0.08, 'fit', cmap)
    plotimage(imgres_pl, zminfl, zmaxfl, 4, row, column, xdim, ydim, 0.44,
              0.08, 'residual', cmap)

    # plot data color bar
    barwin = plt.axes([0.84,0.08,0.06,0.9])
    if imscale == 'linear':
        brange = np.arange(zminfl, zmaxfl, (zmaxfl - zminfl) / 1000)
    elif imscale == 'logarithmic':
        brange = np.arange(10.0 ** zminfl, 10.0 ** zmaxfl,
                           (10.0 ** zmaxfl - 10.0**zminfl) / 1000)
    elif imscale == 'squareroot':
        brange = np.arange(zminfl**2, zmaxfl ** 2,
                           (zmaxfl**2 - zminfl**2) / 1000)
    if imscale == 'linear':
        barimg = np.resize(brange, (1000, 1))
    elif imscale == 'logarithmic':
        barimg = np.log10(np.resize(brange, (1000, 1)))
    elif imscale == 'squareroot':
        barimg = np.sqrt(np.resize(brange, (1000, 1)))
    try:
        nrm = len(str(int(np.nanmax(brange)))) - 1
    except:
        nrm = 0
    brange = brange / 10 ** nrm
    plt.imshow(barimg, aspect='auto', interpolation='nearest', origin='lower',
               vmin=np.nanmin(barimg), vmax=np.nanmax(barimg),
               extent=(0.0, 1.0, brange[0], brange[-1]), cmap=cmap)
    barwin.yaxis.tick_right()
    barwin.yaxis.set_label_position('right')
    barwin.yaxis.set_major_locator(plt.MaxNLocator(7))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().set_autoscale_on(False)
    plt.setp(plt.gca(), xticklabels=[], xticks=[])
    plt.ylabel('Flux (10$^{%d}$ e$^-$ s$^{-1}$)' % nrm)
    barwin.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

    # render plot
    print("Writing output file {}...".format(outfile))
    plt.savefig(outfile)
    if plot:
        plt.draw()
        plt.show()
    # stop time
    kepmsg.clock('\nKEPPRF ended at', logfile, verbose)

# plot channel image
def plotimage(imgflux_pl, zminfl, zmaxfl, plmode, row, column,
              xdim, ydim, winx, winy, tlabel, cmap):
    # pixel limits of the subimage
    ymin = row
    ymax = ymin + ydim
    xmin = column
    xmax = xmin + xdim
    # plot limits for flux image
    ymin = float(ymin) - 0.5
    ymax = float(ymax) - 0.5
    xmin = float(xmin) - 0.5
    xmax = float(xmax) - 0.5
    # plot the image window
    ax = plt.axes([winx,winy,0.37,0.45])
    plt.imshow(imgflux_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zminfl,vmax=zmaxfl,extent=(xmin,xmax,ymin,ymax),cmap=cmap)
    plt.gca().set_autoscale_on(False)
    labels = ax.get_yticklabels()
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    if plmode == 1:
        plt.setp(plt.gca(),xticklabels=[])
    if plmode == 2:
        plt.setp(plt.gca(),xticklabels=[],yticklabels=[])
    if plmode == 4:
        plt.setp(plt.gca(),yticklabels=[])
    if plmode == 3 or plmode == 4:
        plt.xlabel('Pixel Column Number', {'color' : 'k'})
    if plmode == 1 or plmode == 3:
        plt.ylabel('Pixel Row Number', {'color' : 'k'})
    plt.text(0.05, 0.93,tlabel,horizontalalignment='left',verticalalignment='center',
               fontsize=36,fontweight=500,transform=ax.transAxes)

def kepprf_main():
    import argparse

    parser = argparse.ArgumentParser(
             description="Fitting PRF model to Target Pixel image",
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input target pixel file',
                        type=str)
    parser.add_argument('--prfdir',
                        help=("Folder containing Point Response Function "
                              "FITS files"), type=str)
    parser.add_argument('--frameno',
                        help='Cadence number of image stored in infile',
                        type=int)
    parser.add_argument('--columns',
                        help=("Initial guesses for the center of each source "
                              "on the x-axis"),
                        nargs='+', type=float)
    parser.add_argument('--rows',
                        help=("Initial guesses for the center of each source "
                              "on the x-axis"),
                        nargs='+', type=float)
    parser.add_argument('--fluxes',
                        help='Relative flux of each source to be fit',
                        nargs='+', type=float)
    parser.add_argument('--background', action='store_true',
                        help='Fit background?')
    parser.add_argument('--border', help='Order of background polynmial fit',
                        default=1, type=int)
    parser.add_argument('--focus', action='store_true',
                        help='Fit focus changes?', default=False)
    parser.add_argument('--xtol', default=1.0e-4,
                        help='Fit parameter tolerance', dest='xtol',
                        type=float)
    parser.add_argument('--ftol', '-f', default=1.0,
                        help='Fit minimization tolerance', dest='ftol',
                        type=float)
    parser.add_argument('--outfile',
                        help=('Name of output PNG plot file.'
                              ' If None, outfile is infile-kepprf.'),
                        default=None)
    parser.add_argument('--plot', action='store_true',
                        help='Plot fit results?', default=False)
    parser.add_argument('--imscale', help='Type of image intensity scale',
                        default='linear', type=str,
                        choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', help='Image colormap', default='YlOrBr',
                        type=str)
    parser.add_argument('--apercol', help='Aperture color', default='#ffffff',
                        type=str)
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', default='kepprf.log',
                        help='Name of ascii log file', type=str)
    args = parser.parse_args()

    kepprf(args.infile, args.prfdir, args.frameno, args.columns,
           args.rows, args.fluxes, args.background, args.border, args.focus,
           args.xtol, args.ftol, args.outfile, args.plot, args.imscale, args.cmap,
           args.apercol, args.verbose, args.logfile)
