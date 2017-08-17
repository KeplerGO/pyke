from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepplot, kepstat, kepfunc


__all__ = ['keppixseries']


def keppixseries(infile, outfile=None, plotfile=None, plottype='global',
                 filterlc=False, function='boxcar', cutoff=1.0, overwrite=False,
                 verbose=False, logfile='keppixseries.log'):
    """
    keppixseries -- individual time series photometry for all pixels within a
    target mask

    keppixseries plots a light curve for each individual pixel in a target
    mask. Light curves are extracted from a target pixel file obtained from the
    Kepler data archive at MAST. If required, the data can be fed through a
    boxcar, gaussian or sinc function high bandpass filter in order to remove
    low frequency signal from the data. keppixseries is a diagnostic tool for
    identifying source contaminants in the background or foreground of the
    target. It can be employed to identify pixels for inclusion or exclusion
    when re-extracting a Kepler light curve from target pixel files.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler Target
        Pixel data within the first data extension.
    outfile : str
        The name of the output FITS file. This file has two data extensions.
        The first called 'PIXELSERIES' contains a table with columns of
        barycenter-corrected time, barycenter time correction, cadence number,
        cadence quality flag and a series of photometric light curves, one for
        each pixel within the target mask. Each pixel is labeled COLx_ROWy,
        where :math:`x` is the pixel column number and :math:`y` is the pixel
        row number on the CCD module/output. The second extension contains the
        mask definition map copied directly from the input target pixel file.
    plotfile : str
        Name of an optional diagnostic output plot file containing the results
        of keppixseries. An example is provided in Figure 1. Typically this is
        a PNG format file. If no diagnostic file is required, plotfile can be
        'None'. The plot will be generated regardless of the value of this
        field, but the plot will not be saved to a file if ``plotfile='None'``.
    plottype : str
        keppixseries can plot light curves of three types.
        The choice is made using this argument. The options are:

        * local - All individual pixel light curves are scaled separately to
          provide the most dynamic range for each pixel.
        * global - All pixel light curves are scaled between zero and the
          maximum flux attained by the brightest pixel in the mask. This option
          provides the relative contribution to the archived light curve by each
          pixel.
        * full - All pixels light curves are scaled between zero and the
          maximum flux attained by that pixel. This provides the fraction of
          variability within each individual pixel.
    filterlc : bool
        If True, the light curve for each pixel will be treated by a high
        band-pass filter to remove long-term trends from e.g. differential
        velocity aberration.
    function : str
        The functional form of the high pass-band filter:

        * boxcar
        * gauss
        * sinc
    cutoff : float
        The frequency of the high pass-band cutoff in units of :math:`days^{-1}`.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile = str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block :: bash

        $ keppixseries kplr008256049-2010174085026_lpd-targ.fits.gz

    .. image:: ../_static/images/api/keppixseries.png
        :align: center
    """
    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPPIXSERIES -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' plotfile={}'.format(plotfile)
            + ' plottype={}'.format(plottype)
            + ' filterlc={}'.format(filterlc)
            + ' function={}'.format(function)
            + ' cutoff={}'.format(cutoff)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPPIXSERIES started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPPIXSERIES: {} exists. Use --overwrite'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # open TPF FITS file
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, barytime = \
        kepio.readTPF(infile, 'TIME', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, tcorr = \
        kepio.readTPF(infile, 'TIMECORR', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, cadno = \
        kepio.readTPF(infile, 'CADENCENO', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, fluxpixels = \
        kepio.readTPF(infile, 'FLUX', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, errpixels = \
        kepio.readTPF(infile, 'FLUX_ERR', logfile, verbose)
    kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, qual = \
        kepio.readTPF(infile, 'QUALITY', logfile, verbose)

    # read mask defintion data from TPF file
    maskimg, pixcoord1, pixcoord2 = kepio.readMaskDefinition(infile, logfile,
                                                             verbose)
    # print target data
    print('')
    print('      KepID:  {}'.format(kepid))
    print(' RA (J2000):  {}'.format(ra))
    print('Dec (J2000): {}'.format(dec))
    print('     KepMag:  {}'.format(kepmag))
    print('   SkyGroup:    {}'.format(skygroup))
    print('     Season:    {}'.format(season))
    print('    Channel:    {}'.format(channel))
    print('     Module:    {}'.format(module))
    print('     Output:     {}'.format(output))
    print('')
    # how many quality = 0 rows?
    npts = 0
    nrows = len(fluxpixels)
    for i in range(nrows):
        if (qual[i] == 0 and np.isfinite(barytime[i])
            and np.isfinite(fluxpixels[i, ydim * xdim // 2])):
            npts += 1
    time = np.empty((npts))
    timecorr = np.empty((npts))
    cadenceno = np.empty((npts))
    quality = np.empty((npts))
    pixseries = np.empty((ydim, xdim, npts))
    errseries = np.empty((ydim, xdim, npts))

    # construct output light curves
    nptsx = 0
    for i in tqdm(range(ydim)):
        for j in range(xdim):
            npts = 0
            for k in range(nrows):
                if (qual[k] == 0 and np.isfinite(barytime[k])
                    and np.isfinite(fluxpixels[k, int(ydim*xdim/2)])):
                    time[npts] = barytime[k]
                    timecorr[npts] = tcorr[k]
                    cadenceno[npts] = cadno[k]
                    quality[npts] = qual[k]
                    pixseries[i, j, npts] = fluxpixels[k, nptsx]
                    errseries[i, j, npts] = errpixels[k, nptsx]
                    npts += 1
            nptsx += 1
    # define data sampling
    if filterlc:
        tpf = pyfits.open(infile)
        cadence = kepkey.cadence(tpf[1], infile, logfile, verbose)
        tr = 1.0 / (cadence / 86400)
        timescale = 1.0 / (cutoff / tr)

        # define convolution function
        if function == 'boxcar':
            filtfunc = np.ones(int(np.ceil(timescale)))
        elif function == 'gauss':
            timescale /= 2
            dx = np.ceil(timescale * 10 + 1)
            filtfunc = filtfunc([1.0, dx / 2 - 1.0, timescale],
                                np.linspace(0, dx - 1, dx))
        elif function == 'sinc':
            dx = np.ceil(timescale * 12 + 1)
            fx = np.linspace(0, dx - 1, dx)
            fx = fx - dx / 2 + 0.5
            fx /= timescale
            filtfunc = np.sinc(fx)
        filtfunc /= np.sum(filtfunc)

        # pad time series at both ends with noise model
        for i in range(ydim):
            for j in range(xdim):
                ave, sigma = (np.mean(pixseries[i, j, :len(filtfunc)]),
                              np.std(pixseries[i, j, :len(filtfunc)]))
                padded = np.append(kepstat.randarray(np.ones(len(filtfunc)) * ave,
                                   np.ones(len(filtfunc)) * sigma), pixseries[i, j, :])
                ave, sigma = (np.mean(pixseries[i, j, -len(filtfunc):]),
                              np.std(pixseries[i, j, -len(filtfunc):]))
                padded = np.append(padded,
                                   kepstat.randarray(np.ones(len(filtfunc)) * ave,
                                   np.ones(len(filtfunc)) * sigma))
                # convolve data
                convolved = np.convolve(padded, filtfunc, 'same')
                # remove padding from the output array
                outdata = convolved[len(filtfunc): -len(filtfunc)]
                # subtract low frequencies
                outmedian = np.median(outdata)
                pixseries[i, j, :] = pixseries[i, j, :] - outdata + outmedian

    # construct output file
    print("Writing output file {}...".format(outfile))
    if ydim * xdim < 1000:
        instruct = pyfits.open(infile, 'readonly')
        kepkey.history(call, instruct[0], outfile, logfile, verbose)
        hdulist = pyfits.HDUList(instruct[0])
        cols = []
        cols.append(pyfits.Column(name='TIME', format='D',
                                  unit='BJD - 2454833', disp='D12.7',
                                  array=time))
        cols.append(pyfits.Column(name='TIMECORR', format='E', unit='d',
                                  disp='E13.6', array=timecorr))
        cols.append(pyfits.Column(name='CADENCENO', format='J', disp='I10',
                                  array=cadenceno))
        cols.append(pyfits.Column(name='QUALITY', format='J', array=quality))
        for i in range(ydim):
            for j in range(xdim):
                colname = 'COL{}_ROW{}'.format(i + column, j + row)
                cols.append(pyfits.Column(name=colname, format='E',
                                          disp='E13.6',
                                          array=pixseries[i, j, :]))
        hdu1 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))
        try:
            hdu1.header['INHERIT'] = (True, 'inherit the primary header')
        except:
            pass
        try:
            hdu1.header['EXTNAME'] = ('PIXELSERIES', 'name of extension')
        except:
            pass
        try:
            hdu1.header['EXTVER' ] = (instruct[1].header['EXTVER'],
                                      'extension version number (not format version)')
        except:
            pass
        try:
            hdu1.header['TELESCOP'] = (instruct[1].header['TELESCOP'],
                                       'telescope')
        except:
            pass
        try:
            hdu1.header['INSTRUME'] = (instruct[1].header['INSTRUME'],
                                       'detector type')
        except:
            pass
        try:
            hdu1.header['OBJECT' ] = (instruct[1].header['OBJECT'],
                                      'string version of KEPLERID')
        except:
            pass
        try:
            hdu1.header['KEPLERID'] = (instruct[1].header['KEPLERID'],
                                       'unique Kepler target identifier')
        except:
            pass
        try:
            hdu1.header['RADESYS'] = (instruct[1].header['RADESYS'],
                                      'reference frame of celestial coordinates')
        except:
            pass
        try:
            hdu1.header['RA_OBJ' ] = (instruct[1].header['RA_OBJ'],
                                      '[deg] right ascension from KIC')
        except:
            pass
        try:
            hdu1.header['DEC_OBJ'] = (instruct[1].header['DEC_OBJ'],
                                      '[deg] declination from KIC')
        except:
            pass
        try:
            hdu1.header['EQUINOX'] = (instruct[1].header['EQUINOX'],
                                      'equinox of celestial coordinate system')
        except:
            pass
        try:
            hdu1.header['TIMEREF'] = (instruct[1].header['TIMEREF'],
                                      'barycentric correction applied to times')
        except:
            pass
        try:
            hdu1.header['TASSIGN'] = (instruct[1].header['TASSIGN'],
                                      'where time is assigned')
        except:
            pass
        try:
            hdu1.header['TIMESYS'] = (instruct[1].header['TIMESYS'],
                                      'time system is barycentric JD')
        except:
            pass
        try:
            hdu1.header['BJDREFI'] = (instruct[1].header['BJDREFI'],
                                      'integer part of BJD reference date')
        except:
            pass
        try:
            hdu1.header['BJDREFF'] = (instruct[1].header['BJDREFF'],
                                      'fraction of the day in BJD reference date')
        except:
            pass
        try:
            hdu1.header['TIMEUNIT'] = (instruct[1].header['TIMEUNIT'],
                                       'time unit for TIME, TSTART and TSTOP')
        except:
            pass
        try:
            hdu1.header['TSTART'] = (instruct[1].header['TSTART'],
                                     'observation start time in BJD-BJDREF')
        except:
            pass
        try:
            hdu1.header['TSTOP'] = (instruct[1].header['TSTOP'],
                                    'observation stop time in BJD-BJDREF')
        except:
            pass
        try:
            hdu1.header['LC_START'] = (instruct[1].header['LC_START'],
                                       'mid point of first cadence in MJD')
        except:
            pass
        try:
            hdu1.header['LC_END'] = (instruct[1].header['LC_END'],
                                       'mid point of last cadence in MJD')
        except:
            pass
        try:
            hdu1.header['TELAPSE'] = (instruct[1].header['TELAPSE'],
                                       '[d] TSTOP - TSTART')
        except:
            pass
        try:
            hdu1.header['LIVETIME'] = (instruct[1].header['LIVETIME'],
                                       '[d] TELAPSE multiplied by DEADC')
        except:
            pass
        try:
            hdu1.header['EXPOSURE'] = (instruct[1].header['EXPOSURE'],
                                       '[d] time on source')
        except:
            pass
        try:
            hdu1.header['DEADC'] = (instruct[1].header['DEADC'],
                                    'deadtime correction')
        except:
            pass
        try:
            hdu1.header['TIMEPIXR'] = (instruct[1].header['TIMEPIXR'],
                                       'bin time beginning=0 middle=0.5 end=1')
        except:
            pass
        try:
            hdu1.header['TIERRELA'] = (instruct[1].header['TIERRELA'],
                                       '[d] relative time error')
        except:
            pass
        try:
            hdu1.header['TIERABSO'] = (instruct[1].header['TIERABSO'],
                                       '[d] absolute time error')
        except:
            pass
        try:
            hdu1.header['INT_TIME'] = (instruct[1].header['INT_TIME'],
                                       '[s] photon accumulation time per frame')
        except:
            pass
        try:
            hdu1.header['READTIME'] = (instruct[1].header['READTIME'],
                                       '[s] readout time per frame')
        except:
            pass
        try:
            hdu1.header['FRAMETIM'] = (instruct[1].header['FRAMETIM'],
                                       '[s] frame time (INT_TIME + READTIME)')
        except:
            pass
        try:
            hdu1.header['NUM_FRM'] = (instruct[1].header['NUM_FRM'],
                                      'number of frames per time stamp')
        except:
            pass
        try:
            hdu1.header['TIMEDEL'] = (instruct[1].header['TIMEDEL'],
                                      '[d] time resolution of data')
        except:
            pass
        try:
            hdu1.header['DATE-OBS'] = (instruct[1].header['DATE-OBS'],
                                       'TSTART as UTC calendar date')
        except:
            pass
        try:
            hdu1.header['DATE-END'] = (instruct[1].header['DATE-END'],
                                       'TSTOP as UTC calendar date')
        except:
            pass
        try:
            hdu1.header['BACKAPP'] = (instruct[1].header['BACKAPP'],
                                      'background is subtracted')
        except:
            pass
        try:
            hdu1.header['DEADAPP'] = (instruct[1].header['DEADAPP'],
                                      'deadtime applied')
        except:
            pass
        try:
            hdu1.header['VIGNAPP'] = (instruct[1].header['VIGNAPP'],
                                      'vignetting or collimator correction applied')
        except:
            pass
        try:
            hdu1.header['GAIN'] = (instruct[1].header['GAIN'],
                                   '[electrons/count] channel gain')
        except:
            pass
        try:
            hdu1.header['READNOIS'] = (instruct[1].header['READNOIS'],
                                       '[electrons] read noise')
        except:
            pass
        try:
            hdu1.header['NREADOUT'] = (instruct[1].header['NREADOUT'],
                                       'number of read per cadence')
        except:
            pass
        try:
            hdu1.header['TIMSLICE'] = (instruct[1].header['TIMSLICE'],
                                       'time-slice readout sequence section')
        except:
            pass
        try:
            hdu1.header['MEANBLCK'] = (instruct[1].header['MEANBLCK'],
                                       '[count] FSW mean black level')
        except:
            pass
        hdulist.append(hdu1)
        hdulist.writeto(outfile)
        kepkey.new('EXTNAME', 'APERTURE', 'name of extension', instruct[2],
                   outfile, logfile, verbose)
        pyfits.append(outfile, instruct[2].data, instruct[2].header)
        instruct.close()
    else:
        warnmsg = ('WARNING -- KEPPIXSERIES: output FITS file requires > 999'
                   'columns. Non-compliant with FITS convention.')
        kepmsg.warn(logfile, warnmsg, verbose)

    # plot pixel array
    fmin = 1.0e33
    fmax = -1.033
    plt.figure()
    plt.clf()
    dx = 0.93 / xdim
    dy = 0.94 / ydim
    ax = plt.axes([0.06, 0.05, 0.93, 0.94])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xlim(np.min(pixcoord1) - 0.5, np.max(pixcoord1) + 0.5)
    plt.ylim(np.min(pixcoord2) - 0.5, np.max(pixcoord2) + 0.5)
    plt.xlabel('time', {'color' : 'k'})
    plt.ylabel('arbitrary flux', {'color' : 'k'})
    for i in range(ydim):
        for j in range(xdim):
            tmin = np.amin(time)
            tmax = np.amax(time)
            try:
                np.isfinite(np.amin(pixseries[i, j, :]))
                np.isfinite(np.amin(pixseries[i, j, :]))
                fmin = np.amin(pixseries[i, j, :])
                fmax = np.amax(pixseries[i, j, :])
            except:
                ugh = 1
            xmin = tmin - (tmax - tmin) / 40
            xmax = tmax + (tmax - tmin) / 40
            ymin = fmin - (fmax - fmin) / 20
            ymax = fmax + (fmax - fmin) / 20
            if kepstat.bitInBitmap(maskimg[i, j], 2):
                plt.axes([0.06 + float(j) * dx, 0.05 + i * dy, dx, dy],
                         facecolor='lightslategray')
            elif maskimg[i, j] == 0:
                plt.axes([0.06 + float(j) * dx, 0.05 + i * dy, dx, dy],
                         facecolor='black')
            else:
                plt.axes([0.06 + float(j) * dx, 0.05 + i * dy, dx, dy])
            if j == int(xdim / 2) and i == 0:
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
            elif j == 0 and i == int(ydim / 2):
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
            else:
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
            ptime = time * 1.0
            ptime = np.insert(ptime, [0], ptime[0])
            ptime = np.append(ptime, ptime[-1])
            pflux = pixseries[i, j, :] * 1.0
            pflux = np.insert(pflux, [0], -1000.0)
            pflux = np.append(pflux, -1000.0)
            plt.plot(time,pixseries[i, j, :], color='#0000ff', linestyle='-',
                     linewidth=0.5)
            if not kepstat.bitInBitmap(maskimg[i, j], 2):
                plt.fill(ptime, pflux, fc='lightslategray', linewidth=0.0,
                         alpha=1.0)
            plt.fill(ptime, pflux, fc='#FFF380', linewidth=0.0,alpha=1.0)
            if 'loc' in plottype:
                plt.xlim(xmin, xmax)
                plt.ylim(ymin, ymax)
            if 'glob' in plottype:
                plt.xlim(xmin, xmax)
                plt.ylim(1.0e-10, np.nanmax(pixseries) * 1.05)
            if 'full' in plottype:
                plt.xlim(xmin, xmax)
                plt.ylim(1.0e-10, ymax * 1.05)

    # render plot
    plt.show()
    plt.savefig(plotfile)

    # stop time
    kepmsg.clock('KEPPIXSERIES ended at', logfile, verbose)

def keppixseries_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Individual time series photometry for all pixels'
                          ' within a target mask'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-keppixseries.'),
                        default=None)
    parser.add_argument('--plotfile', default='None',
                        help='name of output PNG plot file', type=str)
    parser.add_argument('--plottype', default='global', help='Plotting type',
                        type=str, choices=['local','global','full'])
    parser.add_argument('--filterlc', action='store_true',
                        help='High-pass Filter data?')
    parser.add_argument('--function', default='boxcar', help='Type of filter',
                        type=str, choices=['boxcar','gauss','sinc'])
    parser.add_argument('--cutoff', default=1.0,
                        help='Characteristic frequency cutoff of filter [1/days]',
                        type=float)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='keppixseries.log', dest='logfile', type=str)
    args = parser.parse_args()
    keppixseries(args.infile, args.outfile, args.plotfile, args.plottype,
                 args.filterlc, args.function, args.cutoff, args.overwrite,
                 args.verbose, args.logfile)
