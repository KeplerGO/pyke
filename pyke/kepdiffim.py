from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey, kepplot, kepstat
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits


__all__ = ['kepdiffim']


def kepdiffim(infile, outfile=None, plotfile=None, imscale='logarithmic',
              colmap='PuBu', filterlc=False, function='boxcar', cutoff=1.0,
              overwrite=False, verbose=False, logfile='kepdiffim.log'):
    """
    kepdiffim -- difference imaging of pixels within a target mask

    kepdiffim plots the mean, standard deviation and chi distribution images
    for the mask contained within a target pixel file. The standard deviation
    on each pixel is defined as :math:`[flux - mean]^2 / [N - 1]`. The chi
    distribution is :math:`\sqrt{[mean - flux] ^ 2 / sigma ^ 2}`. If required,
    the data can be fed through a **boxcar**, **gaussian** or **sinc** function
    high bandpass filter in order to remove low frequency signal from the data.
    kepdiffim is a diagnostic tool for identifying source contaminants in the
    background or foreground of the target. It can be employed to identify
    pixels for inclusion or exclusion when re-extracting a Kepler light curve
    from target pixel files.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler Target
        Pixel data within the first data extension.
    outfile : str
        The name of the output FITS file. This file has four data extensions.
        The first called 'FLUX' contains an image of the pixel-by-pixel
        mean flux within target mask. The second contains the pixel variance
        image of the mask pixels over time. The third contains the standard
        deviation image, in this case the variance image is normalized to the
        median 1-:math:`\sigma` error for each pixel. The fourth extension is
        the pixel mask, as defined in the second extension of the target pixel
        file.
    plotfile : str
        Name of an optional diagnostic output plot file containing the results
        of kepdiffim. Typically this is a PNG format file. If no diagnostic
        file is required, plotfile can be **None**. If **plotfile** is **None**
        the plot will be generated but the plot will not be saved to a file.
    imscale : str
        **kepdiffim** can plot images with three choices of image scales. The
        choice is made using this argument.
        The options are:

        * linear
        * logarithmic
        * squareroot
    cmap : str
        color intensity scheme for the image display.
    filter : bool
        If **filter** is **True**, the light curve for each pixel will be
        treated by a high band-pass filter to remove long-term trends from
        e. g. differential velocity aberration.
    function : str
        The functional form of the high pass-band filter. The options are:

        * boxcar
        * gauss
        * sinc
    cutoff : float [days]
        The frequency of the high pass-band cutoff.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepdiffim kplr011390659-2010355172524_lpd-targ.fits.gz
        --filter --function boxcar --cutoff 0.1 --plotfile kepdiffim.png
        --cmap YlOrBr --imscale linear --verbose

    .. image:: ../_static/images/api/kepdiffim.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = ('KEPDIFFIM -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' plotfile={}'.format(plotfile)
            + ' imscale={}'.format(imscale)
            + ' cmap={}'.format(colmap)
            + ' filterlc={}'.format(filterlc)
            + ' function={}'.format(function)
            + ' cutoff={}'.format(cutoff)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPDIFFIM started at: ', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPDIFFIM: {} exists. Use --overwrite'.format(outfile)
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
        kepio.readTPF(infile, 'CADENCENO',logfile, verbose)
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
    print('      KepID:  %s' % kepid)
    print(' RA (J2000):  %s' % ra)
    print('Dec (J2000): %s' % dec)
    print('     KepMag:  %s' % kepmag)
    print('   SkyGroup:    %2s' % skygroup)
    print('     Season:    %2s' % str(season))
    print('    Channel:    %2s' % channel)
    print('     Module:    %2s' % module)
    print('     Output:     %1s' % output)
    print('')

    # how many quality = 0 rows?
    npts = 0
    nrows = len(fluxpixels)
    for i in range(nrows):
        if (qual[i] == 0 and np.isfinite(barytime[i])
            and np.isfinite(fluxpixels[i, int(ydim * xdim / 2)])):
            npts += 1
    time = np.empty((npts))
    timecorr = np.empty((npts))
    cadenceno = np.empty((npts))
    quality = np.empty((npts))
    pixseries = np.empty((ydim * xdim, npts))
    errseries = np.empty((ydim * xdim, npts))

    # construct output light curves
    nptsx = 0
    for i in range(ydim*xdim):
        npts = 0
        for k in range(nrows):
            if (qual[k] == 0
                and np.isfinite(barytime[k])
                and np.isfinite(fluxpixels[k, int(ydim * xdim / 2)])):
                time[npts] = barytime[k]
                timecorr[npts] = tcorr[k]
                cadenceno[npts] = cadno[k]
                quality[npts] = qual[k]
                pixseries[i, npts] = fluxpixels[k, nptsx]
                errseries[i, npts] = errpixels[k, nptsx]
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
            filtfunc = kepfunc.gauss()
            filtfunc = filtfunc([1.0, dx / 2 - 1.0, timescale],
                                linspace(0, dx - 1, dx))
        elif function == 'sinc':
            dx = np.ceil(timescale * 12 + 1)
            fx = linspace(0, dx - 1, dx)
            fx = fx - dx / 2 + 0.5
            fx /= timescale
            filtfunc = np.sinc(fx)
        filtfunc /= np.sum(filtfunc)

    # pad time series at both ends with noise model
        for i in range(ydim * xdim):
            ave, sigma  = (np.mean(pixseries[i, :len(filtfunc)]),
                           np.std(pixseries[i, :len(filtfunc)]))
            padded = np.append(kepstat.randarray(np.ones(len(filtfunc)) * ave,
                               np.ones(len(filtfunc)) * sigma), pixseries[i, :])
            ave, sigma  = (np.mean(pixseries[i, -len(filtfunc):]),
                           np.std(pixseries[i, -len(filtfunc):]))
            padded = np.append(padded,
                               kepstat.randarray(np.ones(len(filtfunc)) * ave,
                               np.ones(len(filtfunc)) * sigma))

            # convolve data
            convolved = np.convolve(padded,filtfunc,'same')
            # remove padding from the output array
            outdata = convolved[len(filtfunc):-len(filtfunc)]
            # subtract low frequencies
            outmedian = np.nanmedian(outdata)
            pixseries[i, :] = pixseries[i, :] - outdata + outmedian

    # sum pixels over cadence
    nptsx = 0
    nrows = len(fluxpixels)
    pixsum = np.zeros((ydim*xdim))
    errsum = np.zeros((ydim*xdim))
    for i in range(npts):
        if quality[i] == 0:
            pixsum += pixseries[:, i]
            errsum += errseries[:, i] **2
            nptsx += 1
    pixsum /= nptsx
    errsum = np.sqrt(errsum) / nptsx

    # calculate standard deviation pixels
    pixvar = np.zeros((ydim*xdim))
    for i in range(npts):
        if quality[i] == 0:
            pixvar += (pixsum - pixseries[:,i] / errseries[:,i])**2
    pixvar = np.sqrt(pixvar)

    # median pixel errors
    errmed = np.empty((ydim*xdim))
    for i in range(ydim*xdim):
        errmed[i] = np.median(errseries[:,i])

    # calculate chi distribution pixels
    pixdev = np.zeros((ydim*xdim))
    for i in range(npts):
        if quality[i] == 0:
            pixdev += ((pixsum - pixseries[:,i]) / pixsum)**2
    pixdev = np.sqrt(pixdev)

    # image scale and intensity limits
    pixsum_pl, zminsum, zmaxsum = kepplot.intScale1D(pixsum, imscale)
    pixvar_pl, zminvar, zmaxvar = kepplot.intScale1D(pixvar, imscale)
    pixdev_pl, zmindev, zmaxdev = kepplot.intScale1D(pixdev, imscale)

    # construct output summed image
    imgsum = np.empty((ydim, xdim))
    imgvar = np.empty((ydim, xdim))
    imgdev = np.empty((ydim, xdim))
    imgsum_pl = np.empty((ydim, xdim))
    imgvar_pl = np.empty((ydim, xdim))
    imgdev_pl = np.empty((ydim, xdim))
    n = 0
    for i in range(ydim):
        for j in range(xdim):
            imgsum[i, j] = pixsum[n]
            imgvar[i, j] = pixvar[n]
            imgdev[i, j] = pixdev[n]
            imgsum_pl[i, j] = pixsum_pl[n]
            imgvar_pl[i, j] = pixvar_pl[n]
            imgdev_pl[i, j] = pixdev_pl[n]
            n += 1

    # construct output file
    print("Writing output file {}...".format(outfile))
    instruct = pyfits.open(infile)
    kepkey.history(call, instruct[0], outfile, logfile, verbose)
    hdulist = pyfits.HDUList(instruct[0])
    hdulist.writeto(outfile)
    kepkey.new('EXTNAME', 'FLUX', 'name of extension', instruct[2], outfile,
               logfile, verbose)
    pyfits.append(outfile, imgsum, instruct[2].header)
    kepkey.new('EXTNAME', 'CHI', 'name of extension', instruct[2], outfile,
               logfile, verbose)
    pyfits.append(outfile, imgvar, instruct[2].header)
    kepkey.new('EXTNAME', 'STDDEV', 'name of extension', instruct[2], outfile,
               logfile, verbose)
    pyfits.append(outfile, imgdev, instruct[2].header)
    kepkey.new('EXTNAME', 'APERTURE', 'name of extension', instruct[2], outfile,
               logfile, verbose)
    pyfits.append(outfile, instruct[2].data,instruct[2].header)
    instruct.close()

    # pixel limits of the subimage
    ymin = row
    ymax = ymin + ydim
    xmin = column
    xmax = xmin + xdim

    # plot limits for summed image
    ymin = float(ymin) - 0.5
    ymax = float(ymax) - 0.5
    xmin = float(xmin) - 0.5
    xmax = float(xmax) - 0.5

    # plot style
    plotimage(imgsum_pl, imgvar_pl, imgdev_pl, zminsum, zminvar, zmindev,
              zmaxsum, zmaxvar, zmaxdev, xmin, xmax, ymin, ymax, colmap,
              plotfile)

    # stop time
    kepmsg.clock('KEPDIFFIM ended at: ',logfile,verbose)

# plot channel image
def plotimage(imgsum_pl, imgvar_pl, imgdev_pl, zminsum, zminvar, zmindev,
              zmaxsum, zmaxvar, zmaxdev, xmin ,xmax, ymin, ymax, colmap,
              plotfile):

    plt.figure(figsize=[15, 6])
    plt.clf()

    # plot the image window
    ax = plt.axes([0.04, 0.11, 0.31, 0.78])
    plt.imshow(imgsum_pl,aspect='auto', interpolation='nearest',
               origin='lower', vmin=zminsum, vmax=zmaxsum,
               extent=(xmin, xmax, ymin, ymax), cmap=colmap)
    plt.gca().set_autoscale_on(False)
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    plt.ylabel('Pixel Row Number', {'color' : 'k'})
    plt.title('Flux', {'color' : 'k', 'fontsize' : '16'})

    # plot the variance window
    plt.axes([0.36, 0.11, 0.31, 0.78])
    plt.imshow(imgvar_pl, aspect='auto', interpolation='nearest',
               origin='lower', vmin=zminvar, vmax=zmaxvar,
               extent=(xmin, xmax, ymin, ymax), cmap=colmap)
    plt.gca().set_autoscale_on(False)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.setp(plt.gca(),yticklabels=[])
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    try:
        plt.title(r'$\chi$ Distribution', {'color' : 'k', 'fontsize' : '16'})
    except:
        plt.title('Chi Distribution', {'color' : 'k', 'fontsize' : '16'})

    # plot the normalized standard deviation window
    plt.axes([0.68, 0.11, 0.31, 0.78])
    plt.imshow(imgdev_pl, aspect='auto', interpolation='nearest',
               origin='lower', vmin=zmindev, vmax=zmaxdev,
               extent=(xmin, xmax, ymin, ymax), cmap=colmap)
    plt.gca().set_autoscale_on(False)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.setp(plt.gca(),yticklabels=[])
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    plt.title('Normalized Standard Deviation', {'color' : 'k', 'fontsize' : '16'})

    # render plot
    plt.show()

    if plotfile is not None:
        plt.savefig(plotfile)

def kepdiffim_main():

    import argparse
    parser = argparse.ArgumentParser(
             description=('Difference imaging of pixels within a target mask'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepdiffim.'),
                        default=None)
    parser.add_argument('--plotfile', default='None',
                        help='name of output PNG plot file', type=str)
    parser.add_argument('--imscale', default='logarithmic',
                        help='type of image intensity scale', type=str,
                        choices=['linear', 'logarithmic', 'squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap', type=str)
    parser.add_argument('--filterlc', action='store_true',
                        help='High-pass Filter data?')
    parser.add_argument('--function', help='filter function', default='boxcar',
                        type=str, choices=['boxcar','gauss','sinc'])
    parser.add_argument('--cutoff',
                        help='Characteristic frequency cutoff of filter [1/days]',
                        type=float, default=1.0)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepdiffim.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepdiffim(args.infile, args.outfile, args.plotfile, args.imscale,
              args.cmap, args.filterlc, args.function, args.cutoff,
              args.overwrite, args.verbose, args.logfile)
