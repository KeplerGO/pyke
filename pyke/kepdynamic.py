from .utils import PyKEArgumentHelpFormatter
import re
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepstat, kepfourier


__all__ = ['kepdynamic']


def kepdynamic(infile, outfile=None, fcol='SAP_FLUX', pmin=0.1, pmax=10., nfreq=100,
               deltat=10., nslice=10, plot=False, plotscale='log', cmap='PuBu',
               overwrite=False, verbose=False, logfile='kepdynamic.log'):
    """
    kepdynamic -- Construct a dynamic (time-dependent) power spectrum from
    time series data

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : string
        The name of the output FITS file with a new image extension containing
        the dynamic power spectrum.
    fcol : str
        The name of the FITS table column in extension 1 of infile upon which
        the sequence power spectra will be calculated.
    pmin : float [day]
        The minimum of the period range over which the power spectra will be
        calculated.
    pmax : float [day]
        The maximum of the period range over which the power spectra will be
        calculated.
    nfreq : int
        The number of uniform frequency steps between :math:`1/pmax` and
        :math:`1/pmin` over which the power spectra will be calculated.
    deltat : float [days]
        The uniform length of each time slice from which a power spectrum will
        be calculated.
    nslice : int
        The number of time slices from which power spectra are calculated.
        These will be distributed uniformly across the input time series. If
        nslice is small there may be data gaps in the dynamic power spectrum.
        If nslice is large, the time slices will overlap. Both cases are valid.
    plot : boolean
        Plot the output Fourier spectrum?
    cmap : str
        A matplotlib colormap scheme.
    overwrite : bool
        Overwrite the output file?
    verbose : boolean
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepdynamic kplr002436324-2009259160929_llc.fits --fcol SAP_FLUX
        --pmin 0.08 --pmax 0.1 --nfreq 500 --deltat 5.0 --nslice 500 --plot --verbose

    .. image:: ../_static/images/api/kepdynamic.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPDYNAMIC -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' fcol={}'.format(fcol)
            + ' pmin={}'.format(pmin)
            + ' pmax={}'.format(pmax)
            + ' nfreq={}'.format(nfreq)
            + ' deltat={}'.format(deltat)
            + ' nslice={}'.format(nslice)
            + ' plot={}'.format(plot)
            + ' plotscale={}'.format(plotscale)
            + ' cmap={}'.format(cmap)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('Start time is', logfile, verbose)

    if pmin >= pmax:
        errmsg = 'ERROR -- KEPDYNAMIC: PMIN must be less than PMAX'
        kepmsg.err(logfile, errmsg, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPDYNAMIC: {} exists. Use overwrite=True'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    # read table columns
    barytime = kepio.readtimecol(infile, instr[1].data, logfile, verbose)
    signal = kepio.readfitscol(infile, instr[1].data, fcol, logfile, verbose)
    barytime = barytime + bjdref
    signal = signal / cadenom

    # remove infinite data from time series
    incols = [barytime, signal]
    outcols = kepstat.removeinfinlc(signal, incols)
    barytime = outcols[0]
    signal = outcols[1]

    # period to frequency conversion
    fmin = 1.0 / pmax
    fmax = 1.0 / pmin
    deltaf = (fmax - fmin) / nfreq

    # determine bounds of time slices
    t1, t2 = [], []
    dt = barytime[-1] - barytime[0]
    dt -= deltat
    if dt < 0:
        message = 'ERROR -- KEPDYNAMIC: time slices are larger than data range'
        kepmsg.err(logfile, message, verbose)
    ds = dt / (nslice - 1)
    for i in range(nslice):
        t1.append(barytime[0] + ds * float(i))
        t2.append(barytime[0] + deltat + ds * float(i))

    # loop through time slices
    dynam = []
    for i in tqdm(range(nslice)):
        x, y = [], []
        for j in range(len(barytime)):
            if (barytime[j] >= t1[i] and barytime[j] <= t2[i]):
                x.append(barytime[j])
                y.append(signal[j])
        x = np.array(x, dtype='float64')
        y = np.array(y, dtype='float')
        y = y - np.median(y)

        # determine FT power
        fr, power = kepfourier.ft(x, y, fmin, fmax, deltaf, False)
        for j in range(len(power)):
            dynam.append(power[j])

    # define shape of results array
    dynam = np.array(dynam, dtype='float64')
    dynam.shape = len(t1), len(power)

    # write output file
    print("Writing output file {}...".format(outfile))
    instr.append(pyfits.ImageHDU())
    instr[-1].data = dynam.transpose()
    instr[-1].header['EXTNAME'] = ('DYNAMIC FT', 'extension name')
    instr[-1].header['WCSAXES'] = (2, 'number of WCS axes')
    instr[-1].header['CRPIX1' ] = (0.5, 'reference pixel along axis 1')
    instr[-1].header['CRPIX2' ] = (0.5, 'reference pixel along axis 2')
    instr[-1].header['CRVAL1' ] = (t1[0], 'time at reference pixel (BJD)')
    instr[-1].header['CRVAL2' ] = (fmin, 'frequency at reference pixel (1/day)')
    instr[-1].header['CDELT1' ] = ((barytime[-1] - barytime[0]) / nslice,
                                   'pixel scale in dimension 1 (days)')
    instr[-1].header['CDELT2'] = (deltaf, 'pixel scale in dimension 2 (1/day)')
    instr[-1].header['CTYPE1'] = ('BJD', 'data type of dimension 1')
    instr[-1].header['CTYPE2'] = ('FREQUENCY', 'data type of dimension 2')
    instr.writeto(outfile)

    # history keyword in output file
    kepkey.history(call,instr[0],outfile,logfile,verbose)

    # close input file
    instr.close()

    # clean up x-axis unit
    time0 = float(int(barytime[0] / 100) * 100.0)
    barytime = barytime - time0
    xlab = 'BJD $-$ {}'.format(time0)

    # image intensity min and max
    if plotscale == 'log':
        dynam = np.log10(dynam)
    elif plotscale == 'squareroot':
        dynam = np.sqrt(dynam)
    elif 'loglog' in plotscale:
        dynam = np.log10(np.abs(np.log10(dynam)))
    nstat = 2; pixels = []
    for i in range(dynam.shape[0]):
        for j in range(dynam.shape[1]):
            pixels.append(dynam[i, j])
    pixels = np.array(np.sort(pixels), dtype='float')
    if int(float(len(pixels)) * 0.1 + 0.5) > nstat:
        nstat = int(float(len(pixels)) * 0.1 + 0.5)
    zmin = np.median(pixels[:nstat])
    zmax = np.median(pixels[-1:])
    if np.isnan(zmax):
        zmax = np.median(pixels[-nstat / 2:])
    if np.isnan(zmax):
        zmax = np.nanmax(pixels)

    # plot power spectrum
    plt.figure()
    plt.clf()
    plt.axes([0.08,0.113,0.91,0.86])
    dynam = dynam.transpose()
    plt.imshow(dynam, origin='lower', aspect='auto', cmap=cmap, vmin=zmin,
               vmax=zmax, extent=[barytime[0], barytime[-1], fmin, fmax],
               interpolation='bilinear')
    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
    plt.grid()
    plt.savefig(re.sub('\.\S+', '.png', outfile), dpi=100)

    plt.show()

    kepmsg.clock('KEPDYNAMIC completed at', logfile, verbose)

def kepdynamic_main():
    import argparse
    parser = argparse.ArgumentParser(
            description=('Construct a dynamic (time-dependent) power spectrum '
                         'from Kepler time series data'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepdynamic.'),
                        default=None)
    parser.add_argument('--fcol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--pmin', default=0.1,
                        help='Minimum search period [days]', type=float)
    parser.add_argument('--pmax', default=10.,
                        help='Maximum search period [days]', type=float)
    parser.add_argument('--nfreq', default=100,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--deltat', default=10.,
                        help='Length of time slice [days]', type=float)
    parser.add_argument('--nslice', default=10,
                        help='Number of time slices', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--plotscale', default='log',
                        help='type of image intensity scale', type=str,
                        choices=['linear', 'log', 'squareroot', 'loglog'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap',
                        type=str)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepdynamic.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepdynamic(args.infile, args.outfile, args.fcol, args.pmin, args.pmax,
               args.nfreq, args.deltat, args.nslice, args.plot, args.plotscale,
               args.cmap, args.overwrite, args.verbose, args.logfile)
