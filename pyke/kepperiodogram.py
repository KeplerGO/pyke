from .utils import PyKEArgumentHelpFormatter
import math
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from . import kepio, kepmsg, kepkey, kepstat
from astropy.stats import LombScargle

__all__ = ['kepperiodogram']


def kepperiodogram(infile, outfile=None, datacol='PDCSAP_FLUX', pmin=0.1, pmax=10., nfreq=2000,
          plot=False, noninteractive=False, overwrite=False, verbose=False,
          logfile='kepperiodogram.log'):
    """
    kepperiodogram -- Calculate and store a Lomb Scargle Periodogram based on a
    Kepler time series. The result is stored in a new FITS file that is a
    direct copy of the input file but with an additional table extension
    containing the periodogram.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing the
        Fourier spectrum.
    datacol : str
        The name of the FITS table column in extension 1 of infile upon which
        the Fourier transform will be calculated.
    pmin : float [day]
        The minimum of the period range over which the Fourier transform will
        be calculated.
    pmax : float [day]
        The maximum of the period range over which the Fourier transform will
        be calculated.
    nfreq : int
        The number of uniform frequency steps between :math:`1/pmax` and
        :math:`1/pmin` that the Fourier transform will be calculated.
    plot : bool
        Plot the output Fourier spectrum?
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

        $ kepperiodogram kplr002436324-2009259160929_llc.fits --pmin 0.5
          --pmax 100 --nfreq 1000 --plot --verbose

    .. image:: ../_static/images/api/kepperiodogram.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('kepperiodogram -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' pmin={}'.format(pmin)
            + ' pmax={}'.format(pmax)
            + ' nfreq={}'.format(nfreq)
            + ' plot={}'.format(plot)
            + ' noninteractive={}'.format(noninteractive)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    ## start time
    kepmsg.clock('Start time is', logfile, verbose)
    ## overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- kepperiodogram: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)
    ## open input file
    instr = pyfits.open(infile)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    ## fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    ## read table columns
    try:
        barytime = instr[1].data.field('barytime') + bjdref
    except:
        barytime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                     verbose) + bjdref
    signal = kepio.readfitscol(infile, instr[1].data, datacol, logfile, verbose)
    ## remove infinite data from time series
    incols = [barytime, signal]
    outcols = kepstat.removeinfinlc(signal, incols)
    barytime = outcols[0]
    signal = outcols[1] - np.median(outcols[1])
    ## period to frequency conversion
    fmin = 1.0 / pmax
    fmax = 1.0 / pmin
    deltaf = (fmax - fmin) / nfreq
    ## loop through frequency steps; determine FT power
    fr = np.linspace(fmin,fmax,nfreq)
    power = LombScargle(barytime, signal, deltaf).power(fr)
    #find highest power period
    period = 1. / fr[power.argmax()]

    ## write output file
    col1 = pyfits.Column(name='FREQUENCY', format='E', unit='1/day',
                         array=fr)
    col2 = pyfits.Column(name='POWER', format='E', array=power)
    cols = pyfits.ColDefs([col1, col2])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    instr[-1].header['EXTNAME'] = ('POWER SPECTRUM', 'extension name')
    instr[-1].header['PERIOD'] = (period, 'most significant trial period [d]')

    kepmsg.log(logfile, "kepperiodogram - best period found: {}".format(period), verbose)
    kepmsg.log(logfile, "Writing output file {}...".format(outfile), verbose)
    instr.writeto(outfile)
    ## history keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    ## close input file
    instr.close()

    ## data limits
    nrm = int(math.log10(power.max()))
    power = power / 10 ** nrm
    ylab = 'Power (x10$^{}$)'.format(nrm)
    xmin = fr.min()
    xmax = fr.max()
    ymin = power.min()
    ymax = power.max()
    xr = xmax - xmin
    yr = ymax - ymin
    fr = np.insert(fr, [0], fr[0])
    fr = np.append(fr, fr[-1])
    power = np.insert(power, [0], 0.0)
    power = np.append(power, 0.0)

    if plot:
        plt.figure()
        plt.clf()
        plt.axes([0.06, 0.113, 0.93, 0.86])
        plt.plot(fr, power, color='#0000ff', linestyle='-', linewidth=1.0)
        plt.fill(fr, power, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin - yr * 0.01 <= 0.0:
            plt.ylim(1.0e-10, ymax + yr * 0.01)
        else:
            plt.ylim(ymin - yr * 0.01, ymax + yr *0.01)
        plt.xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()
        # render plot
        if not noninteractive:
            plt.show()
    ## end time
    kepmsg.clock('kepperiodogram completed at', logfile, verbose)

def kepperiodogram_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Calculate and store a Fourier Transform from a'
                          ' Kepler time series.'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepperiodogram.'),
                        default=None)
    parser.add_argument('--datacol', default='PDCSAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--pmin', default=0.1,
                        help='Minimum search period [days]', type=float)
    parser.add_argument('--pmax', default=10.,
                        help='Maximum search period [days]', type=float)
    parser.add_argument('--nfreq', default=2000,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--non-interactive', action='store_true',
                        help='Pop up matplotlib plot window?',
                        dest='noninteractive')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepperiodogram.log', type=str)
    args = parser.parse_args()

    kepperiodogram(args.infile, args.outfile, args.datacol, args.pmin, args.pmax,
          args.nfreq, args.plot, args.noninteractive, args.overwrite,
          args.verbose, args.logfile)
