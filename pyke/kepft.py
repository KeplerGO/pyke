from .utils import PyKEArgumentHelpFormatter
import math
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from . import kepio, kepmsg, kepkey, kepstat, kepfourier


__all__ = ['kepft']


def kepft(infile, outfile=None, fcol='SAP_FLUX', pmin=0.1, pmax=10., nfreq=100,
          plot=False, overwrite=False, verbose=False, logfile='kepft.log'):
    """
    kepft -- Calculate and store a Fourier Transform from a Kepler time series

    ``kepft`` calculates the discrete Fourier transform for a user-provided
    Kepler time series. The result is stored in a new FITS file that is a
    direct copy of the input file but with an additional table extension
    containing the power spectrum.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing the
        Fourier spectrum.
    fcol : str
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
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepft kplr002436324-2009259160929_llc.fits --pmin 0.5
          --pmax 100 --nfreq 1000 --plot --verbose

    .. image:: ../_static/images/api/kepft.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPFT -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' fcol={}'.format(fcol)
            + ' pmin={}'.format(pmin)
            + ' pmax={}'.format(pmax)
            + ' nfreq={}'.format(nfreq)
            + ' plot={}'.format(plot)
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
        errmsg = 'ERROR -- KEPFT: {} exists. Use --overwrite'.format(outfile)
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
    signal = kepio.readfitscol(infile, instr[1].data, fcol, logfile, verbose)
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
    fr, power = kepfourier.ft(barytime, signal, fmin, fmax, deltaf, True)
    ## write output file
    col1 = pyfits.Column(name='FREQUENCY', format='E', unit='1/day',
                         array=fr)
    col2 = pyfits.Column(name='POWER', format='E', array=power)
    cols = pyfits.ColDefs([col1, col2])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    instr[-1].header['EXTNAME'] = ('POWER SPECTRUM', 'extension name')
    print("Writing output file {}...".format(outfile))
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
        plt.show()
    ## end time
    kepmsg.clock('KEPFT completed at', logfile, verbose)

def kepft_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Calculate and store a Fourier Transform from a'
                          ' Kepler time series.'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepft.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--pmin', default=0.1,
                        help='Minimum search period [days]', type=float)
    parser.add_argument('--pmax', default=10.,
                        help='Maximum search period [days]', type=float)
    parser.add_argument('--nfreq', default=100,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepft.log', type=str)
    args = parser.parse_args()
    kepft(args.infile, args.outfile, args.datacol, args.pmin, args.pmax,
          args.nfreq, args.plot, args.overwrite, args.verbose, args.logfile)
