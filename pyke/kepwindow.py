from .utils import PyKEArgumentHelpFormatter
import numpy as np
import re
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from . import kepio, kepmsg, kepkey, kepstat, kepfourier
from astropy.stats import LombScargle

__all__ = ['kepwindow']


def kepwindow(infile, outfile=None, datacol='SAP_FLUX', nyqfactor=0.01,
              plot=False, noninteractive=False, overwrite=False, verbose=False,
              logfile='kepwindow.log'):
    """
    kepwindow -- Calculate and store the window function for a Kepler time
    series

    Kepler time stamps are not perfectly uniform. There are gaps in the data due
    to operational pauses and issues, and timestamps are corrected to the
    barycenter of the solar system. The size of the barycenter correction is
    time-dependent. kepwindow calculates a discrete window function for a
    user-provided Kepler time series. This is calculated using a Lomb-Scargle
    periodogram. The result is stored in a new FITS file that is a direct copy
    of the input file but with an additional table extension containing the
    window function.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing the
        window function.
    datacol : str
        The name of the FITS table column in extension 1 of infile with which
        the window function should be coupled to. While the window function
        ostensibly requires the timing information, this particular piece of
        information is required so that the task can search the datacol array for
        bad data such as instances of NaN. These will be rejected before the
        window function is calculated.
    nyqfactor : int
        The number of nyquist factors up to which to evaluate. Kepler data is
        fairly regular and so the default will usually encompass most of the window.
    plot : bool
        Plot the output window function?
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

        $ kepwindow kplr002436324-2009259160929_llc.fits --datacol SAP_FLUX
        --nyqfactor 0.02 --plot --verbose

    .. image:: ../_static/images/api/kepwindow.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPWINDOW -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' nyqfactor={}'.format(nyqfactor)
            + ' plot='.format(plot)
            + ' noninteractive='.format(noninteractive)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    ## start time
    kepmsg.clock('KEPWINDOW started at', logfile, verbose)
    ## overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPWINDOW: {} exists. Use overwrite=True'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    ## open input file
    instr = pyfits.open(infile)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    try:
        barytime = instr[1].data.field('barytime')
    except:
        barytime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                     verbose)
    barytime = barytime[np.isfinite(barytime)]
    ls = LombScargle(barytime, 1, center_data=False, fit_mean=False)
    freqW, powerW = ls.autopower(nyquist_factor=nyqfactor)
    freqW = np.append(-freqW[::-1], freqW)
    powerW = np.append(powerW[::-1], powerW)
    freqW, powerW = np.sort(freqW), powerW[np.argsort(freqW)]
    freqW, powerW = freqW[np.isfinite(powerW)], powerW[np.isfinite(powerW)]

    plt.figure()
    plt.axes([0.06, 0.113, 0.93, 0.86])
    plt.plot(freqW, powerW, color='#363636', linestyle='-', linewidth=1.0)
    plt.fill(freqW, powerW, color='#a8a7a7', linewidth=0.0, alpha=0.2)
    plt.xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
    plt.ylabel('Power', {'color' : 'k'})
    plt.grid()
    plt.savefig(re.sub('.fits', '.png', outfile), bbox_inches='tight')
    if not noninteractive:
        plt.show()

    col1 = pyfits.Column(name='FREQUENCY', format='E', unit='days', array=freqW)
    col2 = pyfits.Column(name='POWER', format='E', array=powerW)
    cols = pyfits.ColDefs([col1, col2])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    instr[-1].header['EXTNAME'] = ('WINDOW FUNCTION', 'extension name')
    kepmsg.log(logfile, "Writing output file {}...".format(outfile), verbose)

    ## comment keyword in output file
    kepkey.comment(call, instr[0], outfile, logfile, verbose)
    instr.writeto(outfile)

    ## close input file
    instr.close()
    kepmsg.clock('KEPWINDOW completed at', logfile, verbose)
    return

def kepwindow_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=("Calculate and store the window function for a"
                          " Kepler time series"),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile',
                        help=('The name of a MAST standard format FITS'
                              ' file containing a Kepler light curve'),
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepwindow.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column', type=str,
                        dest='datacol')
    parser.add_argument('--nyqfactor', default=0.01,
                        help='number of nyquist frequencies to evaluate up to', type=float)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--non-interactive', action='store_true',
                        help='Pop up matplotlib plot window?',
                        dest='noninteractive')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepwindow.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepwindow(args.infile, args.outfile, args.datacol, args.nyqfactor,
              args.plot, args.noninteractive, args.overwrite, args.verbose,
              args.logfile)
