from . import kepio, kepmsg, kepkey, kepstat, kepfourier
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt

def kepwindow(infile, outfile, fcol='SAP_FLUX', fmax=1.0, nfreq=100, plot=False,
              overwrite=False, verbose=False, logfile='kepwindow.log'):
    """
    kepwindow -- Calculate and store the window function for a Kepler time
    series

    Kepler time stamps are not perfectly uniform. There are gaps in the data due
    to operational pauses and issues, and timestamps are corrected to the
    barycenter of the solar system. The size of the barycenter correction is
    time-dependent. kepwindow calculates a discrete window function for a
    user-provided Kepler time series. The result is stored in a new FITS file
    that is a direct copy of the input file but with an additional table
    extension containing the window function.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing the
        window function.
    fcol : str
        The name of the FITS table column in extension 1 of infile with which
        the window function should be coupled to. While the window function
        ostensibly requires the timing information, this particular piece of
        information is required so that the task can search the fcol array for
        bad data such as instances of NaN. These will be rejected before the
        window function is calculated.
    fmax : float [1/day]
        The maximum frequency on which the window function will be calculated.
    nfreq : int
        The number of uniform frequency steps between 0 and fmax over which the
        window function will be calculated.
    plot : bool
        Plot the output window function?
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepwindow kplr002436324-2009259160929_llc.fits kepwindow.fits --datacol SAP_FLUX
        --fmax 0.7 --nfreq 500 --plot --verbose

    .. image:: ../_static/images/api/kepwindow.png
        :align: center
    """

    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPWINDOW -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' fcol={}'.format(fcol)
            + ' fmax={}'.format(fmax)
            + ' nfreq={}'.format(nfreq)
            + ' plot='.format(plot)
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
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence

    ## fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    ## read table columns
    try:
        barytime = instr[1].data.field('barytime')
    except:
        barytime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                     verbose)
    signal = kepio.readfitscol(infile, instr[1].data, fcol, logfile, verbose)

    ## remove infinite data from time series
    incols = [barytime, signal]
    outcols = kepstat.removeinfinlc(signal, incols)
    barytime = outcols[0]
    signal = outcols[1]

    ## reset signal data to zero
    signal = np.ones(len(outcols[1]))
    ## frequency steps
    deltaf = fmax / nfreq
    ## loop through frequency steps; determine FT power
    fr, power = kepfourier.ft(barytime, signal, 0.0, fmax, deltaf, True)
    power[0] = 1.0

    ## mirror window function around ordinate
    work1 = []; work2 = []
    for i in range(len(fr)-1, 0, -1):
        work1.append(-fr[i])
        work2.append(power[i])
    for i in range(len(fr)):
        work1.append(fr[i])
        work2.append(power[i])
    fr = np.array(work1, dtype='float32')
    power = np.array(work2, dtype='float32')

    ## write output file
    col1 = pyfits.Column(name='FREQUENCY', format='E', unit='days', array=fr)
    col2 = pyfits.Column(name='POWER', format='E', array=power)
    cols = pyfits.ColDefs([col1, col2])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    instr[-1].header['EXTNAME'] = ('WINDOW FUNCTION', 'extension name')

    ## comment keyword in output file
    kepkey.comment(call, instr[0], outfile, logfile, verbose)
    instr.writeto(outfile)

    ## close input file
    instr.close()

    ## data limits
    nrm = len(str(int(power.max()))) - 1
    power = power / 10 ** nrm
    ylab = 'Power (x10$^%d$)' % nrm
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

    ## plot power spectrum
    if plot:
        plt.figure()
        plt.axes([0.06, 0.113, 0.93, 0.86])
        plt.plot(fr, power, color='#0000ff', linestyle='-', linewidth=1.0)
        plt.fill(fr, power, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.xlim(xmin - xr * 0.01, xmax+xr*0.01)
        if ymin - yr * 0.01 <= 0.0:
            plt.ylim(1.0e-10, ymax + yr * 0.01)
        else:
            plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
        plt.xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
        plt.ylabel('Power', {'color' : 'k'})

        plt.show()
    kepmsg.clock('KEPWINDOW completed at', logfile, verbose)

def kepwindow_main():
    import argparse
    parser = argparse.ArgumentParser(
            description=("Calculate and store the window function for a"
                         " Kepler time series"))
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of output file', type=str)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column', type=str,
                        dest='fcol')
    parser.add_argument('--fmax', default=1.0,
                        help='Minimum search frequency [1/day]',
                        type=float)
    parser.add_argument('--nfreq', default=100,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepwindow.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepwindow(args.infile, args.outfile, args.fcol, args.fmax, args.nfreq,
              args.plot, args.overwrite, args.verbose, args.logfile)
