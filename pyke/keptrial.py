from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepfit, kepfunc, kepstat, kepfourier


__all__ = ['keptrial']


def keptrial(infile, outfile=None, datacol='SAP_FLUX', errcol='SAP_FLUX_ERR',
             fmin=0.1, fmax=50, nfreq=10, method='ft', ntrials=1000,
             plot=False, overwrite=False, verbose=False,
             logfile='keptrial.log'):
    """
    keptrial -- Calculate best period and error estimate from time series

    ``keptrial`` measures the strongest period within the frequency range
    :math:`fmin` to :math:`fmax` and estimates 1-:math:`\sigma` error
    associated with that period. The error estimate is performed by
    constructing ntrial new light curves from the original data provided in
    datacol and adjusting each individual data point according to a random
    number generator and a shot noise model. While a shot noise model is not
    uniformly applicable to all Kepler targets it provides a useful 1st order
    estimate for most. A power spectrum is calculated for each light curve
    using a user-specified method and the highest peak in each power spectrum
    recorded. The distribution of peaks is fit by a normal function, the
    centroid is adopted as the best period and 1-standard deviation error is
    taken from the standard deviation. A confidence limit is recorded as the
    range within which all trial periods fall. While this is termed a '100%'
    confidence limit, this only refers to the total number of trials rather
    than formal confidence.

    The larger the number of **ntrials**, the more robust the result. The
    values of nfreq and ntrial have to be chosen carefully to avoid excessive
    run times. The values of **fmin**, **fmax** and **nfreq** have to be
    chosen carefully in order to provide a sensible measure of period and
    error. It is recommended that ``kepft`` be used to estimate the period and
    error before attempting to use ``keptrial``. An exercise of trial and error
    will most-likely be needed to choose a permutation of :math:`fmin`,
    :math:`fmax` and :math:`nfreq` that resolves the period distribution over a
    significant number of frequency bins. If requested, the distribution and
    normal fit are plotted. The plot updates after every ntrial iteration,
    partly to relieve boredom, and partly for the user to assess whether they
    are using the correct permutation of input parameters.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing the
        results of a Monte Carlo period analysis.
    datacol : str
        The column name containing data stored within extension 1 of infile.
        This data is the input data for a series of Fourier transform
        calculations. Typically this name is SAP_FLUX (Simple Aperture
        Photometry fluxes), but any data column within extension 1 of the FITS
        file can be used provided it is coupled to an error column name using
        errcol.
    errcol : str
        The uncertainty data coupled to datacol. Typically this column is
        called SAP_FLUX_ERR.
    fmin : float [1/day]
        The minimum frequency on which each power spectrum will be calculated.
    fmax : float [1/day]
        The maximum frequency on which each power spectrum will be calculated.
    nfreq : int
        The number of uniform frequency steps between fmin and fmax over which
        the power spectrum will be calculated.
    method : str
        Choose a method for calculating the power spectrum. Currently, only
        'ft', a discrete Fourier transform, is available.
    ntrials : int
        The number of Monte Carlo trials required before calculating the best
        periods, period uncertainty and confidence in the measurement.
    plot : bool
        Plot the output window function?
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPTRIAL -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' errcol={}'.format(errcol)
            + ' fmin={}'.format(fmin)
            + ' fmax={}'.format(fmax)
            + ' nfreq={}'.format(nfreq)
            + ' method={}'.format(method)
            + ' ntrials={}'.format(ntrials)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPTRIAL started at', logfile, verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPTRIAL: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)
    # open input file
    instr = pyfits.open(infile, 'readonly')
    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # input data
    try:
        barytime = instr[1].data.field('barytime')
    except:
        barytime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                     verbose)
    signal = kepio.readfitscol(infile, instr[1].data, datacol, logfile,
                               verbose)
    err = kepio.readfitscol(infile, instr[1].data, errcol, logfile, verbose)
    # remove infinite data from time series
    try:
        nanclean = instr[1].header['NANCLEAN']
    except:
        incols = [barytime, signal, err]
        [barytime, signal, err] = kepstat.removeinfinlc(signal, incols)
    # frequency steps and Monte Carlo iterations
    deltaf = (fmax - fmin) / float(nfreq)
    freq, pmax, trial = [], [], []
    for i in tqdm(range(ntrials)):
        trial.append(i + 1)
        # adjust data within the error bars
        #work1 = kepstat.randarray(signal, err)
        # determine FT power
        fr, power = kepfourier.ft(barytime, signal, fmin, fmax, deltaf, False)
        # determine peak in FT
        pmax.append(-1.0e30)
        for j in range(len(fr)):
            if (power[j] > pmax[-1]):
                pmax[-1] = power[j]
                f1 = fr[j]
        freq.append(f1)
    # plot stop-motion histogram
    plt.figure()
    plt.clf()
    plt.axes([0.08, 0.08, 0.88, 0.89])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    n, bins, patches = plt.hist(freq, bins=nfreq, range=[fmin, fmax],
                                align='mid', rwidth=1, ec='#0000ff',
                                fc='#ffff00', lw=2)
    # fit normal distribution to histogram
    x = np.zeros(len(bins))
    for j in range(1, len(bins)):
        x[j] = (bins[j] + bins[j - 1]) / 2.
    pinit = np.array([float(i), freq[-1], deltaf])
    n = np.array(n, dtype='float32')
    coeffs, errors, covar, sigma, chi2, dof, fit, plotx, ploty = \
            kepfit.leastsquares(kepfunc.gauss, pinit, x[1:], n, None,
                                logfile, verbose)
    f = np.arange(fmin, fmax, (fmax - fmin) / 100.)
    fit = kepfunc.gauss(coeffs, f)
    plt.plot(f, fit, 'r-', linewidth=2)
    plt.xlabel(r'Frequency (1/d)', {'color' : 'k'})
    plt.ylabel('N', {'color' : 'k'})
    plt.xlim(fmin, fmax)
    plt.grid()
    # render plot
    if plot:
        plt.show()
    # period results
    p = 1.0 / coeffs[1]
    perr = p * coeffs[2] / coeffs[1]
    f1 = fmin; f2 = fmax
    gotbin = False
    for i in range(len(n)):
        if n[i] > 0 and not gotbin:
            f1 = bins[i]
            gotbin = True
    gotbin = False
    for i in range(len(n) - 1, 0, -1):
        if n[i] > 0 and not gotbin:
            f2 = bins[i + 1]
            gotbin = True
    powave, powstdev = np.mean(pmax), np.std(pmax)

    # print result
    print('              best period: %.10f days (%.7f min)' % (p, p * 1440.0))
    print('     1-sigma period error: %.10f days (%.7f min)' % (perr, perr * 1440.0))
    print('             search range: %.10f - %.10f days  ' % (1.0 / fmax, 1.0 / fmin))
    print('    100%% confidence range: %.10f - %.10f days  ' % (1.0 / f2, 1.0 / f1))
    print('         number of trials: %d' % ntrials)
    print(' number of frequency bins: %d' % nfreq)

    # history keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)

    ## write output file
    col1 = pyfits.Column(name='TRIAL', format='J', array=trial)
    col2 = pyfits.Column(name='FREQUENCY', format='E', unit='1/day',
                         array=freq)
    col3 = pyfits.Column(name='POWER', format='E', array=pmax)
    cols = pyfits.ColDefs([col1,col2,col3])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    try:
        instr[-1].header['EXTNAME'] = ('TRIALS', 'Extension name')
    except:
        raise KeyError("Could not write EXTNAME to the header of the output"
                       " file")
    try:
        instr[-1].header['SEARCHR1'] = (1.0 / fmax,
                                        'Search range lower bound (days)')
    except:
        raise KeyError("Could not write SEARCHR1 to the header of the output"
                       " file")
    try:
        instr[-1].header['SEARCHR2'] = (1.0 / fmin,
                                        'Search range upper bound (days)')
    except:
        raise KeyError("Could not write SEARCHR2 to the header of the output"
                       " file")
    try:
        instr[-1].header['NFREQ'] = (nfreq, 'Number of frequency bins')
    except:
        raise KeyError("Could not write NFREQ to the header of the output"
                       " file")
    try:
        instr[-1].header['PERIOD'] = (p, 'Best period (days)')
    except:
        raise KeyError("Could not write PERIOD to the header of the output"
                       " file")
    try:
        instr[-1].header['PERIODE'] = (perr, '1-sigma period error (days)')
    except:
        raise KeyError("Could not write PERIODE to the header of the output"
                       " file")
    try:
        instr[-1].header['CONFIDR1'] = (1.0 / f2,
                                        'Trial confidence lower bound (days)')
    except:
        raise KeyError("Could not write CONFIDR1 to the header of the output"
                       " file")
    try:
        instr[-1].header['CONFIDR2'] = (1.0 / f1,
                                        'Trial confidence upper bound (days)')
    except:
        raise KeyError("Could not write CONFIDR2 to the header of the output"
                       " file")
    try:
        instr[-1].header['NTRIALS'] = (ntrials, 'Number of trials')
    except:
        raise KeyError("Could not write NTRIALS to the header of the output"
                       " file")
    print("Writing output file {}...".format(outfile))
    instr.writeto(outfile)
    # close input file
    instr.close()
    ## end time
    kepmsg.clock('KEPTRAIL completed at', logfile, verbose)

def keptrial_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Calculate best period and error estimate from'
                          ' Fourier transform'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-keptrial.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column', type=str)
    parser.add_argument('--errcol', default='SAP_FLUX_ERR',
                        help='Name of data error column', type=str)
    parser.add_argument('--fmin', default=0.1,
                        help='Minimum search frequency [1/day]', type=float)
    parser.add_argument('--fmax', default=50.,
                        help='Maximum search frequency [1/day]', type=float)
    parser.add_argument('--nfreq', default=100,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--method', default='ft',
                        help='Frequency search method', type=str,
                        choices=['ft'])
    parser.add_argument('--ntrials', default=1000,
                        help='Number of search trials', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='keptrial.log', type=str)
    args = parser.parse_args()
    keptrial(args.infile, args.outfile, args.datacol, args.errcol, args.fmin,
             args.fmax, args.nfreq, args.method, args.ntrials, args.plot,
             args.overwrite, args.verbose, args.logfile)
