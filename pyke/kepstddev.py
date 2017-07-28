from . import kepio
from . import kepmsg
from . import kepstat
from . import kepkey
from .utils import PyKEArgumentHelpFormatter
import math
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from scipy import stats
from copy import copy
from tqdm import tqdm


__all__ = ['kepstddev']


def kepstddev(infile, outfile=None, datacol='PDCSAP_FLUX', timescale=6.5,
              overwrite=False, verbose=False, logfile='kepstddev.log'):
    """
    kepstddev -- Calculate Combined Differential Photometric Precision for a
    time series light curve.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler light
        curve data within the first data extension. While the kepstddev task
        will perform a calculation upon any data column in the first FITS
        extension of the input file, the output is only meaningful if the data
        column has first been normalized to a time-dependent model or function.
        For example, the kepflatten task removes astrophysical and systematic
        features in a light curve by fitting and normalizing to a running
        polynomial. The data output of kepflatten is a column named DETSAP_FLUX
        and is a suitable input column for kepstddev.
    outfile : str
        Name of the output FITS file containing the results of kepstddev. The
        output file is a direct copy of infile except for the addition of a new
        column called CDPPnn, where nn refers to the CDPP timescale. If the
        column CDPPnn exists already in the input file, then it will be
        overwritten in the output file.
    datacol : str
        The name of the FITS data column upon which to calculate CDPP. datacol
        must be a column within the FITS table extension of the light curve -
        the first extension of the input file. The time-series within datacol
        must be normalized, by e.g. the task kepflatten.
    timescale : float
        The characteristic timescale over which to calculate CDPP. The units
        are hours.
    overwrite : bool
        Overwrite the output file? if overwrite = False and an existing file has the
        same name as outfile then the task will stop with an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    After applying ``kepflatten`` to ``kplr002437145-2009350155506_llc.fits``,
    one can input it to ``kepstddev``:

    .. code-block:: bash

        $ kepstddev kplr002437145-2009350155506_llc.fits --datacol DETSAP_FLUX
        --verbose
        --------------------------------------------------------------
        KEPSTDDEV --  infile=kplr002437145-2009350155506_llc.fits
        outfile=kplr002437145-2009350155506_llc-kepstddev.fits datacol=DETSAP_FLUX
        timescale=6.5 overwrite=False verbose=True logfile=kepstddev.log

        Standard deviation = 1295.0731328136349 ppm
        Median 6.5hr CDPP = 313 ppm
           RMS 6.5hr CDPP = 329 ppm
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPSTDDEV -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' timescale={}'.format(timescale)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPSTDDEV started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPSTDDEV: {} exists. Use overwrite=True'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)

    # filter input data table
    work1 = np.array([table.field('time'), table.field(datacol)])
    work1 = np.rot90(work1, 3)
    work1 = work1[~np.isnan(work1).any(1)]

    # read table columns
    intime = work1[:, 1] + bjdref
    indata = work1[:, 0]

    # calculate STDDEV in units of ppm
    stddev = kepstat.running_frac_std(intime, indata, timescale / 24) * 1.0e6
    astddev = np.std(indata) * 1.0e6
    cdpp = stddev / math.sqrt(timescale * 3600.0 / cadence)
    print('\nStandard deviation = {} ppm'.format(astddev))

    # filter cdpp
    for i in range(len(cdpp)):
        if cdpp[i] > np.median(cdpp) * 10.0:
            cdpp[i] = cdpp[i - 1]

    # calculate median STDDEV
    med = np.median(cdpp[:])
    medcdpp = np.ones(len(cdpp), dtype='float32') * med
    print('Median %.1fhr CDPP = %d ppm' % (timescale, med))

    # calculate RMS STDDEV
    rms = kepstat.rms(cdpp, np.zeros(len(stddev)), logfile, verbose)
    rmscdpp = np.ones((len(cdpp)), dtype='float32') * rms
    print('   RMS %.1fhr CDPP = %d ppm\n' % (timescale, rms))

    # clean up x-axis unit
    intime0 = float(int(tstart / 100) * 100.0)
    ptime = intime - intime0
    xlab = 'BJD $-$ {}'.format(intime0)

    # clean up y-axis units
    pout = copy(cdpp)
    nrm = math.ceil(math.log10(np.median(cdpp))) - 1.0
    ylab = '%.1fhr $\sigma$ (ppm)' % timescale

    # data limits
    xmin = ptime.min()
    xmax = ptime.max()
    ymin = pout.min()
    ymax = pout.max()
    xr = xmax - xmin
    yr = ymax - ymin
    ptime = np.insert(ptime,[0],[ptime[0]])
    ptime = np.append(ptime,[ptime[-1]])
    pout = np.insert(pout,[0],[0.0])
    pout = np.append(pout,0.0)

    # define size of plot on monitor screen
    plt.figure()
    # delete any fossil plots in the matplotlib window
    plt.clf()
    # position first axes inside the plotting window
    ax = plt.axes([0.07, 0.15, 0.92, 0.83])
    # force tick labels to be absolute rather than relative
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    # plot flux vs time
    ltime = np.array([], dtype='float64')
    ldata = np.array([], dtype='float32')
    dt = 0
    work1 = 2.0 * cadence / 86400
    for i in range(1, len(ptime) - 1):
        dt = ptime[i] - ptime[i - 1]
        if dt < work1:
            ltime = np.append(ltime, ptime[i])
            ldata = np.append(ldata, pout[i])
        else:
            plt.plot(ltime, ldata, color='#0000ff', linestyle='-',
                     linewidth=1.0)
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
    plt.plot(ltime, ldata, color='#0000ff', linestyle='-', linewidth=1.0)

    # plot the fill color below data time series, with no data gaps
    plt.fill(ptime, pout, fc='#ffff00', linewidth=0.0, alpha=0.2)

    # define plot x and y limits
    plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
    if ymin - yr * 0.01 <= 0.0:
        plt.ylim(1.0e-10, ymax + yr * 0.01)
    else:
        plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
    # plot labels
    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab, {'color' : 'k'})
    # make grid on plot
    plt.grid()

    # render plot
    plt.show()

    # add NaNs back into data
    n = 0
    work1 = np.array([], dtype='float32')
    instr = pyfits.open(infile)
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    for i in tqdm(range(len(table.field(0)))):
        if np.isfinite(table.field('time')[i]) and np.isfinite(table.field(datacol)[i]):
            work1 = np.append(work1, cdpp[n])
            n += 1
        else:
            work1 = np.append(work1, np.nan)

    # write output file
    print("Writing output file {}...".format(outfile))
    kepkey.new('MCDPP%d' % (timescale * 10.0), medcdpp[0],
               'Median %.1fhr CDPP (ppm)' % timescale,
               instr[1], outfile, logfile, verbose)
    kepkey.new('RCDPP%d' % (timescale * 10.0), rmscdpp[0],
               'RMS %.1fhr CDPP (ppm)' % timescale,
               instr[1], outfile, logfile, verbose)
    colname = 'CDPP_{}'.format(timescale * 10)
    col1 = pyfits.Column(name=colname, format='E13.7', array=work1)
    cols = instr[1].data.columns + col1
    instr[1] = pyfits.BinTableHDU.from_columns(cols, header=instr[1].header)
    instr.writeto(outfile)
    # comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    # close FITS
    instr.close()

    # end time
    kepmsg.clock('KEPSTDDEV completed at', logfile, verbose)

def kepstddev_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Calculate Combined Differential Photometric'
                          ' Precision for a time series light curve'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input FITS file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepstddev.'),
                        default=None)
    parser.add_argument('--datacol', default='PDCSAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--timescale', '-t', default=6.5,
                        help='CDPP timescale', type=float)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepstddev.log', type=str)
    args = parser.parse_args()
    kepstddev(args.infile, args.outfile, args.datacol, args.timescale,
              args.overwrite,args.verbose,
              args.logfile)
