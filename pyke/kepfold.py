from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey, kepstat, kepfit
import numpy as np
from scipy import stats
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
import re

__all__ = ['kepfold']


def kepfold(infile, outfile=None, period=None, bjd0=None, bindata=False,
            nbins=1000, datacol='SAP_FLUX',noninteractive=False,
            overwrite=False, verbose=False, logfile="kepfold.log"):
    """
    kepfold: Phase-fold light curve data on linear ephemeris.

    kepfold calculates the phase of all time-tagged data points relative to a
    user-supplied linear ephemeris. The relation is:


    :math:`TIME` is the column within the FITS light curve file containing
    barycenter-corrected time stamps. :math:`bjd0` is a user-supplied BJD for
    zero phase. period is a user-supplied period in units of days. PHASE is the
    calculated phase for each time stamp; these results are written to a new
    float column in the LIGHT CURVE extension of the input file before being
    exported as a new file with name defined by the user. Optionally, kepfold
    will plot the data folded on the ephemeris and store it within a new FITS
    extension of the output file called FOLDED.

    Parameters
    ----------
    inile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file with a new extension containing a
        phased light curve.
    period : str
        Period over which to fold the light curve, in units of days.
    bjd0 : float
        Time of zero phase for the folded data, in units of BJD.
    bindata: bool
        Collect the data into discrete bins during the fold?
    nbins : int
        The number of phase bins to calculate.
    datacol : str
        The name of the FITS table column in extension 1 of infile upon which
        the Fourier transform will be calculated.
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

        $ kepfold kplr010544976-2009201121230_slc.fits
          0.350471 2455002.825 --bindata --verbose
    .. image:: ../_static/images/api/kepfold.png
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # check if there is period or BJD0 information
    if np.any([(period is None), (bjd0 is None)]):
        # open input file
        kepmsg.log(logfile,'KEPFOLD -- Searching for periods in headers.',verbose)
        instr = pyfits.open(infile, 'readonly')
        if period is None:
            for i in instr:
                if 'PERIOD' in i.header:
                    period = i.header['PERIOD']
        if bjd0 is None:
            for i in instr:
                if 'BJD0' in i.header:
                    bjd0 = i.header['BJD0']

        if np.all([(period is None), (bjd0 is None)]):
            errmsg = ('ERROR -- KEPFOLD: No period information found. Either specify'
            'period and BJD0 or run KEPBLS.')
            kepmsg.err(logfile, errmsg, verbose)

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPFOLD -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' period={}'.format(period)
            + ' bjd0={}'.format(bjd0)
            + ' bindata={}'.format(bindata)
            + ' datacol={}'.format(datacol)
            + ' noninteractive='.format(noninteractive)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPFOLD started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPFOLD: {} exists. Use --overwrite'
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

    try:
        barytime = instr[-1].data.field('barytime') + bjdref
    except:
        barytime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                     verbose) + bjdref

    signal = kepio.readfitscol(infile, instr[1].data, datacol, logfile, verbose)
    try:
        err = kepio.readfitscol(infile, instr[1].data, '{}_ERR'.format(datacol), logfile, verbose)
    except:
        err = np.copy(signal) / 0 #NaNs
    barytime , signal , err = kepstat.removeinfinlc(signal, [barytime, signal, err])
    if bjd0 is None:
        bjd0 = 0.

    phase = (barytime - bjd0) /period % 1

    if bindata:
        bs = (np.linspace(0,1,nbins))
        db = np.median(bs[1:]-bs[0:-1])
        binned = np.zeros(len(bs))
        binned_err = np.zeros(len(bs))
        for i,b in enumerate(bs):
            p = np.where( (phase > (b - db/2.)) & (phase <= (b + db/2.)) )[0]
            if len(p) == 0:
                binned[i]=0
                binned_err[i]==0
                continue
            else:
                binned[i] = np.mean(signal[p])
                binned_err[i] = (np.sum(err[p]**2)**0.5)/float(len(p))

    # update HDU1 for output file
    col0 = pyfits.Column(name='PHASE', format='E',
                                            array=phase)
    col1 = pyfits.Column(name='FLUX', format='E',
                                            unit='e/s', array=signal)
    col2 = pyfits.Column(name='ERR'.format(datacol), format='E',
                                            unit='e/s', array=err)
    cols = pyfits.ColDefs([col0, col1, col2])
    instr.append(pyfits.BinTableHDU.from_columns(cols))

    instr[-1].header.cards['TTYPE1'].comment = 'column title: phase'
    instr[-1].header.cards['TTYPE2'].comment = 'column title: {}'.format(datacol)
    instr[-1].header.cards['TTYPE3'].comment = 'column title: {} 1-sigma error'.format(datacol)
    instr[-1].header.cards['TFORM1'].comment = 'column type: float32'
    instr[-1].header.cards['TFORM2'].comment = 'column type: float32'
    instr[-1].header.cards['TFORM3'].comment = 'column type: float32'
    instr[-1].header.cards['TUNIT2'].comment = 'column units: electrons per second'
    instr[-1].header.cards['TUNIT3'].comment = 'column units: electrons per second'
    instr[-1].header['PERIOD'] = (period, 'period defining the phase [d]')
    instr[-1].header['BJD0'] = (bjd0, 'time of phase zero [BJD]')
    instr[-1].header['EXTNAME'] = ('FOLDED LC', 'extension name')


    if bindata:
        col3 = pyfits.Column(name='PHASE_BINNED', format='E',
                                                array=bs)
        col4 = pyfits.Column(name='FLUX_BINNED', format='E',
                                                unit='e/s', array=binned)
        col5 = pyfits.Column(name='ERR_BINNED'.format(datacol), format='E',
                                                unit='e/s', array=binned_err)
        cols = pyfits.ColDefs([col3, col4, col5])
        instr.append(pyfits.BinTableHDU.from_columns(cols))
        instr[-1].header.cards['TTYPE1'].comment = 'column title: phase'
        instr[-1].header.cards['TTYPE2'].comment = 'column title: {}'.format(datacol)
        instr[-1].header.cards['TTYPE3'].comment = 'column title: {} 1-sigma error'.format(datacol)
        instr[-1].header.cards['TFORM1'].comment = 'column type: float32'
        instr[-1].header.cards['TFORM2'].comment = 'column type: float32'
        instr[-1].header.cards['TFORM3'].comment = 'column type: float32'
        instr[-1].header.cards['TUNIT2'].comment = 'column units: electrons per second'
        instr[-1].header.cards['TUNIT3'].comment = 'column units: electrons per second'
        instr[-1].header['NBINS'] = (nbins, 'number of bins')
        instr[-1].header['PERIOD'] = (period, 'period defining the phase [d]')
        instr[-1].header['BJD0'] = (bjd0, 'time of phase zero [BJD]')
        instr[-1].header['EXTNAME'] = ('BINNED FOLDED LC', 'extension name')


    # history keyword in output file
    kepmsg.log(logfile, "Writing output file {}...".format(outfile), verbose)
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    instr.writeto(outfile)

    # close input file
    instr.close()

    # plot new light curve
    plt.figure()
    plt.clf()
    xlab = 'Orbital Phase ($\phi$)'
    ylab = 'e$^-$ s$^{-1}$'
    ax = plt.axes([0.06, 0.11, 0.93, 0.86])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    labels = ax.get_yticklabels()
    plt.fill(phase, signal, color='#a8a7a7', linewidth=0.0, alpha=0.2)
    if bindata:
        plt.errorbar(bs, binned, binned_err, color='#363636', linestyle='', linewidth=2.0,
                marker='.')
    else:
        plt.errorbar(phase, signal, err, color='#363636', linestyle='', linewidth=2.0,
                marker='.')

    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab, {'color' : 'k'})
    plt.grid()
    plt.savefig(re.sub('.fits', '.png', outfile),bbox_inches='tight')
    if not noninteractive:
        plt.show()

    # stop time
    kepmsg.clock('KEPFOLD ended at: ', logfile, verbose)
    return

def kepfold_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=("Phase-fold light curve data on linear ephemeris."),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of FITS input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepfold.'),
                        default=None)
    parser.add_argument('--period', help='Period to fold data upon [days]',
                        type=float, default=None)
    parser.add_argument('--bjd0',
                        help='time of zero phase for the folded period [BJD]',
                        type=float, default=None)
    parser.add_argument('--bindata', action='store_true',
                        help='Bin output data?')
    parser.add_argument('--nbins', default=1000, help='Number of period bins',
                        type=int)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--non-interactive', action='store_true',
                        help='Pop up matplotlib plot window?',
                        dest='noninteractive')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepfold.log', dest='logfile', type=str)
    args = parser.parse_args()

    kepfold(args.infile, args.outfile, args.period, args.bjd0, args.bindata,
            args.nbins, args.datacol, args.noninteractive, args.overwrite,
            args.verbose, args.logfile)
