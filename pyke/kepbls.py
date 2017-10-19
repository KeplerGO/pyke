from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey
import math
import numpy as np
from copy import copy
from scipy import stats
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm


__all__ = ['kepbls']


def kepbls(infile, outfile=None, datacol='DETSAP_FLUX',
           errcol='DETSAP_FLUX_ERR', minper=1.0, maxper=30, mindur=0.5,
           maxdur=12, nsearch=1000, nbins=1000, plot=False, overwrite=False,
           verbose=False, logfile='kepbls.log'):
    """
    kepbls -- Perform Box-Least Square searches for periodic exoplanet transits

    Parameters
    ----------
    infile : str
        The name of a standard format FITS file containing a Kepler light
        curve within the first data extension. The data in **infile** will
        typically have been flattened by ``kepflatten``. Multiple quarters can
        be searched by appending light curves within a single file using
        ``kepstitch``.
    outfile : str
        The name of the output FITS file. **outfile** will be a direct copy of
        **infile** but with a new extension called BLS appended containing a
        table of i) trial periods, PERIOD, ii) a reference Barycentric Julian
        Date (BJD) corresponding the center of the most transit-like structure
        in the folded light curve at the trial period, BJD0, iii) a duration
        (in hours) corresponding the width of the most transit-like structure
        in the folded light curve at the trial period, DURATION, and iv) the
        normalized signal residue of the most transit-like structure in the
        folded light curve at the trial period, SIG_RES. The definition of
        SIG_RES is provided in equation 5 of Kovacs, Zucker and Mazeh (2002).
        The maximum calculated value of SIG_RES and the corresponding trial
        period, BJD epoch and transit duration are stored as keywords in the
        BLS extension called SIGNRES, PERIOD, BJD0, TRANSDUR.
    datacol : str
        The column name containing data stored within FITS extension 1 of
        **infile**. This data will be searched for outliers. Typically this
        name is DETSAP_FLUX (Detrended Simple Aperture Photometry fluxes). This
        version of the data is computed by the task ``kepflatten``.
        Other flux data will be accepted - SAP_FLUX (Simple Aperture
        Photometry), PDCSAP_FLUX (Pre-search Data Conditioning fluxes) or
        CBVSAP_FLUX (SAP_FLUX corrected for systematic artifacts by the PyKE
        tool kepcotrend). However neither of these three options are
        recommended because the flux data contain either astrophysical
        variability, systematic variability, or both.
    errcol : str
        The column name containing photometric 1-sigma errors stored within
        extension 1 of infile. Typically this name is DETSAP_FLUX_ERR.
    minper : float [days]
        The shortest trial period on which to search for transits.
    maxper : float [days]
        The longest trial period on which to search for transits.
    mindur : float [hours]
        For each trial period, the BLS function will be fit to the data by
        i) iterating upon the epoch of mid-transit in the model, and
        ii) adjusting the width of the modeled transit. The width is adjusted
        systematically in step sizes equaling the cadence of the input data.
        **mindur** provides a lower limit to the range of transit widths
        tested.
    maxdur : float [hours]
        Provides an upper limit to the range of transit widths tested over each
        trial period.
    nsearch : int
        The number of trial periods to search between the lower bound
        **minper** and the upper bound **maxper**.
    nbins : int
        Before the BLS transit model is fit to the data, data are folded upon
        the trail orbital period and then phase binned by calculating the mean
        flux level within each bin interval. **nbins** is the number of phase bins
        in which to store the data before each fit.
    plot : bool
        Plot the calculated Normalized Signal Residue as a function of trial
        orbital period?
    overwrite : bool
        Overwrite the output file? If overwrite is False and an existing file
        has the same name as outfile then the task will stop with an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------

    After using ``kepflatten`` to remove low frequency variability in ``kplr011904151-2009350155506_llc.fits``,
    we can use the output, ``kepflatten.fits``, into ``kepbls``, i.e.,

    .. code-block:: bash

        $ kepbls kepflatten.fits --datacol DETSAP_FLUX --errcol DETSAP_FLUX_ERR
        --minper 0.8 --maxper 1.0 --mindur 1.0 --maxdur 12.0 --nsearch 1000
        --nbins 1000 --plot --verbose

             Best trial period = 0.8375062346458435 days
           Time of mid-transit = BJD 2455093.5457274746
              Transit duration = 1.7099086232483387 hours
        Maximum signal residue = 4.487271046981536e-06

    .. image:: ../_static/images/api/kepbls.png
        :align: center
    """
    # log the call
    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPBLS -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' errcol={}'.format(errcol)
            + ' minper={}'.format(minper)
            + ' maxper={}'.format(maxper)
            + ' mindur={}'.format(mindur)
            + ' maxdur={}'.format(maxdur)
            + ' nsearch={}'.format(nsearch)
            + ' nbins={}'.format(nbins)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPBLS started at', logfile, verbose)

    # is duration greater than one bin in the phased light curve?
    if not nbins * maxdur / 24.0 / maxper > 1.0:
        message = ('WARNING -- KEPBLS: {}'
                   ' hours transit duration < 1 phase bin when P = {} days.'
                   .format(maxdur, maxper))
        kepmsg.warn(logfile, message, verbose)

    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPBLS: {} exists. Use overwrite=True'.format(outfile)
        kepmsg.err(logfile, message, verbose)

    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # read table columns
    intime = np.array(table.field('time')) + bjdref
    indata = np.array(table.field(datacol))
    inerr = np.array(table.field(errcol))
    # filter input data table
    good_data_mask = (intime == intime) & (indata == indata)
    indata = indata[good_data_mask]
    intime = intime[good_data_mask]
    inerr = inerr[good_data_mask]

    # test whether the period range is sensible
    tr = intime[-1] - intime[0]
    if maxper > tr:
        message = ('ERROR -- KEPBLS: maxper is larger than the time range of'
                   ' the input data')
        kepmsg.err(logfile, message, verbose)

    # prepare time series
    time_arr = intime - intime[0]
    flux_arr = indata - np.nanmean(indata)

    # start period search
    srMax = np.array([], dtype='float32')
    transitDuration = np.array([], dtype='float32')
    transitPhase = np.array([], dtype='float32')
    dPeriod = (maxper - minper) / float(nsearch)
    trialPeriods = np.arange(minper, maxper + dPeriod, dPeriod, dtype='float32')
    print(' ')
    for trialPeriod in tqdm(trialPeriods):
        srMax = np.append(srMax, 0.0)
        transitDuration = np.append(transitDuration, np.nan)
        transitPhase = np.append(transitPhase, np.nan)
        trialFrequency = 1.0 / trialPeriod

        # minimum and maximum transit durations in quantized phase units
        duration1 = max(int(nbins * mindur / 24.0 / trialPeriod), 2)
        duration2 = max(int(nbins * maxdur / 24.0 / trialPeriod) + 1,
                        duration1 + 1)

        # 30minutes in quantized phase units
        halfHour = int(0.02083333 / trialPeriod * nbins) + 1

        # compute folded time series with trial period
        work4 = np.zeros(nbins, dtype='float32')
        work5 = np.zeros(nbins, dtype='float32')
        phase = np.array(((time_arr * trialFrequency)
                         - np.floor(time_arr * trialFrequency)) * nbins,
                         dtype='int')
        ptuple = np.array([phase, flux_arr, inerr])
        ptuple = np.rot90(ptuple, 3)
        phsort = np.array(sorted(ptuple, key=lambda ph: ph[2]))
        for i in range(nbins):
            elements = np.nonzero(phsort[:, 2] == float(i))[0]
            work4[i] = np.nanmean(phsort[elements, 1])
            work5[i] = (math.sqrt(np.sum(np.power(phsort[elements, 0], 2))
                        / len(elements)))

        # extend the work arrays beyond nbins by wrapping
        work4 = np.append(work4, work4[:duration2])
        work5 = np.append(work5, work5[:duration2])

        # calculate weights of folded light curve points
        sigmaSum = np.nansum(np.power(work5, -2))
        omega = np.power(work5,-2) / sigmaSum

        # calculate weighted phased light curve
        s = omega * work4

        # iterate through trial period phase
        for i1 in range(nbins):
            # iterate through transit durations
            for duration in range(duration1, duration2 + 1, int(halfHour)):
                # calculate maximum signal residue
                i2 = i1 + duration
                sr1 = np.sum(np.power(s[i1:i2], 2))
                sr2 = np.sum(omega[i1:i2])
                sr = math.sqrt(sr1 / (sr2 * (1.0 - sr2)))
                if sr > srMax[-1]:
                    srMax[-1] = sr
                    transitDuration[-1] = float(duration)
                    transitPhase[-1] = float((i1 + i2) / 2)

    # normalize maximum signal residue curve
    bestSr = np.max(srMax)
    bestTrial = np.nonzero(srMax == bestSr)[0][0]
    srMax /= bestSr
    transitDuration *= trialPeriods / 24.0
    BJD0 = (np.array(transitPhase * trialPeriods / nbins,dtype='float64')
            + intime[0] - 2454833.0)
    print('\n')

    # clean up x-axis unit
    ptime = copy(trialPeriods)
    xlab = 'Trial Period (days)'

    # clean up y-axis units
    pout = copy(srMax)
    ylab = 'Normalized Signal Residue'

    # data limits
    xmin = ptime.min()
    xmax = ptime.max()
    ymin = pout.min()
    ymax = pout.max()
    xr = xmax - xmin
    yr = ymax - ymin
    ptime = np.insert(ptime, [0], [ptime[0]])
    ptime = np.append(ptime, [ptime[-1]])
    pout = np.insert(pout, [0], [0.0])
    pout = np.append(pout, 0.0)

    # plot light curve
    if plot:
        plt.figure(figsize=[16, 8])
        plt.clf()
        # plot data
        ax = plt.axes([0.06,0.10,0.93,0.87])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90)
        plt.plot(ptime[1:-1], pout[1:-1], color='#0000ff', linestyle='-',
                 linewidth=1.0)
        plt.fill(ptime, pout, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()
        plt.xlim(xmin-xr*0.01, xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01, ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10, ymax+yr*0.01)

        # render plot
        plt.show()

    # append new BLS data extension to the output file
    col1 = pyfits.Column(name='PERIOD',format='E', unit='days',
                         array=trialPeriods)
    col2 = pyfits.Column(name='BJD0', format='D', unit='BJD - 2454833',
                         array=BJD0)
    col3 = pyfits.Column(name='DURATION', format='E', unit='hours',
                         array=transitDuration)
    col4 = pyfits.Column(name='SIG_RES', format='E', array=srMax)
    cols = pyfits.ColDefs([col1, col2, col3, col4])
    instr.append(pyfits.BinTableHDU.from_columns(cols))
    instr[-1].header.cards['TTYPE1'].comment = 'column title: trial period'
    instr[-1].header.cards['TTYPE2'].comment = 'column title: trial mid-transit zero-point'
    instr[-1].header.cards['TTYPE3'].comment = 'column title: trial transit duration'
    instr[-1].header.cards['TTYPE4'].comment = 'column title: normalized signal residue'
    instr[-1].header.cards['TFORM1'].comment = 'column type: float32'
    instr[-1].header.cards['TFORM2'].comment = 'column type: float64'
    instr[-1].header.cards['TFORM3'].comment = 'column type: float32'
    instr[-1].header.cards['TFORM4'].comment = 'column type: float32'
    instr[-1].header.cards['TUNIT1'].comment = 'column units: days'
    instr[-1].header.cards['TUNIT2'].comment = 'column units: BJD - 2454833'
    instr[-1].header.cards['TUNIT3'].comment = 'column units: hours'
    instr[-1].header['EXTNAME' ] = ('BLS', 'extension name')
    instr[-1].header['PERIOD'  ] = (trialPeriods[bestTrial], 'most significant trial period [d]')
    instr[-1].header['BJD0'    ] = (BJD0[bestTrial] + 2454833.0, 'time of mid-transit [BJD]')
    instr[-1].header['TRANSDUR'] = (transitDuration[bestTrial], 'transit duration [hours]')
    instr[-1].header['SIGNRES' ] = (srMax[bestTrial] * bestSr, 'maximum signal residue')
    # history keyword in output file
    print("Writing output file {}...".format(outfile))
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    instr.writeto(outfile)
    # close input file
    instr.close()
    # print best trial period results
    print('      Best trial period = {} days'.format(trialPeriods[bestTrial]))
    print('    Time of mid-transit = BJD {}'.format(BJD0[bestTrial] + 2454833.0))
    print('       Transit duration = {} hours'.format(transitDuration[bestTrial]))
    print(' Maximum signal residue = {} \n'.format(srMax[bestTrial] * bestSr))

    # end time
    kepmsg.clock('KEPBLS completed at', logfile, verbose)

def kepbls_main():
    import argparse
    parser = argparse.ArgumentParser(
            description=('Perform Box-Least Square searches for periodic'
                         ' exoplanet transits'),
            formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepbls.'),
                        default=None)
    parser.add_argument('--datacol', default='DETSAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--errcol', default='DETSAP_FLUX_ERR',
                        help='Name of data error column to plot', type=str)
    parser.add_argument('--minper', default=1.0,
                        help='Minimum search period [days]', type=float)
    parser.add_argument('--maxper', default=30.0,
                        help='Maximum search period [days]', type=float)
    parser.add_argument('--mindur', default=0.5,
                        help='Minimum transit duration [hours]', type=float)
    parser.add_argument('--maxdur', default=12.0,
                        help='Maximum transit duration [hours]', type=float)
    parser.add_argument('--nsearch', default=1000,
                        help='Number of test periods between minper and maxper',
                        type=int)
    parser.add_argument('--nbins', default=1000,
                        help='Number of bins in the folded time series at any test period',
                        type=int)
    parser.add_argument('--plot', action='store_true',
                        help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepbls.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepbls(args.infile, args.outfile, args.datacol, args.errcol, args.minper,
           args.maxper, args.mindur, args.maxdur, args.nsearch, args.nbins,
           args.plot, args.overwrite, args.verbose, args.logfile)
