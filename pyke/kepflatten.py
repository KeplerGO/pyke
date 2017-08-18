from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey, kepfit, kepstat, kepfunc
import re
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
from astropy.io import fits as pyfits
from tqdm import tqdm


__all__ = ['kepflatten']


def kepflatten(infile, outfile=None, datacol='PDCSAP_FLUX',
               errcol='PDCSAP_FLUX_ERR', nsig=3., stepsize=0.5, winsize=5.0,
               npoly=3, niter=1, ranges='0,0', plot=False, overwrite=False,
               verbose=False, logfile='kepflatten.log'):
    """
    kepflatten -- Remove low frequency variability from time-series, preserve
    transits and flares

    kepflatten detrends data for low-frequency photometric structure by
    dividing by the mean of best-fit sliding polynomials over a sequential
    series of small time ranges across the data. For example, a typical
    timestamp is fit three times of ``stepsize=1.0`` and ``winsize=3.0``. The
    adopted fit to the timestamp will be the mean of the three values. Outliers
    are iteratively-clipped from the fit, therefore structure in e.g.
    short-lived transits or flares are better preserved compared to e.g.
    bandpass filtering methods (``kepfilter``). Optionally, input data, best
    fits, fit outliers and output data are rendered to a plot window. In many
    respects kepflatten performs the opposite task to kepoutlier which removes
    statistical outliers while preserving low-frequency structure in light
    curves.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file. outfile will be a direct copy of
        infile but with NaN timestamps removed and two new columns in the 1st
        extension - DETSAP_FLUX (a flattened or detrended for low-frequency
        variations version of the data) and DETSAP_FLUX_ERR (the associated
        1-:math:`\sigma` error).
    datacol : str
        The column name containing data stored within extension 1 of infile.
        Typically this name is SAP_FLUX (Simple Aperture Photometry fluxes),
        PDCSAP_FLUX (Pre-search Data Conditioning fluxes) or CBVSAP_FLUX
        (SAP_FLUX corrected for systematic artifacts by the PyKE tool
        kepcotrend).
    errcol : str
        The column name containing photometric 1-:math:`\sigma` errors
        stored within extension 1 of infile. Typically this name is
        SAP_FLUX_ERR (Simple Aperture Photometry fluxes), PDCSAP_FLUX_ERR
        (Pre-search Data Conditioning fluxes). The error column coupled to
        CBVSAP_FLUX data is SAP_FLUX_ERR. kepflatten normalizes datacol and
        errcol consistently using a series of best fit polynomials.
    nsig : float
        The sigma clipping threshold in units of standard deviation. Data
        deviating from a best fit function by more than the threshold will
        ignored during subsequent fit iterations.
    stepsize : float
        The data within datacol is unlikely to be well represented by a single
        polynomial function. stepsize splits the data up into a series of time
        blocks, each is fit independently by a separate function. The user can
        provide an informed choice of stepsize after inspecting the data with
        the kepdraw tool. Units are days.
    winsize : float
        The size of the window to be fit during each step. Units are days.
        winsize must be greater or equal to stepsize. winsize >> stepsize is
        recommended.
    npoly : integer
        The order of each piecemeal polynomial function.
    niter : integer
        If outliers outside of nsig are found in a particular data section,
        that data will be removed temporarily and the time series fit again.
        This will be iterated niter times before freezing upon the best current
        fit.
    ranges : str
        The user can choose specific time ranges of data on which to work. This
        could, for example, avoid removing known stellar flares from a dataset.
        Time ranges are supplied as comma-separated pairs of Barycentric Julian
        Dates (BJDs). Multiple ranges are separated by a semi-colon. An example
        containing two time ranges is:

        ``2455012.48517,2455014.50072;2455022.63487,2455025.08231``.

        If the user wants to correct the entire time series then providing
        ``ranges='0,0'`` will tell the task to operate on the whole time
        series.
    plot : bool
        Plot the data, fit, outliers and result?
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepflatten kplr012557548-2011177032512_llc.fits
        --nsig 3 --stepsize 1.0 --winsize 3.0 --npoly 3 --niter 10 --plot
        --overwrite --verbose

    .. image:: ../_static/images/api/kepflatten.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPFLATTEN -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' errcol={}'.format(errcol)
            + ' nsig={}'.format(nsig)
            + ' stepsize={}'.format(stepsize)
            + ' winsize={}'.format(winsize)
            + ' npoly={}'.format(npoly)
            + ' niter={}'.format(niter)
            + ' ranges={}'.format(ranges)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPFLATTEN started at', logfile, verbose)

    # test winsize > stepsize
    if winsize < stepsize:
        errmsg = 'ERROR -- KEPFLATTEN: winsize must be greater than stepsize'
        kepmsg.err(logfile, errmsg, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPFLATTEN: {} exists. Use overwrite=True'
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
    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # filter input data table
    try:
        datac = table.field(datacol)
    except:
        errmsg = ('ERROR -- KEPFLATTEN: cannot find or read data column {}'
                  .format(datacol))
        kepmsg.err(logfile, message, verbose)
    try:
        err = table.field(errcol)
    except:
        errmsg = ('WARNING -- KEPFLATTEN: cannot find or read error column {}'
                  .format(errcol))
        kepmsg.warn(logfile, errmsg, verbose)
        errcol = 'None'
    if errcol.lower() == 'none' or errcol == 'PSF_FLUX_ERR':
        err = datac * cadence
        err = np.sqrt(np.abs(err)) / cadence
        work1 = np.array([table.field('time'), datac, err])
    else:
        work1 = np.array([table.field('time'), datac, err])
    work1 = np.rot90(work1, 3)
    work1 = work1[~np.isnan(work1).any(1)]

    # read table columns
    intime = work1[:, 2] + bjdref
    indata = work1[:, 1]
    inerr = work1[:, 0]
    if len(intime) == 0:
         message = 'ERROR -- KEPFLATTEN: one of the input arrays is all NaN'
         kepmsg.err(logfile, message, verbose)

    # time ranges for region to be corrected
    t1, t2 = kepio.timeranges(ranges, logfile, verbose)
    cadencelis = kepstat.filterOnRange(intime, t1, t2)
    # find limits of each time step
    tstep1, tstep2 = [], []
    work = intime[0]
    while work <= intime[-1]:
        tstep1.append(work)
        tstep2.append(np.array([work + winsize, intime[-1]],
                               dtype='float64').min())
        work += stepsize

    # find cadence limits of each time step
    cstep1, cstep2 = [], []
    for n in range(len(tstep1)):
        for i in range(len(intime)-1):
            if intime[i] <= tstep1[n] and intime[i+1] > tstep1[n]:
                for j in range(i, len(intime)-1):
                    if intime[j] < tstep2[n] and intime[j+1] >= tstep2[n]:
                        cstep1.append(i)
                        cstep2.append(j + 1)

    # comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    # clean up x-axis unit
    intime0 = tstart // 100 * 100.0
    ptime = intime - intime0
    xlab = 'BJD $-$ {}'.format(intime0)

    # clean up y-axis units
    pout = copy(indata)
    nrm = len(str(int(pout.max()))) - 1
    pout = pout / 10 ** nrm
    ylab = '10$^{}$'.format(nrm) + 'e$^-$ s$^{-1}$'

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

    if plot:
        plt.figure()
        plt.clf()
        # plot data
        ax = plt.axes([0.06, 0.54, 0.93, 0.43])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(plt.gca(), xticklabels=[])
        plt.plot(ptime[1:-1], pout[1:-1], color='#363636', linestyle='-',
                 linewidth=1.0)
        plt.fill(ptime, pout, color='#a8a7a7', linewidth=0.0, alpha=0.2)
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

    # loop over each time step, fit data, determine rms
    fitarray = np.zeros((len(indata), len(cstep1)), dtype='float32')
    sigarray = np.zeros((len(indata), len(cstep1)), dtype='float32')
    fitarray[:, :] = np.nan
    sigarray[:, :] = np.nan
    masterfit = indata * 0.0
    mastersigma = np.zeros(len(masterfit))
    functype = getattr(kepfunc, 'poly' + str(npoly))
    for i in tqdm(range(len(cstep1))):
        timeSeries = intime[cstep1[i]:cstep2[i]+1] - intime[cstep1[i]]
        dataSeries = indata[cstep1[i]:cstep2[i]+1]
        fitTimeSeries = np.array([], dtype='float32')
        fitDataSeries = np.array([], dtype='float32')
        pinit = [dataSeries.mean()]
        if npoly > 0:
            for j in range(npoly):
                pinit.append(0.0)
        pinit = np.array(pinit, dtype='float32')
        try:
            if len(fitarray[cstep1[i]:cstep2[i]+1,i]) > len(pinit):
                coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                    kepfit.lsqclip(functype, pinit, timeSeries, dataSeries,
                                   None, nsig, nsig, niter, logfile, verbose)
                fitarray[cstep1[i]:cstep2[i]+1, i] = 0.0
                sigarray[cstep1[i]:cstep2[i]+1, i] = sigma
                for j in range(len(coeffs)):
                    fitarray[cstep1[i]:cstep2[i]+1, i] += coeffs[j] * timeSeries ** j
        except:
            message  = ('WARNING -- KEPFLATTEN: could not fit range '
                        + str(intime[cstep1[i]]) + '-' + str(intime[cstep2[i]]))
            kepmsg.warn(logfile, message, verbose)

    # find mean fit for each timestamp
    for i in range(len(indata)):
        masterfit[i] = np.nanmean(fitarray[i, :])
        mastersigma[i] = np.nanmean(sigarray[i, :])
    masterfit[-1] = masterfit[-4] #fudge
    masterfit[-2] = masterfit[-4] #fudge
    masterfit[-3] = masterfit[-4] #fudge
    plt.plot(intime - intime0, masterfit / 10 ** nrm, 'b')

    # reject outliers
    rejtime, rejdata = [], []
    naxis2 = 0
    for i in range(len(masterfit)):
        if (abs(indata[i] - masterfit[i]) > nsig * mastersigma[i]
            and i in cadencelis):
            rejtime.append(intime[i])
            rejdata.append(indata[i])
    rejtime = np.array(rejtime, dtype='float64')
    rejdata = np.array(rejdata, dtype='float32')
    if plot:
        plt.plot(rejtime - intime0, rejdata / 10 ** nrm, 'ro', markersize=2)
    # new data for output file
    outdata = indata / masterfit
    outerr = inerr / masterfit
    pout = copy(outdata)
    ylab = 'Normalized Flux'
    # plot ranges
    if plot:
        plt.xlim(xmin-xr*0.01, xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01, ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10, ymax+yr*0.01)
        # plot residual data
        ax = plt.axes([0.06,0.09,0.93,0.43])
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()

        ymin = pout.min()
        ymax = pout.max()
        yr = ymax - ymin
        pout = np.insert(pout, [0], [0.0])
        pout = np.append(pout, 0.0)
        plt.plot(ptime[1:-1], pout[1:-1], color='#363636', linestyle='-',
                 linewidth=1.0)
        plt.fill(ptime, pout, color='#a8a7a7', linewidth=0.0, alpha=0.2)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()
        # plot ranges
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin >= 0.0:
            plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
        else:
            plt.ylim(1.0e-10, ymax + yr * 0.01)
        # render plot
        plt.savefig(re.sub('.fits', '.png', outfile))
        plt.show()
    # add NaNs back into data
    n = 0
    work1 = np.array([], dtype='float32')
    work2 = np.array([], dtype='float32')
    instr = pyfits.open(infile, 'readonly')
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    tn = table.field('time')
    dn = table.field(datacol)
    for i in range(len(table.field(0))):
        if np.isfinite(tn[i]) and np.isfinite(dn[i]) and np.isfinite(err[i]):
            try:
                work1 = np.append(work1, outdata[n])
                work2 = np.append(work2, outerr[n])
                n += 1
            except:
                pass
        else:
            work1 = np.append(work1, np.nan)
            work2 = np.append(work2, np.nan)

    # history keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)

    # write output file
    try:
        print("Writing output file {}...".format(outfile))
        col1 = pyfits.Column(name='DETSAP_FLUX',format='E13.7',array=work1)
        col2 = pyfits.Column(name='DETSAP_FLUX_ERR',format='E13.7',array=work2)
        cols = instr[1].data.columns + col1 + col2
        instr[1] = pyfits.BinTableHDU.from_columns(cols, header=instr[1].header)
        instr.writeto(outfile)
    except ValueError:
        try:
            instr[1].data.field('DETSAP_FLUX')[:] = work1
            instr[1].data.field('DETSAP_FLUX_ERR')[:] = work2
            instr.writeto(outfile)
        except:
            message = ('ERROR -- KEPFLATTEN: cannot add DETSAP_FLUX data to '
                       'FITS file')
            kepmsg.err(logfile, message, verbose)

    # close input file
    instr.close()
    ## end time
    kepmsg.clock('KEPFLATTEN completed at', logfile, verbose)

def kepflatten_main():
    import argparse

    parser = argparse.ArgumentParser(
             description=('Remove low frequency variability from time-series,'
                          'preserve transits and flares'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepflatten.'),
                        default=None)
    parser.add_argument('--datacol', default='PDCSAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--errcol', default='PDCSAP_FLUX_ERR',
                        help='Name of data error column to plot', type=str)
    parser.add_argument('--nsig', default=3.,
                        help='Sigma clipping threshold for outliers',
                        type=float)
    parser.add_argument('--stepsize', default=0.5,
                        help='Stepsize on which to fit data [days]',
                        type=float)
    parser.add_argument('--winsize', default=5.0,
                        help=('Window size of data to fit after each step'
                              ' (>= stepsize) [days]'),
                        type=float)
    parser.add_argument('--npoly', default=3,
                        help='Polynomial order for each fit', type=int)
    parser.add_argument('--niter', default=1,
                        help='Maximum number of clipping iterations', type=int)
    parser.add_argument('--ranges', default='0,0',
                        help='Time ranges of regions to filter',
                        type=str)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepflatten.log', dest='logfile', type=str)
    args = parser.parse_args()

    kepflatten(args.infile, args.outfile, args.datacol, args.errcol, args.nsig,
               args.stepsize, args.winsize, args.npoly, args.niter,
               args.ranges, args.plot, args.overwrite, args.verbose,
               args.logfile)
