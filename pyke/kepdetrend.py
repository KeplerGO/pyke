from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from . import kepfunc
from . import kepio
from . import kepmsg
from . import kepkey
from . import kepfit
from . import kepstat


__all__ = ['kepdetrend']


def kepdetrend(infile, ranges1, ranges2, npoly1, npoly2, nsig1, nsig2,
               niter1, niter2, outfile=None, datacol='SAP_FLUX',
               errcol='SAP_FLUX_ERR', popnans=False, plot=False,
               overwrite=False, verbose=False, logfile='kepdetrend.log'):
    """
    kepdetrend -- Detrend aperture photometry data

    Simple Aperture Photometry (SAP) data can contain a series of systematic
    trends associated with the spacecraft, detector and environment rather than
    the target. Within the Kepler pipeline these contaminants are treated
    during Pre-search Data Conditioning (PDC) and cleaned data are provided in
    the archived files as the PDCSAP_FLUX data. See the Kepler Data
    Characteristics Handbook for more precise descriptions of systematics. The
    Kepler pipeline attempts to remove systematics with a combination of data
    detrending and cotrending against weighted cotrending basis vectors derived
    from the time series structure most-common to all neighbors of the
    scientific target. These processes are imperfect but tackled within the
    pipeline in the spirit of correcting as many targets as possible with
    enough accuracy for the mission to meet exoplanet detection specifications.
    This approach is, however, not optimized for individual targets. Users of
    the Kepler archive may well find that individually-tailored detrending and
    cotrending yields corrected light curves more suited to their science. The
    purpose of kepdetrend is to provide a detrending algorithm that can be
    tailored for individual targets. We stress that an absolute correction is
    often impossible for Kepler data. We also recommend the use of kepcotrend
    instead of kepdetrend for systematic removal in most Kepler targets.

    The current version of this task asks the user to define data ranges that
    are free of a systematic feature that needs to be removed and can be
    well-characterized by a single polynomial function. This function is what
    the task attempts to correct data to. The user then defines a data range
    that needs correcting and fits this with a second polynomial. The
    correction is the subtraction of one polynomial from the other in the data
    range to be detrended. The examples plotted below show the piecemeal
    correction of three systematic features within a Q2 light curve. These
    three corrections are provided in the task examples below.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file. outfile will be an amended version of
        infile with specified time ranges detrended by subtraction of a
        polynomial fit.
    datacol : str
        The column name containing data stored within extension 1 of infile.
        This data will be detrended. Typically this name is SAP_FLUX (Simple
        Aperture Photometry fluxes), but any data column within extension 1 of
        the FITS file can be corrected.
    errcol : str
        The uncertainty data coupled to datacol. Typically this column is
        called SAP_FLUX_ERR. If no errors are associated with datacol then
        use ``errcol=None``.
    ranges1, ranges2 : list, list
        Time ranges are supplied as comma-separated pairs of Barycentric Julian
        Dates (BJDs). Multiple ranges are separated by a semi-colon. An example
        containing two time ranges is::

            '2455012.48517,2455014.50072;2455022.63487,2455025.08231'

        Data within the range **ranges1** will be detrended by subtracting the
        difference between the best fit to data in that range and the best fit
        function in the range **ranges2** extrapolated into ranges1.
    npoly1, npoly2 : int
        The polynomial order for the function that fits the data ranges to be
        detrended.
    nsig1, nsig2 : float, float
        The data to be detrended is fit by a polynomial using an iterative
        scheme. After a best fit is found, those data points deviating from the
        fit by more than this specified amount are rejected and the remaining
        data are fit again, etc, until there are no further rejections. This
        threshold is in units of the standard deviation of the data about the
        best fit function.
    niter1, niter2 : int, int
        The polynomial fit over the data to be detrended will be iterated until
        there are no further outlier rejections or the number of iterations
        exceeds **niter1** (**niter2**).
    popnans : bool
        Keep NaN flux values (times without a flux measurement) in the output
        FITS file. If set to no, any rows in the input FITS file containing
        NaNs will no be in the output file.
    plot : boolean
        Plot the data, the fits and the correction?
    overwrite : bool
        Overwrite the output file? if **overwrite** is **False** and an
        existing file has the same name as outfile then the task will stop with
        an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepdetrend kplr002436324-2009259160929_llc.fits --datacol SAP_FLUX
        --errcol SAP_FLUX_ERR --ranges1 '2455063.59357,2455066.47292' --npoly1 5
        --nsig1 3 --niter1 10 --ranges2 '2455060.57026,2455063.59357;2455066.9768,2455068.99235'
        --npoly2 5 --nsig2 3 --niter2 10 --plot --verbose

    .. image:: ../_static/images/api/kepdetrend.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPDETREND -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' errcol={}'.format(errcol)
            + ' ranges1={}'.format(ranges1)
            + ' npoly1={}'.format(npoly1)
            + ' nsig1={}'.format(nsig1)
            + ' niter1={}'.format(niter1)
            + ' ranges2={}'.format(ranges2)
            + ' npoly2={}'.format(npoly2)
            + ' nsig2={}'.format(nsig2)
            + ' niter2={}'.format(niter2)
            + ' popnans={}'.format(popnans)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPDETREND started at',logfile,verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPDETREND: {} exists. Use --overwrite'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # filter input data table
    work1 = np.array([table.field('time'), table.field(datacol),
                      table.field(errcol)])
    work1 = np.rot90(work1, 3)
    work1 = work1[~np.isnan(work1).any(1)]

    # read table columns
    intime = work1[:, 2] + bjdref
    indata = work1[:, 1]
    inerr = work1[:, 0]

    # time ranges for region 1 (region to be corrected)
    time1, data1, err1 = [], [], []
    t1start, t1stop = kepio.timeranges(ranges1, logfile, verbose)
    cadencelis1 = kepstat.filterOnRange(intime, t1start, t1stop)
    for i in range(len(cadencelis1)):
        time1.append(intime[cadencelis1[i]])
        data1.append(indata[cadencelis1[i]])
        if errcol.lower() != 'none':
            err1.append(inerr[cadencelis1[i]])
    t0 = time1[0]
    time1 = np.array(time1, dtype='float64') - t0
    data1 = np.array(data1, dtype='float32')
    if errcol.lower() != 'none':
        err1 = np.array(err1, dtype='float32')
    else:
        err1 = None

    # fit function to range 1
    fit_func = getattr(kepfunc, 'poly' + str(npoly1))
    pinit = [data1.mean()]
    if npoly1 > 0:
        for i in range(npoly1):
            pinit.append(0)
    pinit = np.array(pinit, dtype='float32')
    coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx1, ploty1 = \
        kepfit.lsqclip(fit_func, pinit, time1, data1, err1, nsig1, nsig1,
                       niter1, logfile, verbose)
    fit1 = indata * 0.0
    for i in range(len(coeffs)):
        fit1 += coeffs[i] * (intime - t0) ** i
    for i in range(len(intime)):
        if i not in cadencelis1:
            fit1[i] = 0.0
    plotx1 += t0

    # time ranges for region 2 (region that is correct)
    time2, data2, err2 = [], [], []
    t2start, t2stop = kepio.timeranges(ranges2, logfile, verbose)
    cadencelis2 = kepstat.filterOnRange(intime, t2start, t2stop)
    for i in range(len(cadencelis2)):
        time2.append(intime[cadencelis2[i]])
        data2.append(indata[cadencelis2[i]])
        if errcol.lower() != 'none':
            err2.append(inerr[cadencelis2[i]])
    t0 = time2[0]
    time2 = np.array(time2, dtype='float64') - t0
    data2 = np.array(data2, dtype='float32')
    if errcol.lower() != 'none':
        err2 = np.array(err2, dtype='float32')
    else:
        err2 = None

    # fit function to range 2
    fit_func = getattr(kepfunc, 'poly' + str(npoly2))
    pinit = [data2.mean()]
    if npoly2 > 0:
        for i in range(npoly2):
            pinit.append(0)
    pinit = np.array(pinit, dtype='float32')
    coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx2, ploty2 = \
        kepfit.lsqclip(fit_func, pinit, time2, data2, err2, nsig2, nsig2,
                       niter2, logfile, verbose)
    fit2 = indata * 0.0
    for i in range(len(coeffs)):
        fit2 += coeffs[i] * (intime - t0) ** i
    for i in range(len(intime)):
        if i not in cadencelis1:
            fit2[i] = 0.0
    plotx2 += t0

    # normalize data
    outdata = indata - fit1 + fit2
    if errcol.lower() != 'none':
        outerr = inerr * 1.0

    # comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)

    # clean up x-axis unit
    intime0 = tstart // 100 * 100.0
    if intime0 < 2.4e6:
        intime0 += 2.4e6
    ptime = intime - intime0
    plotx1 = plotx1 - intime0
    plotx2 = plotx2 - intime0
    xlab = 'BJD $-$ {}'.format(intime0)

    # clean up y-axis units
    pout = outdata
    ploty1
    ploty2
    nrm = len(str(int(np.nanmax(indata))))-1
    indata = indata / 10 ** nrm
    pout = pout / 10 ** nrm
    ploty1 = ploty1 / 10 ** nrm
    ploty2 = ploty2 / 10 ** nrm
    ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

    # data limits
    xmin = ptime.min()
    xmax = ptime.max()
    ymin = indata.min()
    ymax = indata.max()
    omin = pout.min()
    omax = pout.max()
    xr = xmax - xmin
    yr = ymax - ymin
    oo = omax - omin
    ptime = np.insert(ptime, [0], [ptime[0]])
    ptime = np.append(ptime, [ptime[-1]])
    indata = np.insert(indata, [0], [0.0])
    indata = np.append(indata, [0.0])
    pout = np.insert(pout, [0], [0.0])
    pout = np.append(pout, 0.0)

    # plot light curve
    if plot:
        plt.figure()
        plt.clf()
        # plot original data
        ax = plt.axes([0.06,0.523,0.93,0.45])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.plot(ptime, indata, color='#0000ff', linestyle='-', linewidth=1.0)
        plt.fill(ptime, indata, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.plot(plotx1, ploty1, color='r', linestyle='-', linewidth=2.0)
        plt.plot(plotx2, ploty2, color='g', linestyle='-', linewidth=2.0)
        plt.xlim(xmin-xr*0.01, xmax+xr*0.01)
        if ymin > 0.0:
            plt.ylim(ymin-yr*0.01, ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10, ymax+yr*0.01)
            plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()
        # plot detrended data
        ax = plt.axes([0.06, 0.073, 0.93, 0.45])

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.plot(ptime, pout, color='#0000ff', linestyle='-', linewidth=1.0)
        plt.fill(ptime, pout, color=fcolor, linewidth=0.0, alpha=0.2)
        plt.xlim(xmin-xr*0.01, xmax+xr*0.01)
        if ymin > 0.0:
            plt.ylim(omin-oo*0.01, omax+oo*0.01)
        else:
            plt.ylim(1.0e-10, omax+oo*0.01)
        plt.xlabel(xlab, {'color' : 'k'})
        try:
            plt.ylabel(ylab, {'color' : 'k'})
        except:
            ylab = '10**{} e-/s'.format(nrm)
            plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

    # render plot
    plt.show()
    # write output file
    print("Writing output file {}...".format(outfile))
    if popnans:
        instr[1].data.field(datacol)[good_data] = outdata
        instr[1].data.field(errcol)[good_data] = outerr
        instr[1].data.field(datacol)[bad_data] = None
        instr[1].data.field(errcol)[bad_data] = None
        instr.writeto(outfile)
    elif not popnans:
        for i in range(len(outdata)):
            instr[1].data.field(datacol)[i] = outdata[i]
            if errcol.lower() != 'none':
                instr[1].data.field(errcol)[i] = outerr[i]
        instr.writeto(outfile)
    # close input file
    instr.close()
    ## end time
    kepmsg.clock('KEPDETREND completed at', logfile, verbose)

def kepdetrend_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Detrend systematic features from Simple Aperture'
                          ' Photometry (SAP) data'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepdetrend.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column', type=str)
    parser.add_argument('--errcol', default='SAP_FLUX_ERR',
                        help='Name of data error column', type=str)
    parser.add_argument('--ranges1', help='Time ranges of region 1', type=str)
    parser.add_argument('--npoly1', help='Polynomial order for region 1',
                        type=int)
    parser.add_argument('--nsig1',
                        help='Sigma clipping threshold for region 1', type=int)
    parser.add_argument('--niter1',
                        help='Maximum number of clipping iterations for region 1',
                        type=int)
    parser.add_argument('--ranges2', help='Time ranges of region 2', type=str)
    parser.add_argument('--npoly2', help='Polynomial order for region 2',
                        type=int)
    parser.add_argument('--nsig2', help='Sigma clipping threshold for region 2',
                        type=int)
    parser.add_argument('--niter2',
                        help='Maximum number of clipping iterations for region 2',
                        type=int)
    parser.add_argument('--popnans', action='store_true',
                        help='Keep cadences with no flux value?')
    parser.add_argument('--plot', action='store_true',
                        help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='kepdetrend.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepdetrend(args.infile, args.ranges1, args.ranges2,
               args.npoly1, args.npoly2, args.nsig1, args.nsig2,
               args.niter1, args.niter2, args.outfile, args.datacol, args.errcol,
               args.popnans, args.plot, args.overwrite, args.verbose,
               args.logfile)
