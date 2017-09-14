from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepfit, kepstat, kepfunc


__all__ = ['kepoutlier']


def kepoutlier(infile, outfile=None, datacol='SAP_FLUX', nsig=3.0, stepsize=1.0,
               npoly=3, niter=1, operation='remove', ranges='0,0', plot=False,
               plotfit=False, overwrite=False, verbose=False,
               logfile='kepoutlier.log'):
    """
    kepoutlier -- Remove or replace statistical outliers from time series data

    kepoutlier identifies data outliers relative to piecemeal best-fit
    polynomials. Outliers are either removed from the output time series or
    replaced by a noise-treated value defined by the polynomial fit. Identified
    outliers and the best fit functions are optionally plotted for inspection
    purposes.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output FITS file. ``outfile`` will be direct copy of
        infile with either data outliers removed (i.e. the table will have
        fewer rows) or the outliers will be corrected according to a best-fit
        function and a noise model.
    datacol : str
        The column name containing data stored within extension 1 of infile.
        This data will be searched for outliers. Typically this name is
        SAP_FLUX (Simple Aperture Photometry fluxes) or PDCSAP_FLUX (Pre-search
        Data Conditioning fluxes).
    nsig : float
        The sigma clipping threshold. Data deviating from a best fit function
        by more than the threshold will be either removed or corrected
        according to the user selection of operation.
    stepsize : float
        The data within datacol is unlikely to be well represented by a single
        polynomial function. stepsize splits the data up into a series of time
        blocks, each is fit independently by a separate function. The user can
        provide an informed choice of stepsize after inspecting the data with
        the kepdraw tool. Units are days.
    npoly : int
        The polynomial order of each best-fit function.
    niter : int
        If outliers are found in a particular data section, that data will be
        removed temporarily and the time series fit again. This will be
        iterated niter times before freezing upon the best available fit.
    operation : str

        * ``remove`` throws away outliers. The output data table will smaller
          or equal in size to the input table.
        * ``replace`` replaces outliers with a value that is consistent with
          the best-fit polynomial function and a random component defined by the
          rms of the data relative to the fit and calculated using the inverse
          normal cumulative function and a random number generator.
    ranges : str
        The user can choose specific time ranges of data on which to work. This
        could, for example, avoid removing known stellar flares from a dataset.
        Time ranges are supplied as comma-separated pairs of Barycentric Julian
        Dates (BJDs). Multiple ranges are separated by a semi-colon. An example
        containing two time ranges is::

            '2455012.48517,2455014.50072;2455022.63487,2455025.08231'

        If the user wants to correct the entire time series then providing
        ``ranges = '0,0'`` will tell the task to operate on the whole time series.
    plot : bool
        Plot the data and outliers?
    plotfit : bool
        Overlay the polynomial fits upon the plot?
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepoutlier kplr002437329-2010355172524_llc.fits --datacol SAP_FLUX
        --nsig 4 --stepsize 5 --npoly 2 --niter 10 --operation replace
        --verbose --plot --plotfit

    .. image:: ../_static/images/api/kepoutlier.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPOUTLIER -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' nsig={}'.format(nsig)
            + ' stepsize={}'.format(stepsize)
            + ' npoly={}'.format(npoly)
            + ' niter={}'.format(niter)
            + ' operation={}'.format(operation)
            + ' ranges={}'.format(ranges)
            + ' plot={}'.format(plot)
            + ' plotfit={}'.format(plotfit)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPOUTLIER started at', logfile, verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPOUTLIER: {} exists. Use overwrite=True'
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
    try:
        nanclean = instr[1].header['NANCLEAN']
    except:
        time = kepio.readtimecol(infile, table, logfile, verbose)
        flux = kepio.readfitscol(infile, table, datacol, logfile, verbose)
        finite_data_mask = np.isfinite(time) & np.isfinite(flux) & (flux != 0)
        table = table[finite_data_mask]
        instr[1].data = table
        comment = 'NaN cadences removed from data'
        kepkey.new('NANCLEAN', True, comment, instr[1], outfile, logfile,
                   verbose)

    # read table columns
    try:
        intime = instr[1].data.field('barytime') + 2.4e6
    except:
        intime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                   verbose)
    indata = kepio.readfitscol(infile, instr[1].data, datacol, logfile,
                               verbose)
    intime = intime + bjdref
    indata = indata / cadenom
    # time ranges for region to be corrected
    t1, t2 = kepio.timeranges(ranges, logfile, verbose)
    cadencelis = kepstat.filterOnRange(intime, t1, t2)

    # find limits of each time step
    tstep1, tstep2 = [], []
    work = intime[0]
    while work < intime[-1]:
        tstep1.append(work)
        tstep2.append(np.array([work + stepsize, intime[-1]],
                               dtype='float64').min())
        work += stepsize

    # find cadence limits of each time step
    cstep1, cstep2 = [], []
    work1 = 0
    work2 = 0
    for i in range(len(intime)):
        if intime[i] >= intime[work1] and intime[i] < intime[work1] + stepsize:
            work2 = i
        else:
            cstep1.append(work1)
            cstep2.append(work2)
            work1 = i
            work2 = i
    cstep1.append(work1)
    cstep2.append(work2)

    outdata = indata * 1.0
    # comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    # clean up x-axis unit
    intime0 = (tstart // 100) * 100.0
    ptime = intime - intime0
    xlab = 'BJD $-$ {}'.format(intime0)

    # clean up y-axis units
    pout = indata * 1.0
    nrm = len(str(int(pout.max())))-1
    pout = pout / 10**nrm
    ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

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
        plt.figure()
        plt.clf()
        # plot data
        ax = plt.axes([0.06, 0.1, 0.93, 0.87])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.plot(ptime, pout, color='#0000ff', linestyle='-', linewidth=1.0)
        plt.fill(ptime, pout, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()
    # loop over each time step, fit data, determine rms
    masterfit = indata * 0.0
    mastersigma = np.zeros(len(masterfit))
    functype = getattr(kepfunc, 'poly' + str(npoly))
    for i in range(len(cstep1)):
        pinit = [indata[cstep1[i]:cstep2[i]+1].mean()]
        if npoly > 0:
            for j in range(npoly):
                pinit.append(0.0)
        pinit = np.array(pinit, dtype='float32')
        try:
            coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit,
                               intime[cstep1[i]:cstep2[i]+1] - intime[cstep1[i]],
                               indata[cstep1[i]:cstep2[i]+1], None, nsig,
                               nsig, niter, logfile, verbose)
            for j in range(len(coeffs)):
                masterfit[cstep1[i]: cstep2[i] + 1] += (coeffs[j]
                        * (intime[cstep1[i]:cstep2[i]+1] - intime[cstep1[i]]) ** j)
            for j in range(cstep1[i], cstep2[i] + 1):
                mastersigma[j] = sigma
            if plotfit:
                plt.plot(plotx + intime[cstep1[i]] - intime0, ploty / 10 ** nrm,
                         'g', lw=3)
        except:
            for j in range(cstep1[i], cstep2[i] + 1):
                masterfit[j] = indata[j]
                mastersigma[j] = 1.0e10
            message = ('WARNING -- KEPOUTLIER: could not fit range '
                       + str(intime[cstep1[i]]) + '-' + str(intime[cstep2[i]]))
            kepmsg.warn(logfile, message, verbose)

    # reject outliers
    rejtime, rejdata = [], []
    naxis2 = 0
    for i in tqdm(range(len(masterfit))):
        if (abs(indata[i] - masterfit[i]) > nsig * mastersigma[i]
            and i in cadencelis):
            rejtime.append(intime[i])
            rejdata.append(indata[i])
            if operation == 'replace':
                [rnd] = kepstat.randarray([masterfit[i]], [mastersigma[i]])
                table[naxis2] = table[i]
                table.field(datacol)[naxis2] = rnd
                naxis2 += 1
        else:
            table[naxis2] = table[i]
            naxis2 += 1
    instr[1].data = table[:naxis2]

    if plot:
        rejtime = np.array(rejtime, dtype='float64')
        rejdata = np.array(rejdata, dtype='float32')
        plt.plot(rejtime - intime0, rejdata / 10 ** nrm, 'ro')
        # plot ranges
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin >= 0.0:
            plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
        else:
            plt.ylim(1.0e-10, ymax + yr * 0.01)

        # render plot
        plt.show()
    # write output file
    print("Writing output file {}...".format(outfile))
    instr.writeto(outfile)
    # close input file
    instr.close()
    kepmsg.clock('KEPOUTLIER completed at', logfile, verbose)

def kepoutlier_main():
    import argparse
    parser = argparse.ArgumentParser(
             description='Remove or replace data outliers from a time series',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepoutlier.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--nsig', default=3.,
                        help='Sigma clipping threshold for outliers',
                        type=float)
    parser.add_argument('--stepsize', default=1.0,
                        help='Stepsize on which to fit data [days]',
                        type=float)
    parser.add_argument('--npoly', default=3,
                        help='Polynomial order for each fit', type=int)
    parser.add_argument('--niter', default=1,
                        help='Maximum number of clipping iterations', type=int)
    parser.add_argument('--operation', default='remove',
                        help='Remove or replace outliers?', type=str,
                        choices=['replace','remove'])
    parser.add_argument('--ranges', default='0,0',
                        help='Time ranges of regions to filter', type=str)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--plotfit', action='store_true',
                        help='Plot fit over results?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepoutlier.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepoutlier(args.infile, args.outfile, args.datacol, args.nsig,
               args.stepsize, args.npoly,args.niter, args.operation,
               args.ranges, args.plot, args.plotfit, args.overwrite,
               args.verbose, args.logfile)
