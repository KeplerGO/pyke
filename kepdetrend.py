import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from . import kepfunc
from . import kepio
from . import kepmsg
from . import kepkey
from . import kepfit
from . import kepstat

def kepdetrend(infile, outfile, datacol, errcol, ranges1, npoly1, nsig1,
               niter1, ranges2, npoly2, nsig2, niter2, popnans, plot,
               clobber, verbose, logfile):
    """
    """

    # startup parameters
    labelsize = 24
    ticksize = 16
    xsize = 16
    ysize = 9
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

    # log the call
    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPDETREND -- '
            + 'infile= {}'.format(infile)
            + 'outfile= {}'.format(outfile)
            + 'datacol= {}'.format(datacol)
            + 'errcol= {}'.format(errcol)
            + 'ranges1= {}'.format(ranges1)
            + 'npoly1= {}'.format(npoly1)
            + 'nsig1= {}'.format(nsig1)
            + 'niter1= {}'.format(niter1)
            + 'ranges2= {}'.format(ranges2)
            + 'npoly2= {}'.format(npoly2)
            + 'nsig2= {}'.format(nsig2)
            + 'niter2= {}'.format(niter2)
            + 'popnans= {}'.format(popnans)
            + 'plot= {}'.format(plot)
            + 'clobber= {}'.format(clobber)
            + 'verbose= {}'.format(verbose)
            + 'logfile= {}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPDETREND started at',logfile,verbose)
    # clobber output file
    if clobber:
        kepio.clobber(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        message = ('ERROR -- KEPDETREND: {} exists. Use clobber=True'
                   .format(outfile))
        status = kepmsg.err(logfile, message, verbose)

    # open input file
    instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
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
    intime = work1[:,2] + bjdref
    indata = work1[:,1]
    inerr = work1[:,0]

    # time ranges for region 1 (region to be corrected)
    time1 = []; data1 = []; err1 = []
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
        fit2 += coeffs[i] * (intime - t0)**i
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
    indata = indata / 10**nrm
    pout = pout / 10**nrm
    ploty1 = ploty1 / 10**nrm
    ploty2 = ploty2 / 10**nrm
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
        plt.figure(figsize=[xsize, ysize])
        plt.clf()
        # plot original data
        ax = plt.axes([0.06,0.523,0.93,0.45])
        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=12)
        plt.plot(ptime, indata, color=lcolor, linestyle='-', linewidth=lwidth)
        plt.fill(ptime, indata, color=fcolor, linewidth=0.0, alpha=falpha)
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
        ax = plt.axes([0.06,0.073,0.93,0.45])

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=12)
        plt.plot(ptime, pout, color=lcolor, linestyle='-', linewidth=lwidth)
        plt.fill(ptime, pout, color=fcolor, linewidth=0.0, alpha=falpha)
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

    # render plot
    plt.show()
    # write output file
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
    kepio.closefits(instr, logfile, verbose)
    ## end time
    kepmsg.clock('KEPDETREND completed at', logfile, verbose)

def kepdetrend_main():
    import argparse
    parser = argparse.ArgumentParser(
            description=('Detrend systematic features from Simple Aperture '
                         'Photometry (SAP) data'))
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column', type=str)
    parser.add_argument('--errcol', default='SAP_FLUX_ERR', help='Name of data error column', type=str)

    parser.add_argument('--ranges1', help='Time ranges of region 1', type=str)
    parser.add_argument('--npoly1', help='Polynomial order for region 1', type=int)
    parser.add_argument('--nsig1', help='Sigma clipping threshold for region 1', type=int)
    parser.add_argument('--niter1', help='Maximum number of clipping iterations for region 1', type=int)

    parser.add_argument('--ranges2', help='Time ranges of region 2', type=str)
    parser.add_argument('--npoly2', help='Polynomial order for region 2', type=int)
    parser.add_argument('--nsig2', help='Sigma clipping threshold for region 2', type=int)
    parser.add_argument('--niter2', help='Maximum number of clipping iterations for region 2', type=int)

    parser.add_argument('--popnans', action='store_true', help='Keep cadences with no flux value?')
    parser.add_argument('--plot', action='store_true', help='Plot result?')


    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepdetrend.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepdetrend(args.infile, args.outfile, args.datacol, args.errcol,
               args.ranges1, args.npoly1, args.nsig1, args.niter1,
               args.ranges2, args.npoly2, args.nsig2, args.niter2,
               args.popnans, args.plot, args.clobber, args.verbose,
               args.logfile)
