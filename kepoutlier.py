import numpy as np
import sys, time, re
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from math import *
import kepio, kepmsg, kepkey, kepfit, kepstat

def kepoutlier(infile, outfile, datacol, nsig, stepsize, npoly, niter,
               operation, ranges, plot, plotfit, clobber, verbose, logfile):

    # startup parameters
    labelsize = 24
    ticksize = 16
    xsize = 16
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

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
            + ' clobber={}'.format(clobber)
            + ' verbose={}'.format(chatter)
            + ' logfile={}'.format(logfile)
    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPOUTLIER started at', logfile, verbose)
    # clobber output file
    if clobber:
        kepio.clobber(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPOUTLIER: {} exists. Use clobber=True'
                  .format(outfile))
        kepmsg.err(logfile, message, verbose)

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
        nanclean = instr[1].header['NANCLEAN']
    except:
        naxis2 = 0
        try:
            for i in range(len(table.field(0))):
                if (np.isfinite(table.field('barytime')[i])
                    and np.isfinite(table.field(datacol)[i])):
                    table[naxis2] = table[i]
                    naxis2 += 1
                    instr[1].data = table[:naxis2]
        except:
            for i in range(len(table.field(0))):
                if (np.isfinite(table.field('time')[i])
                    and np.isfinite(table.field(datacol)[i])):
                    table[naxis2] = table[i]
                    naxis2 += 1
                    instr[1].data = table[:naxis2]
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
        plt.figure(figsize=[xsize,ysize])
        plt.clf()

        # plot data
        ax = plt.axes([0.06, 0.1, 0.93, 0.87])

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=12)
        plt.plot(ptime, pout, color=lcolor, linestyle='-', linewidth=lwidth)
        plt.fill(ptime, pout, color=fcolor, linewidth=0.0, alpha=falpha)
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
                         'g', lw='3')
        except:
            for j in range(cstep1[i], cstep2[i] + 1):
                masterfit[j] = indata[j]
                mastersigma[j] = 1.0e10
            message = ('WARNING -- KEPOUTLIER: could not fit range '
                       + str(intime[cstep1[i]]) + '-' + str(intime[cstep2[i]]))
            kepmsg.warn(None, message)

    # reject outliers
    rejtime, rejdata = [], []
    naxis2 = 0
    for i in range(len(masterfit)):
        if abs(indata[i] - masterfit[i]) > nsig * mastersigma[i]
           and i in cadencelis:
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
    rejtime = np.array(rejtime, dtype='float64')
    rejdata = np.array(rejdata, dtype='float32')
    plt.plot(rejtime - intime0, rejdata / 10**nrm, 'ro')

    # plot ranges
    plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
    if ymin >= 0.0:
        plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
    else:
        plt.ylim(1.0e-10, ymax + yr * 0.01)

    # render plot
    plt.show()

    # write output file
    instr.writeto(outfile)
    # close input file
    kepio.closefits(instr, logfile, verbose)

    kepmsg.clock('KEPOUTLIER completed at', logfile, verbose)

# main
def kepoutiler_main():
    import argparse
    parser = argparse.ArgumentParser(
            description='Remove or replace data outliers from a time series')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output',
                        type=str)
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
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepoutlier.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepoutlier(args.infile, args.outfile, args.datacol, args.nsig,
               args.stepsize, args.npoly,args.niter, args.operation,
               args.ranges, args.plot, args.plotfit, args.clobber,
               args.verbose, args.logfile)
