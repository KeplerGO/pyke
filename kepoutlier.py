import numpy as np
import sys, time, re
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from math import *
import kepio, kepmsg, kepkey, kepfit, kepstat

def kepoutlier(infile,outfile,datacol,nsig,stepsize,npoly,niter,
               operation,ranges,plot,plotfit,clobber,verbose,logfile,status, cmdLine=False):

# startup parameters

    status = 0
    labelsize = 24
    ticksize = 16
    xsize = 16
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPOUTLIER -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'nsig='+str(nsig)+' '
    call += 'stepsize='+str(stepsize)+' '
    call += 'npoly='+str(npoly)+' '
    call += 'niter='+str(niter)+' '
    call += 'operation='+str(operation)+' '
    call += 'ranges='+str(ranges)+' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    plotf = 'n'
    if (plotfit): plotf = 'y'
    call += 'plotfit='+plotf+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPOUTLIER started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
	message = 'ERROR -- KEPOUTLIER: ' + outfile + ' exists. Use clobber=yes'
	status = kepmsg.err(logfile,message,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# read table structure

    if status == 0:
	table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# filter input data table

    if status == 0:
        try:
            nanclean = instr[1].header['NANCLEAN']
        except:
            naxis2 = 0
            try:
                for i in range(len(table.field(0))):
                    if np.isfinite(table.field('barytime')[i]) and \
                            np.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
            except:
                for i in range(len(table.field(0))):
                    if np.isfinite(table.field('time')[i]) and \
                            np.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
            comment = 'NaN cadences removed from data'
            status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)

# read table columns

    if status == 0:
	try:
            intime = instr[1].data.field('barytime') + 2.4e6
	except:
            intime, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
	indata, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)
    if status == 0:
        intime = intime + bjdref
        indata = indata / cadenom

# time ranges for region to be corrected

    if status == 0:
        t1, t2, status = kepio.timeranges(ranges,logfile,verbose)
        cadencelis, status = kepstat.filterOnRange(intime,t1,t2)

# find limits of each time step

    if status == 0:
        tstep1 = []; tstep2 = []
        work = intime[0]
        while work < intime[-1]:
            tstep1.append(work)
            tstep2.append(np.array([work+stepsize,intime[-1]],dtype='float64').min())
            work += stepsize

# find cadence limits of each time step

    if status == 0:
        cstep1 = []; cstep2 = []
        work1 = 0; work2 = 0
        for i in range(len(intime)):
            if intime[i] >= intime[work1] and intime[i] < intime[work1] + stepsize:
                work2 = i
            else:
                cstep1.append(work1)
                cstep2.append(work2)
                work1 = i; work2 = i
        cstep1.append(work1)
        cstep2.append(work2)

        outdata = indata * 1.0

# comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
	ptime = intime - intime0
	xlab = 'BJD $-$ %d' % intime0

# clean up y-axis units

    if status == 0:
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
        ptime = np.insert(ptime,[0],[ptime[0]])
        ptime = np.append(ptime,[ptime[-1]])
        pout = np.insert(pout,[0],[0.0])
        pout = np.append(pout,0.0)

# plot light curve

    if status == 0 and plot:
        plotLatex = True
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': labelsize,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': ticksize,
                      'ytick.labelsize': ticksize}
            rcParams.update(params)
        except:
            plotLatex = False
    if status == 0 and plot:
        plt.figure(figsize=[xsize,ysize])
        plt.clf()

# plot data

        ax = plt.axes([0.06,0.1,0.93,0.87])

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=12)

        plt.plot(ptime,pout,color=lcolor,linestyle='-',linewidth=lwidth)
        plt.fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
	plt.xlabel(xlab, {'color' : 'k'})
        if not plotLatex:
            ylab = '10**%d electrons/sec' % nrm
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

# loop over each time step, fit data, determine rms

    if status == 0:
        masterfit = indata * 0.0
        mastersigma = np.zeros(len(masterfit))
        functype = 'poly' + str(npoly)
        for i in range(len(cstep1)):
            pinit = [indata[cstep1[i]:cstep2[i]+1].mean()]
            if npoly > 0:
                for j in range(npoly):
                    pinit.append(0.0)
            pinit = np.array(pinit,dtype='float32')
            try:
                coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
                    kepfit.lsqclip(functype,pinit,intime[cstep1[i]:cstep2[i]+1]-intime[cstep1[i]],
                                   indata[cstep1[i]:cstep2[i]+1],None,nsig,nsig,niter,logfile,
                                   verbose)
                for j in range(len(coeffs)):
                    masterfit[cstep1[i]:cstep2[i]+1] += coeffs[j] * \
                        (intime[cstep1[i]:cstep2[i]+1] - intime[cstep1[i]])**j
                for j in range(cstep1[i],cstep2[i]+1):
                    mastersigma[j] = sigma
                if plotfit:
                    plt.plot(plotx+intime[cstep1[i]]-intime0,ploty / 10**nrm,
                               'g',lw='3')
            except:
                for j in range(cstep1[i],cstep2[i]+1):
                    masterfit[j] = indata[j]
                    mastersigma[j] = 1.0e10
                message  = 'WARNING -- KEPOUTLIER: could not fit range '
                message += str(intime[cstep1[i]]) + '-' + str(intime[cstep2[i]])
                kepmsg.warn(None,message)

# reject outliers

    if status == 0:
        rejtime = []; rejdata = []; naxis2 = 0
        for i in range(len(masterfit)):
            if abs(indata[i] - masterfit[i]) > nsig * mastersigma[i] and i in cadencelis:
                rejtime.append(intime[i])
                rejdata.append(indata[i])
                if operation == 'replace':
                    [rnd] = kepstat.randarray([masterfit[i]],[mastersigma[i]])
                    table[naxis2] = table[i]
                    table.field(datacol)[naxis2] = rnd
                    naxis2 += 1
            else:
                table[naxis2] = table[i]
                naxis2 += 1
        instr[1].data = table[:naxis2]
        rejtime = np.array(rejtime,dtype='float64')
        rejdata = np.array(rejdata,dtype='float32')
        plt.plot(rejtime-intime0,rejdata / 10**nrm,'ro')

# plot ranges

        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10,ymax+yr*0.01)

# render plot
        plt.ion()
        plt.show()

# write output file

    if status == 0:
        instr.writeto(outfile)
# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)

# end time

    if (status == 0):
	    message = 'KEPOUTLIER completed at'
    else:
	    message = '\nKEPOUTLIER aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Remove or replace data outliers from a time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)

    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column to plot', type=str)

    parser.add_argument('--nsig', default=3., help='Sigma clipping threshold for outliers', type=float)
    parser.add_argument('--stepsize', default=1.0, help='Stepsize on which to fit data [days]', type=float)
    parser.add_argument('--npoly', default=3, help='Polynomial order for each fit', type=int)
    parser.add_argument('--niter', default=1, help='Maximum number of clipping iterations', type=int)

    parser.add_argument('--operation', default='remove', help='Remove or replace outliers?',
        type=str, choices=['replace','remove'])

    parser.add_argument('--ranges', default='0,0', help='Time ranges of regions to filter', type=str)

    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--plotfit', action='store_true', help='Plot fit over results?')

    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()

    cmdLine=True

    kepoutlier(args.infile,args.outfile,args.datacol,args.nsig,args.stepsize,args.npoly,args.niter,
               args.operation,args.ranges,args.plot,args.plotfit,args.clobber,args.verbose,args.logfile,args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepoutlier.par")
    t = iraf.IrafTaskFactory(taskname="kepoutlier", value=parfile, function=kepoutlier)
