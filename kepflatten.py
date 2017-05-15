import numpy as np
import scipy, sys, time, math, re
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from copy import copy
import kepio, kepmsg, kepkey, kepfit, kepstat

def kepflatten(infile,outfile,datacol,errcol,nsig,stepsize,winsize,npoly,
               niter,ranges,plot,clobber,verbose,logfile,status,
               cmdLine=False):

# startup parameters

    status = 0
    labelsize = 32
    ticksize = 18
    xsize = 16
    ysize = 10
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFLATTEN -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'errcol='+str(errcol)+' '
    call += 'nsig='+str(nsig)+' '
    call += 'stepsize='+str(stepsize)+' '
    call += 'winsize='+str(winsize)+' '
    call += 'npoly='+str(npoly)+' '
    call += 'niter='+str(niter)+' '
    call += 'ranges='+str(ranges)+' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPFLATTEN started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# test winsize > stepsize

    if winsize < stepsize:
        message = 'ERROR -- KEPFLATTEN: winsize must be greater than stepsize'
        status = kepmsg.err(logfile,message,verbose)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPFLATTEN: ' + outfile + ' exists. Use clobber=yes'
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
            datac = table.field(datacol)
        except:
             message = 'ERROR -- KEPFLATTEN: cannot find or read data column ' + datacol
             status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        try:
            err = table.field(errcol)
        except:
             message = 'WARNING -- KEPFLATTEN: cannot find or read error column ' + errcol
             errcol = 'None'
    if status == 0:
        if errcol.lower() == 'none' or errcol == 'PSF_FLUX_ERR':
            err = datac * cadence
            err = np.sqrt(np.abs(err)) / cadence
            work1 = np.array([table.field('time'), datac, err])
        else:
            work1 = np.array([table.field('time'), datac, err])
        work1 = np.rot90(work1,3)
        work1 = work1[~np.isnan(work1).any(1)]

# read table columns

    if status == 0:
        intime = work1[:,2] + bjdref
        indata = work1[:,1]
        inerr = work1[:,0]
        if len(intime) == 0:
             message = 'ERROR -- KEPFLATTEN: one of the input arrays is all NaN'
             status = kepmsg.err(logfile,message,verbose)

# time ranges for region to be corrected

    if status == 0:
        t1, t2, status = kepio.timeranges(ranges,logfile,verbose)
        cadencelis = kepstat.filterOnRange(intime,t1,t2)

# find limits of each time step

    if status == 0:
        tstep1 = []; tstep2 = []
        work = intime[0]
        while work <= intime[-1]:
            tstep1.append(work)
            tstep2.append(np.array([work+winsize,intime[-1]],dtype='float64').min())
            work += stepsize

# find cadence limits of each time step

    if status == 0:
        cstep1 = []; cstep2 = []
        for n in range(len(tstep1)):
            for i in range(len(intime)-1):
                if intime[i] <= tstep1[n] and intime[i+1] > tstep1[n]:
                    for j in range(i,len(intime)-1):
                        if intime[j] < tstep2[n] and intime[j+1] >= tstep2[n]:
                            cstep1.append(i)
                            cstep2.append(j+1)

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
        pout = copy(indata)
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

        ax = plt.axes([0.06,0.54,0.93,0.43])

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90)
        plt.setp(plt.gca(),xticklabels=[])

        plt.plot(ptime[1:-1],pout[1:-1],color=lcolor,linestyle='-',linewidth=lwidth)
        plt.fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
        if not plotLatex:
            ylab = '10**%d electrons/sec' % nrm
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

# loop over each time step, fit data, determine rms

    if status == 0:
        fitarray = np.zeros((len(indata),len(cstep1)),dtype='float32')
        sigarray = np.zeros((len(indata),len(cstep1)),dtype='float32')
        fitarray[:,:] = np.nan
        sigarray[:,:] = np.nan
        masterfit = indata * 0.0
        mastersigma = np.zeros(len(masterfit))
        functype = 'poly' + str(npoly)
        for i in range(len(cstep1)):
            timeSeries = intime[cstep1[i]:cstep2[i]+1]-intime[cstep1[i]]
            dataSeries = indata[cstep1[i]:cstep2[i]+1]
            fitTimeSeries = np.array([],dtype='float32')
            fitDataSeries = np.array([],dtype='float32')
            pinit = [dataSeries.mean()]
            if npoly > 0:
                for j in range(npoly):
                    pinit.append(0.0)
            pinit = np.array(pinit,dtype='float32')
            try:
                if len(fitarray[cstep1[i]:cstep2[i]+1,i]) > len(pinit):
                    coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
                        kepfit.lsqclip(functype,pinit,timeSeries,dataSeries,None,nsig,nsig,niter,
                                       logfile,verbose)
                    fitarray[cstep1[i]:cstep2[i]+1,i] = 0.0
                    sigarray[cstep1[i]:cstep2[i]+1,i] = sigma
                    for j in range(len(coeffs)):
                        fitarray[cstep1[i]:cstep2[i]+1,i] += coeffs[j] * timeSeries**j
            except:
                for j in range(cstep1[i],cstep2[i]+1):
                    fitarray[cstep1[i]:cstep2[i]+1,i] = 0.0
                    sigarray[cstep1[i]:cstep2[i]+1,i] = 1.0e-10
                message  = 'WARNING -- KEPFLATTEN: could not fit range '
                message += str(intime[cstep1[i]]) + '-' + str(intime[cstep2[i]])
                kepmsg.warn(None,message)

# find mean fit for each timestamp

    if status == 0:
        for i in range(len(indata)):
            masterfit[i] = np.nanmean(fitarray[i,:])
            mastersigma[i] = np.nanmean(sigarray[i,:])
        masterfit[-1] = masterfit[-4] #fudge
        masterfit[-2] = masterfit[-4] #fudge
        masterfit[-3] = masterfit[-4] #fudge
        plt.plot(intime-intime0, masterfit / 10**nrm,'g',lw='3')

# reject outliers

    if status == 0:
        rejtime = []; rejdata = []; naxis2 = 0
        for i in range(len(masterfit)):
            if abs(indata[i] - masterfit[i]) > nsig * mastersigma[i] and i in cadencelis:
                rejtime.append(intime[i])
                rejdata.append(indata[i])
        rejtime = np.array(rejtime,dtype='float64')
        rejdata = np.array(rejdata,dtype='float32')
        if plot:
            plt.plot(rejtime-intime0,rejdata / 10**nrm,'ro')

# new data for output file

    if status == 0:
        outdata = indata / masterfit
        outerr = inerr / masterfit

# plot ranges

    if status == 0 and plot:
        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10,ymax+yr*0.01)

# plot residual data

    if status == 0 and plot:
        ax = plt.axes([0.06,0.09,0.93,0.43])

# force tick labels to be absolute rather than relative

    if status == 0 and plot:
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90)

# clean up y-axis units

    if status == 0:
        pout = copy(outdata)
        ylab = 'Normalized Flux'

# data limits

    if status == 0 and plot:
        ymin = pout.min()
        ymax = pout.max()
        yr = ymax - ymin
        pout = np.insert(pout,[0],[0.0])
        pout = np.append(pout,0.0)

        plt.plot(ptime[1:-1],pout[1:-1],color=lcolor,linestyle='-',linewidth=lwidth)
        plt.fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

# plot ranges

    if status == 0 and plot:
        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10,ymax+yr*0.01)

# render plot

    if status == 0 and plot:
        plt.savefig(re.sub('.fits','.png',outfile))
        plt.ion()
        plt.show()

# add NaNs back into data

    if status == 0:
        n = 0
        work1 = np.array([],dtype='float32')
        work2 = np.array([],dtype='float32')
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)
        tn = table.field('time')
        dn = table.field(datacol)
        for i in range(len(table.field(0))):
            if np.isfinite(tn[i]) and np.isfinite(dn[i]) and np.isfinite(err[i]):
                try:
                    work1 = np.append(work1,outdata[n])
                    work2 = np.append(work2,outerr[n])
                    n += 1
                except:
                    pass
            else:
                work1 = np.append(work1,np.nan)
                work2 = np.append(work2,np.nan)

# history keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# write output file

        try:
            col1 = pyfits.Column(name='DETSAP_FLUX',format='E13.7',array=work1)
            col2 = pyfits.Column(name='DETSAP_FLUX_ERR',format='E13.7',array=work2)
            cols = instr[1].data.columns + col1 + col2
            instr[1] = pyfits.BinTableHDU.from_columns(cols,header=instr[1].header)
            instr.writeto(outfile)
        except ValueError:
            try:
                instr[1].data.field('DETSAP_FLUX')[:] = work1
                instr[1].data.field('DETSAP_FLUX_ERR')[:] = work2
                instr.writeto(outfile)
            except:
                message = 'ERROR -- KEPFLATTEN: cannot add DETSAP_FLUX data to FITS file'
                status = kepmsg.err(logfile,message,verbose)

# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)

## end time

    if status == 0:
        message = 'KEPFLATTEN completed at'
    else:
        message = '\nKEPFLATTEN aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Remove or replace data outliers from a time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)
    parser.add_argument('--datacol', default='PDCSAP_FLUX', help='Name of data column to plot', type=str)
    parser.add_argument('--errcol', default='PDCSAP_FLUX_ERR', help='Name of data error column to plot', type=str)
    parser.add_argument('--nsig', default=3., help='Sigma clipping threshold for outliers', type=float)
    parser.add_argument('--stepsize', default=0.5, help='Stepsize on which to fit data [days]', type=float)
    parser.add_argument('--winsize', default=5.0,
                        help='Window size of data to fit after each step (>= stepsize) [days]', type=float)
    parser.add_argument('--npoly', default=3, help='Polynomial order for each fit', type=int)
    parser.add_argument('--niter', default=1, help='Maximum number of clipping iterations', type=int)
    parser.add_argument('--ranges', default='0,0', help='Time ranges of regions to filter', type=str)
    parser.add_argument('--plot', action='store_true', help='Plot result?', default=False)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepflatten.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    kepflatten(args.infile,args.outfile,args.datacol,args.errcol,args.nsig,args.stepsize,
               args.winsize,args.npoly,args.niter,args.ranges,args.plot,args.clobber,
               args.verbose,args.logfile,args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepflatten.par")
    t = iraf.IrafTaskFactory(taskname="kepflatten", value=parfile, function=kepflatten)
