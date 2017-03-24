import numpy, sys, time, pylab, math, re
from astropy.io import fits as pyfits
from pylab import *
from matplotlib import *
from math import *
import kepio, kepmsg, kepkey, kepfit, kepstat
from numpy import logical_and

def kepdetrend(infile,outfile,datacol,errcol,ranges1,npoly1,nsig1,niter1,
               ranges2,npoly2,nsig2,niter2,popnans,plot,clobber,verbose,logfile,
               status,cmdLine=False):

# startup parameters

    status = 0
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
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDETREND -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'errcol='+str(errcol)+' '
    call += 'ranges1='+str(ranges1)+' '
    call += 'npoly1='+str(npoly1)+' '
    call += 'nsig1='+str(nsig1)+' '
    call += 'niter1='+str(niter1)+' '
    call += 'ranges2='+str(ranges2)+' '
    call += 'npoly2='+str(npoly2)+' '
    call += 'nsig2='+str(nsig2)+' '
    call += 'niter2='+str(niter2)+' '
    popn = 'n'
    if (popnans): popn = 'y'
    call += 'popnans='+popn+ ' '
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

    kepmsg.clock('KEPDETREND started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPDETREND: ' + outfile + ' exists. Use clobber=yes'
	    status = kepmsg.err(logfile,message,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# read table structure

    if status == 0:
	table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# filter input data table

    if status == 0:
        work1 = numpy.array([table.field('time'), table.field(datacol), table.field(errcol)])
        work1 = numpy.rot90(work1,3)
        work1 = work1[~numpy.isnan(work1).any(1)]            
 
# read table columns

    if status == 0:
        intime = work1[:,2] + bjdref
        indata = work1[:,1]
        inerr = work1[:,0]
        print intime

# time ranges for region 1 (region to be corrected)

    if status == 0:
        time1 = []; data1 = []; err1 = []
        t1start, t1stop, status = kepio.timeranges(ranges1,logfile,verbose)
    if status == 0:
        cadencelis1, status = kepstat.filterOnRange(intime,t1start,t1stop)
    if status == 0:
        for i in range(len(cadencelis1)):
            time1.append(intime[cadencelis1[i]])
            data1.append(indata[cadencelis1[i]])
            if errcol.lower() != 'none':
                err1.append(inerr[cadencelis1[i]])
        t0 = time1[0]
        time1 = array(time1,dtype='float64') - t0
        data1 = array(data1,dtype='float32')
        if errcol.lower() != 'none':
            err1 = array(err1,dtype='float32')
        else:
            err1 = None

# fit function to range 1

    if status == 0:
        functype = 'poly' + str(npoly1)
        pinit = [data1.mean()]
        if npoly1 > 0:
            for i in range(npoly1):
                pinit.append(0)
        pinit = array(pinit,dtype='float32')
        coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx1, ploty1, status = \
            kepfit.lsqclip(functype,pinit,time1,data1,err1,nsig1,nsig1,niter1,
                           logfile,verbose)
        fit1 = indata * 0.0
        for i in range(len(coeffs)):
            fit1 += coeffs[i] * (intime - t0)**i
        for i in range(len(intime)):
            if i not in cadencelis1:
                fit1[i] = 0.0
        plotx1 += t0
        print coeffs

# time ranges for region 2 (region that is correct)

    if status == 0:
        time2 = []; data2 = []; err2 = []
        t2start, t2stop, status = kepio.timeranges(ranges2,logfile,verbose)
        cadencelis2, status = kepstat.filterOnRange(intime,t2start,t2stop)
        for i in range(len(cadencelis2)):
            time2.append(intime[cadencelis2[i]])
            data2.append(indata[cadencelis2[i]])
            if errcol.lower() != 'none':
                err2.append(inerr[cadencelis2[i]])
        t0 = time2[0]
        time2 = array(time2,dtype='float64') - t0
        data2 = array(data2,dtype='float32')
        if errcol.lower() != 'none':
            err2 = array(err2,dtype='float32')
        else:
            err2 = None

# fit function to range 2

    if status == 0:
        functype = 'poly' + str(npoly2)
        pinit = [data2.mean()]
        if npoly2 > 0:
            for i in range(npoly2):
                pinit.append(0)
        pinit = array(pinit,dtype='float32')
        coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx2, ploty2, status = \
            kepfit.lsqclip(functype,pinit,time2,data2,err2,nsig2,nsig2,niter2,
                           logfile,verbose)
        fit2 = indata * 0.0
        for i in range(len(coeffs)):
            fit2 += coeffs[i] * (intime - t0)**i
        for i in range(len(intime)):
            if i not in cadencelis1:
                fit2[i] = 0.0
        plotx2 += t0

# normalize data

    if status == 0:
        outdata = indata - fit1 + fit2
        if errcol.lower() != 'none':
            outerr = inerr * 1.0

# comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
        if intime0 < 2.4e6: intime0 += 2.4e6
	ptime = intime - intime0
	plotx1 = plotx1 - intime0
	plotx2 = plotx2 - intime0
	xlab = 'BJD $-$ %d' % intime0

# clean up y-axis units

    if status == 0:
        pout = outdata
        ploty1
        ploty2
	nrm = len(str(int(numpy.nanmax(indata))))-1
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
        ptime = insert(ptime,[0],[ptime[0]]) 
        ptime = append(ptime,[ptime[-1]])
        indata = insert(indata,[0],[0.0]) 
        indata = append(indata,[0.0])
        pout = insert(pout,[0],[0.0]) 
        pout = append(pout,0.0)

# plot light curve

    if status == 0 and plot:
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
            pass

        pylab.figure(figsize=[xsize,ysize])
        pylab.clf()

# plot original data

        ax = pylab.axes([0.06,0.523,0.93,0.45])

# force tick labels to be absolute rather than relative

        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        pylab.setp(labels, 'rotation', 90, fontsize=12)

        pylab.plot(ptime,indata,color=lcolor,linestyle='-',linewidth=lwidth)
        pylab.fill(ptime,indata,color=fcolor,linewidth=0.0,alpha=falpha)
        pylab.plot(plotx1,ploty1,color='r',linestyle='-',linewidth=2.0)
        pylab.plot(plotx2,ploty2,color='g',linestyle='-',linewidth=2.0)
        pylab.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin > 0.0: 
            pylab.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            pylab.ylim(1.0e-10,ymax+yr*0.01)
	    pylab.ylabel(ylab, {'color' : 'k'})
        pylab.grid()

# plot detrended data

        ax = pylab.axes([0.06,0.073,0.93,0.45])

# force tick labels to be absolute rather than relative

        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        pylab.setp(labels, 'rotation', 90, fontsize=12)

        pylab.plot(ptime,pout,color=lcolor,linestyle='-',linewidth=lwidth)
        pylab.fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
        pylab.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin > 0.0: 
            pylab.ylim(omin-oo*0.01,omax+oo*0.01)
        else:
            pylab.ylim(1.0e-10,omax+oo*0.01)
	pylab.xlabel(xlab, {'color' : 'k'})
        try:
            pylab.ylabel(ylab, {'color' : 'k'})
        except:
            ylab = '10**%d e-/s' % nrm
            pylab.ylabel(ylab, {'color' : 'k'})

# render plot

    if status == 0:
        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
# write output file
    if status == 0 and popnans:
	    instr[1].data.field(datacol)[good_data] = outdata
	    instr[1].data.field(errcol)[good_data] = outerr
	    instr[1].data.field(datacol)[bad_data] = None
	    instr[1].data.field(errcol)[bad_data] = None
	    instr.writeto(outfile)
    elif status == 0 and not popnans:
        for i in range(len(outdata)):
            instr[1].data.field(datacol)[i] = outdata[i]
            if errcol.lower() != 'none':
                instr[1].data.field(errcol)[i] = outerr[i]
        instr.writeto(outfile)
    
# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

## end time

    if (status == 0):
	    message = 'KEPDETREND completed at'
    else:
	    message = '\nKEPDETREND aborted at'
    kepmsg.clock(message,logfile,verbose)

# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Detrend systematic features from Simple Aperture Photometry (SAP) data')
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
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()

    cmdLine=True
    
    kepdetrend(args.infile, args.outfile, args.datacol, args.errcol, args.ranges1, args.npoly1, args.nsig1, args.niter1,
               args.ranges2, args.npoly2, args.nsig2, args.niter2, args.popnans, args.plot, args.clobber, 
               args.verbose, args.logfile, args.status, cmdLine)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepdetrend.par")
    t = iraf.IrafTaskFactory(taskname="kepdetrend", value=parfile, function=kepdetrend)
