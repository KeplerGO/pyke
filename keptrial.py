
import numpy, sys, time, pyfits, pylab, math, re
from pyfits import *
from pylab import *
from matplotlib import *
import matplotlib.mlab as mlab
from math import *
from numpy import *
import kepio, kepmsg, kepkey, kepfit, kepfunc, kepstat, kepfourier

def keptrial(infile,outfile,datacol,errcol,fmin,fmax,nfreq,method,
             ntrials,plot,clobber,verbose,logfile,status,cmdLine=False): 

# startup parameters

    status = 0
    labelsize = 24
    ticksize = 16
    xsize = 18
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPTRIAL -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+datacol+' '
    call += 'errcol='+errcol+' '
    call += 'fmin='+str(fmin)+' '
    call += 'fmax='+str(fmax)+' '
    call += 'nfreq='+str(nfreq)+' '
    call += 'method='+method+' '
    call += 'ntrials='+str(ntrials)+' '
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

    kepmsg.clock('KEPTRIAL started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPTRIAL: ' + outfile + ' exists. Use clobber=yes'
	    kepmsg.err(logfile,message,verbose)
	    status = 1

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# input data

    if status == 0:
	try:
            barytime = instr[1].data.field('barytime')
	except:
            barytime, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
    if status == 0:
        signal, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)
    if status == 0:
        err, status = kepio.readfitscol(infile,instr[1].data,errcol,logfile,verbose)

# remove infinite data from time series

    if status == 0:
        try:
            nanclean = instr[1].header['NANCLEAN']
        except:
	    incols = [barytime, signal, err]
	    [barytime, signal, err] = kepstat.removeinfinlc(signal, incols)

# set up plot

    if status == 0:
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
            print 'WARNING: install latex for scientific plotting'
            plotLatex = False

# frequency steps and Monte Carlo iterations

    if status == 0:
        deltaf = (fmax - fmin) / nfreq
        freq = []; pmax = []; trial = []
        for i in range(ntrials):
            trial.append(i+1)

# adjust data within the error bars

            work1 = kepstat.randarray(signal,err)

# determine FT power
            fr, power = kepfourier.ft(barytime,work1,fmin,fmax,deltaf,False)

# determine peak in FT

            pmax.append(-1.0e30)
            for j in range(len(fr)):
                if (power[j] > pmax[-1]):
                    pmax[-1] = power[j]
                    f1 = fr[j]
            freq.append(f1)

# plot stop-motion histogram

            pylab.ion()
	    pylab.figure(1,figsize=[7,10])
            clf()
	    pylab.axes([0.08,0.08,0.88,0.89])
            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            n,bins,patches = pylab.hist(freq,bins=nfreq,range=[fmin,fmax],
                                        align='mid',rwidth=1,ec='#0000ff',
                                        fc='#ffff00',lw=2)

# fit normal distribution to histogram

            x = zeros(len(bins))
            for j in range(1,len(bins)):
                x[j] = (bins[j] + bins[j-1]) / 2
            pinit = numpy.array([float(i),freq[-1],deltaf])
            if i > 3:
                n = array(n,dtype='float32')
                coeffs, errors, covar, sigma, chi2, dof, fit, plotx, ploty, status = \
                    kepfit.leastsquare('gauss',pinit,x[1:],n,None,logfile,verbose)
                fitfunc = kepfunc.gauss()
                f = arange(fmin,fmax,(fmax-fmin)/100)
                fit = fitfunc(coeffs,f)
                pylab.plot(f,fit,'r-',linewidth=2)
            if plotLatex:
                xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
            else:
                xlabel(r'Frequency (1/d)', {'color' : 'k'})
            ylabel('N', {'color' : 'k'})
            xlim(fmin,fmax)
	    grid()

# render plot

        if plot:
            if cmdLine: 
                pylab.show()
            else: 
                pylab.ion()
                pylab.plot([])
                pylab.ioff()

# period results

    if status == 0:
        p = 1.0 / coeffs[1]
        perr = p * coeffs[2] / coeffs[1]
        f1 = fmin; f2 = fmax
        gotbin = False
        for i in range(len(n)):
            if n[i] > 0 and not gotbin:
                f1 = bins[i]
                gotbin = True
        gotbin = False
        for i in range(len(n)-1,0,-1):
            if n[i] > 0 and not gotbin:
                f2 = bins[i+1]
                gotbin = True
        powave, powstdev = kepstat.stdev(pmax)

# print result

    if status == 0:
        print '              best period: %.10f days (%.7f min)' % (p, p * 1440.0)
        print '     1-sigma period error: %.10f days (%.7f min)' % (perr, perr * 1440.0)
        print '             search range: %.10f - %.10f days  ' % (1.0 / fmax, 1.0 / fmin)
        print '    100%% confidence range: %.10f - %.10f days  ' % (1.0 / f2, 1.0 / f1)
#        print '     detection confidence: %.2f sigma' % (powave / powstdev)
        print '         number of trials: %d' % ntrials
        print ' number of frequency bins: %d' % nfreq

# history keyword in output file

    if status == 0:
	    status = kepkey.history(call,instr[0],outfile,logfile,verbose)

## write output file

    if status == 0:
        col1 = Column(name='TRIAL',format='J',array=trial)
        col2 = Column(name='FREQUENCY',format='E',unit='1/day',array=freq)
        col3 = Column(name='POWER',format='E',array=pmax)
        cols = ColDefs([col1,col2,col3])
        instr.append(new_table(cols))
        try:
            instr[-1].header.update('EXTNAME','TRIALS','Extension name')
        except:
            status = 1
        try:
            instr[-1].header.update('SEARCHR1',1.0 / fmax,'Search range lower bound (days)')
        except:
            status = 1
        try:
            instr[-1].header.update('SEARCHR2',1.0 / fmin,'Search range upper bound (days)')
        except:
            status = 1
        try:
            instr[-1].header.update('NFREQ',nfreq,'Number of frequency bins')
        except:
            status = 1
        try:
            instr[-1].header.update('PERIOD',p,'Best period (days)')
        except:
            status = 1
        try:
            instr[-1].header.update('PERIODE',perr,'1-sigma period error (days)')
        except:
            status = 1
#        instr[-1].header.update('DETNCONF',powave/powstdev,'Detection significance (sigma)')
        try:
            instr[-1].header.update('CONFIDR1',1.0 / f2,'Trial confidence lower bound (days)')
        except:
            status = 1
        try:
            instr[-1].header.update('CONFIDR2',1.0 / f1,'Trial confidence upper bound (days)')
        except:
            status = 1
        try:
            instr[-1].header.update('NTRIALS',ntrials,'Number of trials')
        except:
            status = 1
        instr.writeto(outfile)
    
# close input file

    if status == 0:
	    status = kepio.closefits(instr,logfile,verbose)	    

## end time

    if (status == 0):
	    message = 'KEPTRAIL completed at'
    else:
	    message = '\nKEPTRIAL aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate best period and error estimate from Fourier transform')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)

    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column', type=str)
    parser.add_argument('--errcol', default='SAP_FLUX_ERR', help='Name of data error column', type=str)

    parser.add_argument('--fmin', default=0.1, help='Minimum search frequency [1/day]', type=float)
    parser.add_argument('--fmax', default=50., help='Minimum search frequency [1/day]', type=float)
    parser.add_argument('--nfreq', default=100, help='Number of frequency intervals', type=int)

    parser.add_argument('--method', default='ft', help='Frequency search method', type=int, choices=['ft'])
    parser.add_argument('--ntrials', default=1000, help='Number of search trials', type=int)

    parser.add_argument('--plot', action='store_true', help='Plot result?')

    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    cmdLine=True

    keptrial(args.infile,args.outfile,args.datacol,args.errcol,args.fmin,args.fmax,args.nfreq,args.method,
             args.ntrials,args.plot,args.clobber,args.verbose,args.logfile,args.status, cmdLine)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$keptrial.par")
    t = iraf.IrafTaskFactory(taskname="keptrial", value=parfile, function=keptrial)
