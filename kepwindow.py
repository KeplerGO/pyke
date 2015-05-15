
import numpy, sys, time, pyfits, pylab, math, re
from pyfits import *
from pylab import *
from matplotlib import *
from math import *
import kepio, kepmsg, kepkey, kepstat, kepfourier

def kepwindow(infile,outfile,fcol,fmax,nfreq,plot,clobber,verbose,logfile,status, cmdLine=False): 

## startup parameters

    status = 0
    labelsize = 24
    ticksize = 16
    xsize = 18
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

## log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPWINDOW -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'fcol='+fcol+' '
    call += 'fmax='+str(fmax)+' '
    call += 'nfreq='+str(nfreq)+' '
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

## start time

    kepmsg.clock('KEPWINDOW started at',logfile,verbose)

## test log file

    logfile = kepmsg.test(logfile)

## clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPWINDOW: ' + outfile + ' exists. Use clobber=yes'
        status = kepmsg.err(logfile,message,verbose)

## open input file

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

## fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

## read table columns

    if status == 0:
	try:
            barytime = instr[1].data.field('barytime')
	except:
            barytime, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
	signal, status = kepio.readfitscol(infile,instr[1].data,fcol,logfile,verbose)

## remove infinite data from time series

    if status == 0:
        incols = [barytime, signal]
        outcols = kepstat.removeinfinlc(signal, incols)
        barytime = outcols[0]
        signal = outcols[1]

## reset signal data to zero

    if status == 0:
        signal = ones(len(outcols[1]))

## frequency steps

    if status == 0:
        deltaf = fmax / nfreq

## loop through frequency steps; determine FT power

    if status == 0:
        fr, power = kepfourier.ft(barytime,signal,0.0,fmax,deltaf,True)
        power[0] = 1.0
        
## mirror window function around ordinate

    if status == 0:
        work1 = []; work2 = []
        for i in range(len(fr)-1, 0, -1):
            work1.append(-fr[i])
            work2.append(power[i])
        for i in range(len(fr)):
            work1.append(fr[i])
            work2.append(power[i])
        fr = array(work1,dtype='float32')
        power = array(work2,dtype='float32')

## write output file

    if status == 0:
        col1 = Column(name='FREQUENCY',format='E',unit='days',array=fr)
        col2 = Column(name='POWER',format='E',array=power)
        cols = ColDefs([col1,col2])
        instr.append(new_table(cols))
        instr[-1].header.update('EXTNAME','WINDOW FUNCTION','extension name')
        
## comment keyword in output file

    if status == 0:
        status = kepkey.comment(call,instr[0],outfile,logfile,verbose)
        instr.writeto(outfile)

## close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

## data limits

    if status == 0:
        nrm = len(str(int(power.max())))-1
        power = power / 10**nrm
        ylab = 'Power (x10$^%d$)' % nrm
	xmin = fr.min()
	xmax = fr.max()
	ymin = power.min()
	ymax = power.max()
	xr = xmax - xmin
	yr = ymax - ymin
        fr = insert(fr,[0],fr[0])
        fr = append(fr,fr[-1])
        power = insert(power,[0],0.0) 
        power = append(power,0.0)

## plot power spectrum

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
            print 'ERROR -- KEPWINDOW: install latex for scientific plotting'
            status = 1
    if status == 0 and plot:
        pylab.figure(1,figsize=[xsize,ysize])
        pylab.axes([0.06,0.113,0.93,0.86])
        pylab.plot(fr,power,color=lcolor,linestyle='-',linewidth=lwidth)
        fill(fr,power,color=fcolor,linewidth=0.0,alpha=falpha)
        xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin-yr*0.01 <= 0.0:
            ylim(1.0e-10,ymax+yr*0.01)
        else:
            ylim(ymin-yr*0.01,ymax+yr*0.01)
        xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
        ylabel('Power', {'color' : 'k'})

# render plot

        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
		
## end time

    if (status == 0):
	    message = 'KEPWINDOW completed at'
    else:
	    message = '\nKEPWINDOW aborted at'
    kepmsg.clock(message,logfile,verbose)

## main
if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate and store the window function for a Kepler time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)

    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column', type=str, dest='fcol')

    parser.add_argument('--fmax', default=1.0, help='Minimum search frequency [1/day]', type=float)
    parser.add_argument('--nfreq', default=100, help='Number of frequency intervals', type=int)


    parser.add_argument('--plot', action='store_true', help='Plot result?')

    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    cmdLine=True

    kepwindow(args.infile,args.outfile,args.fcol,args.fmax,args.nfreq,args.plot,args.clobber,args.verbose, 
        args.logfile,args.status, cmdLine)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepwindow.par")
    t = iraf.IrafTaskFactory(taskname="kepwindow", value=parfile, function=kepwindow)
