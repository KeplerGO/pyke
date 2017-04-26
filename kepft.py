import sys, time, math, re
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey, kepstat, kepfourier

def kepft(infile,outfile,fcol,pmin,pmax,nfreq,plot,clobber,verbose,logfile,status, cmdLine=False):

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
    call = 'KEPFT -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'fcol='+fcol+' '
    call += 'pmin='+str(pmin)+' '
    call += 'pmax='+str(pmax)+' '
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

    kepmsg.clock('Start time is', logfile, verbose)

## test log file

    logfile = kepmsg.test(logfile)

## clobber output file

    if clobber: status = kepio.clobber(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPFT: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile, message, verbose)

## open input file

    if status == 0:
        instr, status = kepio.openfits(infile, 'readonly', logfile, verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr, infile,
                                                                logfile,
                                                                verbose,
                                                                status)

## fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

## read table columns

    if status == 0:
        try:
            barytime = instr[1].data.field('barytime')
        except:
            barytime, status = kepio.readfitscol(infile, instr[1].data,
                                                 'time', logfile, verbose)
        signal, status = kepio.readfitscol(infile, instr[1].data, fcol,
                                           logfile, verbose)
    if status == 0:
        barytime = barytime + bjdref

## remove infinite data from time series

    if status == 0:
        incols = [barytime, signal]
        outcols = kepstat.removeinfinlc(signal, incols)
        barytime = outcols[0]
        signal = outcols[1] - np.median(outcols[1])

## period to frequency conversion

    fmin = 1.0 / pmax
    fmax = 1.0 / pmin
    deltaf = (fmax - fmin) / nfreq

## loop through frequency steps; determine FT power

    if status == 0:
        fr, power = kepfourier.ft(barytime,signal,fmin,fmax,deltaf,True)

## write output file

    if status == 0:
        col1 = pyfits.Column(name='FREQUENCY',format='E',unit='1/day',array=fr)
        col2 = pyfits.Column(name='POWER',format='E',array=power)
        cols = pyfits.ColDefs([col1,col2])
        instr.append(pyfits.BinTableHDU.from_columns(cols))
        instr[-1].header['EXTNAME'] = ('POWER SPECTRUM','extension name')
        instr.writeto(outfile)
## history keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

## close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)

## data limits

    if status == 0:
        nrm = int(math.log10(power.max()))
        power = power / 10**nrm
        ylab = 'Power (x10$^{%d}$)' % nrm
        xmin = fr.min()
        xmax = fr.max()
        ymin = power.min()
        ymax = power.max()
        xr = xmax - xmin
        yr = ymax - ymin
        fr = np.insert(fr,[0],fr[0])
        fr = np.append(fr,fr[-1])
        power = np.insert(power,[0],0.0)
        power = np.append(power,0.0)

    if status == 0 and plot:
        plt.figure(1,figsize=[xsize,ysize])
        plt.clf()
        plt.axes([0.06,0.113,0.93,0.86])
        plt.plot(fr,power,color=lcolor,linestyle='-',linewidth=lwidth)
        plt.fill(fr,power,color=fcolor,linewidth=0.0,alpha=falpha)
        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin-yr*0.01 <= 0.0:
            plt.ylim(1.0e-10,ymax+yr*0.01)
        else:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        plt.xlabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

# render plot
        plt.ion()
        plt.show()
## end time

    if (status == 0):
        message = 'KEPFT completed at'
    else:
        message = '\nKEPFT aborted at'
    kepmsg.clock(message,logfile,verbose)

## main
if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description=('Calculate and store a '
                                                  'Fourier Transform from a '
                                                  'Kepler time series'))
    parser.add_argument('--shell', action='store_true',
                        help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output',
                        type=str)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str, dest='fcol')
    parser.add_argument('--pmin', default=0.1, help='Minimum search period [days]', type=float)
    parser.add_argument('--pmax', default=10.,
                        help='Maximum search period [days]', type=float)
    parser.add_argument('--nfreq', default=100,
                        help='Number of frequency intervals', type=int)
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    kepft(args.infile,args.outfile, args.fcol, args.pmin, args.pmax, args.nfreq,
        args.plot, args.clobber, args.verbose, args.logfile, args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepft.par")
    t = iraf.IrafTaskFactory(taskname="kepft", value=parfile, function=kepft)
