import numpy as np
import sys
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import kepio, kepmsg, kepkey, kepstat

def kepclip(infile, outfile, ranges, plot, plotcol, clobber, verbose, logfile,
            status, cmdLine=False):

# startup parameters

    status = 0
    labelsize = 32
    ticksize = 24
    xsize = 18
    ysize = 10
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPCLIP -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'ranges='+ranges + ' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    call += 'plotcol='+plotcol+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPCLIP started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPCLIP: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# time ranges for region

    if status == 0:
        t1 = []; t2 = []
        t1, t2, status = kepio.timeranges(ranges,logfile,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
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

# input data

    if status == 0:
        table = instr[1].data

# read time and flux columns

    if status == 0:
        barytime, status = kepio.readtimecol(infile,table,logfile,verbose)
    if status == 0:
        flux, status = kepio.readfitscol(infile,table,plotcol,logfile,verbose)
    if status == 0:
        barytime = barytime + bjdref
        if 'flux' in plotcol.lower():
            flux = flux / cadenom

# filter input data table

    if status == 0:
        naxis2 = 0
        work1 = np.array([],'float64')
        work2 = np.array([],'float32')
        for i in range(len(barytime)):
            if (np.isfinite(barytime[i]) and np.isfinite(flux[i]) and flux[i] != 0.0):
                reject = False
                for j in range(len(t1)):
                    if (barytime[i] >= t1[j] and barytime[i] <= t2[j]):
                        reject = True
                if not reject:
                    table[naxis2] = table[i]
                    work1 = np.append(work1,barytime[i])
                    work2 = np.append(work2,flux[i])
                    naxis2 += 1

# comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# write output file

    if status == 0:
        instr[1].data = table[:naxis2]
        comment = 'NaN cadences removed from data'
        status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)
        instr.writeto(outfile)

# clean up x-axis unit

    if status == 0:
        barytime0 = float(int(tstart / 100) * 100.0)
        barytime = work1 - barytime0
        xlab = 'BJD $-$ %d' % barytime0

# clean up y-axis units

    if status == 0:
        try:
            nrm = len(str(int(work2.max())))-1
        except:
            nrm = 0
        flux = work2 / 10**nrm
        ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

# data limits

        xmin = barytime.min()
        xmax = barytime.max()
        ymin = flux.min()
        ymax = flux.max()
        xr = xmax - xmin
        yr = ymax - ymin

# clear window, plot box

    if status == 0 and plot:
        plt.figure(figsize=[xsize,ysize])
        plt.clf()
        ax = plt.axes([0.05,0.1,0.94,0.88])

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=12)

# plot line data

        ltime = [barytime[0]]; ldata = [flux[0]]
        for i in range(1,len(flux)):
            if (barytime[i-1] > barytime[i] - 0.025):
                ltime.append(barytime[i])
                ldata.append(flux[i])
            else:
                ltime = np.array(ltime, dtype=np.float)
                ldata = np.array(ldata, dtype=np.float)
                plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)
                ltime = []; ldata = []
        ltime = np.array(ltime, dtype=np.float)
        ldata = np.array(ldata, dtype=np.float)
        plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)

# plot fill data

        barytime = np.insert(barytime,[0],[barytime[0]])
        barytime = np.append(barytime,[barytime[-1]])
        flux = np.insert(flux,[0],[0.0])
        flux = np.append(flux,[0.0])
        plt.fill(barytime,flux,fc=fcolor,linewidth=0.0,alpha=falpha)
        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin-yr*0.01 <= 0.0:
            plt.ylim(1.0e-10,ymax+yr*0.01)
        else:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.grid()

# render plot

    if status == 0 and plot:
        plt.ion()
        plt.show()

# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)

# end time

    if (status == 0):
        message = 'KEPCLIP completed at'
    else:
        message = '\nKEPCLIP aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Remove unwanted time '
                                     'ranges from Kepler time series data')
    parser.add_argument('--shell', action='store_true',
                        help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output',
                        type=str)
    parser.add_argument('ranges',
                        help='List of time domain ranges to be excluded',
                        type=str)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--plotcol', '-p',help='Data column to plot',
                        default='SAP_FLUX', dest='plotcol', type=str)
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepclip.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine = True
    kepclip(args.infile, args.outfile, args.ranges, args.plot, args.plotcol,
            args.clobber, args.verbose, args.logfile, args.status,cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepclip.par")
    t = iraf.IrafTaskFactory(taskname="kepclip", value=parfile, function=kepclip)
