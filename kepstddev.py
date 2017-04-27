import sys
import math
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import kepio, kepmsg, kepstat, kepkey
from scipy import stats
from kepstat import savitzky_golay, running_frac_std
from copy import copy

def kepstddev(infile,outfile,datacol,timescale,clobber,verbose,logfile,status,cmdLine=False):

# startup parameters

    status = 0
    labelsize = 44
    ticksize = 36
    xsize = 16
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPSTDDEV -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'timescale='+str(timescale)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPSTDDEV started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPSTDDEV: ' + outfile + ' exists. Use clobber=yes'
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
        work1 = np.array([table.field('time'), table.field(datacol)])
        work1 = np.rot90(work1,3)
        work1 = work1[~np.isnan(work1).any(1)]

# read table columns

    if status == 0:
        intime = work1[:,1] + bjdref
        indata = work1[:,0]

# calculate STDDEV in units of ppm

    if status == 0:
        stddev = running_frac_std(intime,indata,timescale/24) * 1.0e6
        astddev = np.std(indata) * 1.0e6
        cdpp = stddev / math.sqrt(timescale * 3600.0 / cadence)
        print '\nStandard deviation = %d ppm' % astddev

# filter cdpp

    if status == 0:
        for i in range(len(cdpp)):
            if cdpp[i] > np.median(cdpp) * 10.0: cdpp[i] = cdpp[i-1]

# calculate median STDDEV

    if status == 0:
        med = np.median(cdpp[:])
        medcdpp = np.ones((len(cdpp)),dtype='float32') * med
        print 'Median %.1fhr CDPP = %d ppm' % (timescale, med)

# calculate RMS STDDEV

    if status == 0:
        rms, status = kepstat.rms(cdpp,np.zeros(len(stddev)),logfile,verbose)
        rmscdpp = np.ones((len(cdpp)),dtype='float32') * rms
        print '   RMS %.1fhr CDPP = %d ppm\n' % (timescale, rms)

# clean up x-axis unit

    if status == 0:
        intime0 = float(int(tstart / 100) * 100.0)
        ptime = intime - intime0
        xlab = 'BJD $-$ %d' % intime0

# clean up y-axis units

    if status == 0:
        pout = copy(cdpp)
        nrm = math.ceil(math.log10(np.median(cdpp))) - 1.0
        ylab = '%.1fhr $\sigma$ (ppm)' % timescale

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

# define size of plot on monitor screen

        plt.figure(figsize=[xsize,ysize])

# delete any fossil plots in the matplotlib window

        plt.clf()

# position first axes inside the plotting window

        ax = plt.axes([0.07,0.15,0.92,0.83])

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90,fontsize=36)

# plot flux vs time

        ltime = np.array([],dtype='float64')
        ldata = np.array([],dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for i in range(1,len(ptime)-1):
            dt = ptime[i] - ptime[i-1]
            if dt < work1:
                ltime = np.append(ltime,ptime[i])
                ldata = np.append(ldata,pout[i])
            else:
                plt.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)
                ltime = np.array([],dtype='float64')
                ldata = np.array([],dtype='float32')
        plt.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)

# plot the fill color below data time series, with no data gaps

        plt.fill(ptime,pout,fc='#ffff00',linewidth=0.0,alpha=0.2)

# define plot x and y limits

        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin - yr * 0.01 <= 0.0:
            plt.ylim(1.0e-10, ymax + yr * 0.01)
        else:
            plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
# plot labels

        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})

# make grid on plot

        plt.grid()

# render plot

    if status == 0:
        plt.ion()
        plt.show(block=True)

# add NaNs back into data

    if status == 0:
        n = 0
        work1 = np.array([],dtype='float32')
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)
        for i in range(len(table.field(0))):
            if isfinite(table.field('time')[i]) and isfinite(table.field(datacol)[i]):
                work1 = np.append(work1,cdpp[n])
                n += 1
            else:
                work1 = np.append(work1,nan)

# write output file
    if status == 0:
        status = kepkey.new('MCDPP%d' % (timescale * 10.0),medcdpp[0],
                            'Median %.1fhr CDPP (ppm)' % timescale,
                            instr[1],outfile,logfile,verbose)
        status = kepkey.new('RCDPP%d' % (timescale * 10.0),rmscdpp[0],
                            'RMS %.1fhr CDPP (ppm)' % timescale,
                            instr[1],outfile,logfile,verbose)
        colname = 'CDPP_%d' % (timescale * 10)
        col1 = pyfits.Column(name=colname,format='E13.7',array=work1)
        cols = instr[1].data.columns + col1
        instr[1] = pyfits.new_table(cols,header=instr[1].header)
        instr.writeto(outfile)
# comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# close FITS

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)

# end time

    if status == 0:
        message = 'KEPSTDDEV completed at'
    else:
        message = '\nKEPSTDDEV aborted at'
    kepmsg.clock(message,logfile,verbose)

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Calculate CDPP from a time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input FITS file', type=str)
    parser.add_argument('outfile', help='Name of output FITS file', type=str)
    parser.add_argument('--datacol', default='PDCSAP_FLUX', help='Name of data column to plot', type=str)
    parser.add_argument('--timescale', '-t', default=6.5, help='CDPP timescale', dest='timescale', type=float)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepstddev.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    kepstddev(args.infile,args.outfile,args.datacol,args.timescale,args.clobber,args.verbose,
           args.logfile,args.status,cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepstddev.par")
    t = iraf.IrafTaskFactory(taskname="kepstddev", value=parfile, function=kepstddev)
