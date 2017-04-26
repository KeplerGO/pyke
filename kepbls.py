import sys, time, math, re
import numpy as np
from copy import copy
from scipy import stats
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import kepio, kepmsg, kepkey

def kepbls(infile,outfile,datacol,errcol,minper,maxper,mindur,maxdur,nsearch,
           nbins,plot,clobber,verbose,logfile,status,cmdLine=False):

# startup parameters

    np.seterr(all="ignore")
    status = 0
    labelsize = 32
    ticksize = 18
    xsize = 16
    ysize = 8
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPBLS -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'errcol='+str(errcol)+' '
    call += 'minper='+str(minper)+' '
    call += 'maxper='+str(maxper)+' '
    call += 'mindur='+str(mindur)+' '
    call += 'maxdur='+str(maxdur)+' '
    call += 'nsearch='+str(nsearch)+' '
    call += 'nbins='+str(nbins)+' '
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

    kepmsg.clock('KEPBLS started at',logfile,verbose)

# is duration greater than one bin in the phased light curve?

    if nbins * maxdur / 24.0 / maxper <= 1.0:
        message = ('WARNING -- KEPBLS: ' + str(maxdur) +
                   ' hours transit duration < 1 phase bin when P = ')
        message += str(maxper) + ' days'
        kepmsg.warn(logfile,message)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPBLS: ' + outfile + ' exists. Use clobber=yes'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# read table structure

    if status == 0:
        table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# filter input data table

    if status == 0:
        work1 = np.array([table.field('time'), table.field(datacol), table.field(errcol)])
        work1 = np.rot90(work1,3)
        work1 = work1[~np.isnan(work1).any(1)]

# read table columns

    if status == 0:
        intime = work1[:,2] + bjdref
        indata = work1[:,1]
        inerr = work1[:,0]

# test whether the period range is sensible

    if status == 0:
        tr = intime[-1] - intime[0]
        if maxper > tr:
            message = 'ERROR -- KEPBLS: maxper is larger than the time range of the input data'
            status = kepmsg.err(logfile,message,verbose)

# prepare time series

    if status == 0:
        work1 = intime - intime[0]
        work2 = indata - np.mean(indata)

# start period search

    if status == 0:
        srMax = np.array([],dtype='float32')
        transitDuration = np.array([],dtype='float32')
        transitPhase = np.array([],dtype='float32')
        dPeriod = (maxper - minper) / nsearch
        trialPeriods = np.arange(minper,maxper+dPeriod,dPeriod,dtype='float32')
        complete = 0
        print ' '
        for trialPeriod in trialPeriods:
            fracComplete = float(complete) / float(len(trialPeriods) - 1) * 100.0
            txt  = '\r'
            txt += 'Trial period = '
            txt += str(int(trialPeriod))
            txt += ' days ['
            txt += str(int(fracComplete))
            txt += '% complete]'
            txt += ' ' * 20
            sys.stdout.write(txt)
            sys.stdout.flush()
            complete += 1
            srMax = np.append(srMax,0.0)
            transitDuration = np.append(transitDuration,np.nan)
            transitPhase = np.append(transitPhase,np.nan)
            trialFrequency = 1.0 / trialPeriod

# minimum and maximum transit durations in quantized phase units

            duration1 = max(int(float(nbins) * mindur / 24.0 / trialPeriod),2)
            duration2 = max(int(float(nbins) * maxdur / 24.0 / trialPeriod) + 1,duration1 + 1)

# 30 minutes in quantized phase units

            halfHour = int(0.02083333 / trialPeriod * nbins + 1)

# compute folded time series with trial period

            work4 = np.zeros((nbins),dtype='float32')
            work5 = np.zeros((nbins),dtype='float32')
            phase = np.array(((work1 * trialFrequency) - np.floor(work1 * trialFrequency)) * float(nbins),dtype='int')
            ptuple = np.array([phase, work2, inerr])
            ptuple = np.rot90(ptuple,3)
            phsort = np.array(sorted(ptuple,key=lambda ph: ph[2]))
            for i in range(nbins):
                elements = np.nonzero(phsort[:,2] == float(i))[0]
                work4[i] = np.mean(phsort[elements,1])
                work5[i] = math.sqrt(np.sum(np.power(phsort[elements,0], 2)) / len(elements))

# extend the work arrays beyond nbins by wrapping

            work4 = np.append(work4,work4[:duration2])
            work5 = np.append(work5,work5[:duration2])

# calculate weights of folded light curve points

            sigmaSum = np.nansum(np.power(work5,-2))
            omega = np.power(work5,-2) / sigmaSum

# calculate weighted phased light curve

            s = omega * work4

# iterate through trial period phase

            for i1 in range(nbins):

# iterate through transit durations

                for duration in range(duration1,duration2+1,int(halfHour)):

# calculate maximum signal residue

                    i2 = i1 + duration
                    sr1 = np.sum(np.power(s[i1:i2],2))
                    sr2 = np.sum(omega[i1:i2])
                    sr = math.sqrt(sr1 / (sr2 * (1.0 - sr2)))
                    if sr > srMax[-1]:
                        srMax[-1] = sr
                        transitDuration[-1] = float(duration)
                        transitPhase[-1] = float((i1 + i2) / 2)

# normalize maximum signal residue curve

        bestSr = np.max(srMax)
        bestTrial = np.nonzero(srMax == bestSr)[0][0]
        srMax /= bestSr
        transitDuration *= trialPeriods / 24.0
        BJD0 = np.array(transitPhase * trialPeriods / nbins,dtype='float64') + intime[0] - 2454833.0
        print '\n'

# clean up x-axis unit

    if status == 0:
        ptime = copy(trialPeriods)
        xlab = 'Trial Period (days)'

# clean up y-axis units

    if status == 0:
        pout = copy(srMax)
        ylab = 'Normalized Signal Residue'

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

        ax = plt.axes([0.06,0.10,0.93,0.87])

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90)

# plot curve

    if status == 0 and plot:
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
            plt.ion()
            plt.show()

# append new BLS data extension to the output file

    if status == 0:
        col1 = pyfits.Column(name='PERIOD',format='E',unit='days',array=trialPeriods)
        col2 = pyfits.Column(name='BJD0',format='D',unit='BJD - 2454833',array=BJD0)
        col3 = pyfits.Column(name='DURATION',format='E',unit='hours',array=transitDuration)
        col4 = pyfits.Column(name='SIG_RES',format='E',array=srMax)
        cols = pyfits.ColDefs([col1,col2,col3,col4])
        instr.append(pyfits.BinTableHDU.from_columns(cols))
        instr[-1].header.cards['TTYPE1'].comment = 'column title: trial period'
        instr[-1].header.cards['TTYPE2'].comment = 'column title: trial mid-transit zero-point'
        instr[-1].header.cards['TTYPE3'].comment = 'column title: trial transit duration'
        instr[-1].header.cards['TTYPE4'].comment = 'column title: normalized signal residue'
        instr[-1].header.cards['TFORM1'].comment = 'column type: float32'
        instr[-1].header.cards['TFORM2'].comment = 'column type: float64'
        instr[-1].header.cards['TFORM3'].comment = 'column type: float32'
        instr[-1].header.cards['TFORM4'].comment = 'column type: float32'
        instr[-1].header.cards['TUNIT1'].comment = 'column units: days'
        instr[-1].header.cards['TUNIT2'].comment = 'column units: BJD - 2454833'
        instr[-1].header.cards['TUNIT3'].comment = 'column units: hours'
        instr[-1].header['EXTNAME' ] = ('BLS','extension name')
        instr[-1].header['PERIOD'  ] = (trialPeriods[bestTrial],'most significant trial period [d]')
        instr[-1].header['BJD0'    ] = (BJD0[bestTrial] + 2454833.0,'time of mid-transit [BJD]')
        instr[-1].header['TRANSDUR'] = (transitDuration[bestTrial],'transit duration [hours]')
        instr[-1].header['SIGNRES' ] = (srMax[bestTrial] * bestSr,'maximum signal residue')

# history keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)
        instr.writeto(outfile)

# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)


# print best trial period results

    if status == 0:
        print '      Best trial period = %.5f days' % trialPeriods[bestTrial]
        print '    Time of mid-transit = BJD %.5f' % (BJD0[bestTrial] + 2454833.0)
        print '       Transit duration = %.5f hours' % transitDuration[bestTrial]
        print ' Maximum signal residue = %.4g \n' % (srMax[bestTrial] * bestSr)

# end time

    if status == 0:
        message = 'KEPBLS completed at'
    else:
        message = '\nKEPBLS aborted at'
    kepmsg.clock(message,logfile,verbose)

# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Remove or replace data outliers from a time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)
    parser.add_argument('--datacol', default='DETSAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--errcol', default='DETSAP_FLUX_ERR', help='Name of data error column to plot', type=str)
    parser.add_argument('--minper', default=1.0, help='Minimum search period [days]', type=float)
    parser.add_argument('--maxper', default=30.0, help='Maximum search period [days]', type=float)
    parser.add_argument('--mindur', default=0.5, help='Minimum transit duration [hours]', type=float)
    parser.add_argument('--maxdur', default=12.0, help='Maximum transit duration [hours]', type=float)
    parser.add_argument('--nsearch', default=1000, help='Number of test periods between minper and maxper', type=int)
    parser.add_argument('--nbins', default=1000, help='Number of bins in the folded time series at any test period', type=int)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    kepbls(args.infile,args.outfile,args.datacol,args.errcol,args.minper,args.maxper,args.mindur,
           args.maxdur,args.nsearch,args.nbins,args.plot,args.clobber,args.verbose,args.logfile,
           args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepbls.par")
    t = iraf.IrafTaskFactory(taskname="kepbls", value=parfile, function=kepbls)
