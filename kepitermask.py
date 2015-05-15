import pylab, numpy, pyfits
from pylab import *
from matplotlib import *
from numpy import *
from pyfits import *
from scipy.stats import *
import kepio, kepmsg, kepkey, kepplot, kepstat, kepfunc, kepfit
import sys, time, re, math

# -----------------------------------------------------------
# core code

def kepitermask(infile,outfile,plotfile,column,row,timescale,nsig,stepsize,winsize,npoly,niter,
                clobber,verbose,logfile,status,cmdLine=False): 

# input arguments

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPITERMASK -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'plotfile='+plotfile+' '
    call += 'column='+str(column)+' '
    call += 'row='+str(row)+' '
    call += 'timescale='+str(timescale)+' '
    call += 'nsig='+str(nsig)+' '
    call += 'stepsize='+str(stepsize)+' '
    call += 'winsize='+str(winsize)+' '
    call += 'npoly='+str(npoly)+' '
    call += 'niter='+str(niter)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPITERMASK started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPITERMASK: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# open TPF FITS file

    if status == 0:
        try:
            kepid, channel, skygroup, module, output, quarter, season, \
                ra, dec, co, ro, kepmag, xdim, ydim, work1, status = \
                kepio.readTPF(infile,'TIME',logfile,verbose)
        except:
            message = 'ERROR -- KEPPRF: is %s a Target Pixel File? ' % infile
            status = kepmsg.err(logfile,message,verbose)

# print target data

    if status == 0:
        print ''
        print '      KepID:  %s' % kepid
        print ' RA (J2000):  %s' % ra
        print 'Dec (J2000): %s' % dec
        print '     KepMag:  %s' % kepmag
        print '   SkyGroup:    %2s' % skygroup
        print '     Season:    %2s' % str(season)
        print '    Channel:    %2s' % channel
        print '     Module:    %2s' % module
        print '     Output:     %1s' % output
        print ''

# read mask defintion data from TPF file

    if status == 0:
        maskmap, pixcoordx, pixcoordy, status = kepio.readMaskDefinition(infile,logfile,verbose)
        pixcoordx = rot90(pixcoordx)
        pixcoordy = flipud(rot90(pixcoordy))
        maskmap[:,:] = 0.0

# which pixel does the target reside on?

    if status == 0:
        x = where(pixcoordx == float(column))[1][0]
        y = where(pixcoordy == float(row))[0][0]
        maskmap[y,x] = 1.0

# read time series data

    if status == 0:
        instr = pyfits.open(infile,mode='readonly',memmap=True)
        work1 = instr[1].data.field('TIME')[:]
        work2 = instr[1].data.field('FLUX')[:]
        work3 = instr[1].data.field('QUALITY')[:]

# how many quality = 0 rows?

    if status == 0:
        npts = 0
        nrows = len(work1)
        for i in range(nrows):
            if work3[i] == 0 and numpy.isfinite(work1[i]):
                npts += 1
        time = empty((npts))
        flux = empty((npts,ydim,xdim))
        quality = empty((npts))

# construct pixel light curves from quality = 0 data

    if status == 0:
        n = 0
        for i in range(nrows):
            if work3[i] == 0 and numpy.isfinite(work1[i]):
                time[n] = work1[i]
                flux[n] = work2[i,:,:]
                quality[n] = work3[i]
                n +=1

# light curves from central pixel

    if status == 0:
        (pr, pc) = where(maskmap == 1.0)
        best_lc = flux[:,pr[0],pc[0]]

# calculate median CDPP

    if status == 0:
        best_median_cdpp, best_cdpp, status = \
            GetCDPP(time,best_lc,npoly,nsig,niter,winsize,stepsize,
                    timescale,logfile,verbose,status)

# does another pixel improve CDPP of the target?

    if status == 0:
        trial_med = best_median_cdpp
#        while best_median_cdpp == trial_med:
        for i in range(70):
            trial_lc, trial_cdpp, trial_med, xpix, ypix, status = \
                AddPixelToAperture(time,flux,maskmap,best_lc,npoly,nsig,niter,winsize,
                                   stepsize,timescale,logfile,verbose)
#            if trial_med < best_median_cdpp:
            if trial_med < 1e10:
                best_lc = trial_lc
                best_cdpp = trial_cdpp
                best_median_cdpp = trial_med
                maskmap[ypix,xpix] = 1.0
            print maskmap
            print i, best_median_cdpp

# plot style

    if status == 0:
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.0,
                      'axes.labelsize': 32,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 8,
                      'legend.fontsize': 8,
                      'xtick.labelsize': 12,
                      'ytick.labelsize': 12}
            pylab.rcParams.update(params)
        except:
            pass

# tmp

        pylab.plot(time,best_lc,color='#0000ff',linestyle='-',linewidth=1.0)        

# render plot

        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()	
        if plotfile.lower() != 'none':
            pylab.savefig(plotfile)

# stop time

    if status == 0:
        kepmsg.clock('KEPITERMASK ended at',logfile,verbose)

    return


# -----------------------------------------------------------
# calculate median and time-series CDPP

def GetCDPP(time,trial_lc,npoly,nsig,niter,winsize,stepsize,timescale,logfile,verbose,status):

# detrend data: find limits of each time step

    if status == 0:
        npts = len(time)
        tstep1 = []; tstep2 = []
        work = time[0]
        while work <= time[-1]:
            tstep1.append(work)
            tstep2.append(array([work+winsize,time[-1]],dtype='float64').min())
            work += stepsize

# detrend data: find cadence limits of each time step

    if status == 0:
        cstep1 = []; cstep2 = []
        for n in range(len(tstep1)):
            for i in range(len(time)-1):
                if time[i] <= tstep1[n] and time[i+1] > tstep1[n]:
                    for j in range(i,len(time)-1):
                        if time[j] < tstep2[n] and time[j+1] >= tstep2[n]:
                            cstep1.append(i)
                            cstep2.append(j+1)

# detrend data: loop over each time step, fit data, determine rms

    if status == 0:
        fitarray = zeros((npts,len(cstep1)),dtype='float32')
        fitarray[:,:] = numpy.nan
        masterfit = trial_lc * 0.0
        functype = 'poly' + str(npoly)
        for i in range(len(cstep1)):
            timeSeries = time[cstep1[i]:cstep2[i]+1]-time[cstep1[i]]
            dataSeries = trial_lc[cstep1[i]:cstep2[i]+1]
            pinit = [dataSeries.mean()]
            if npoly > 0:
                for j in range(npoly):
                    pinit.append(0.0)
            pinit = array(pinit,dtype='float32')
            try:
                coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
                    kepfit.lsqclip(functype,pinit,timeSeries,dataSeries,None,nsig,nsig,niter,
                                   logfile,verbose)
                fitarray[cstep1[i]:cstep2[i]+1,i] = 0.0
                for j in range(len(coeffs)):
                    fitarray[cstep1[i]:cstep2[i]+1,i] += coeffs[j] * timeSeries**j
            except:
                for j in range(cstep1[i],cstep2[i]+1):
                    fitarray[cstep1[i]:cstep2[i]+1,i] = 0.0
#                message  = 'WARNING -- KEPFLATTEN: could not fit range '
#                message += str(time[cstep1[i]]) + '-' + str(time[cstep2[i]])
#                kepmsg.warn(None,message)

# detrend data: find mean fit for each timestamp

    if status == 0:
        for i in range(npts):
            masterfit[i] = nanmean(fitarray[i,:])
        masterfit[-1] = masterfit[-4] #fudge
        masterfit[-2] = masterfit[-4] #fudge
        masterfit[-3] = masterfit[-4] #fudge

# detrend data: normalize light curve

    if status == 0:
        trial_lc = trial_lc / masterfit

# calculate STDDEV in units of ppm

    if status == 0:
        stddev = kepstat.running_frac_std(time,trial_lc,timescale/24) * 1.0e6

# calculate median STDDEV

    if status == 0:
        medstddev = ones((len(stddev)),dtype='float32') * median(stddev)
#        print '\nMedian %.1fhr STDDEV = %d ppm' % (timescale, median(stddev))

    return median(stddev), stddev, status


# -----------------------------------------------------------
# What is the best CDPP if I add a contiguous pixel to the aperture?

def AddPixelToAperture(time,flux,maskmap,lc,npoly,nsig,niter,winsize,stepsize,timescale,logfile,verbose):

    npts = len(time)
    work = zeros(shape(maskmap))

# where are the current aperture's next-door neighbors?

    (pr, pc) = where(maskmap == 1.0)
    for i in range(len(pc)):
        try:
            if pr[i]-1 >= 0:
                work[pr[i]-1,pc[i]] = 1.0
        except:
            pass
        try:
            work[pr[i]+1,pc[i]] = 1.0
        except:
            pass
        try:
            if pc[i]-1 >= 0:
                work[pr[i],pc[i]-1] = 1.0
        except:
            pass
        try:
            work[pr[i],pc[i]+1] = 1.0
        except:
            pass

    print shape(maskmap)[0]
    print shape(maskmap)[1]
    for j in range(shape(maskmap)[0]):
        for i in range(shape(maskmap)[1]):
            try:
                if pr[j]-1 >= 0 and pc[i]-1 >= 0:
                    if work[pr[j]-1,pc[i]] == 1.0 and \
                            work[pr[j]+1,pc[i]] == 1.0 and \
                            work[pr[j],pc[i]-1] == 1.0 and \
                            work[pr[j],pc[i]+1] == 1.0:
                        work[pr[j],pc[i]] = 1.0
            except:
                pass

    work = work - maskmap
    xpix = where(work == 1.0)[1]
    ypix = where(work == 1.0)[0]
    trial_lcs = zeros((len(xpix),npts))
    trial_cdpps = zeros((len(xpix),npts))
    trial_meds = zeros((len(xpix)))

# what's the best CDPP after adding one new pixel?

    for i in range(len(xpix)):
        trial_lcs[i,:] = lc + flux[:,ypix[i],xpix[i]]
        trial_meds[i], trial_cdpps[i,:], status = GetCDPP(time,trial_lcs[i,:],npoly,nsig,niter,winsize,
                                                              stepsize,timescale,logfile,verbose,0)
    n = where(trial_meds == min(trial_meds))[0][0]
        
    return trial_lcs[n,:], trial_cdpps[n,:], trial_meds[n], xpix[n], ypix[n], status


# -----------------------------------------------------------
# main

if '--shell' in sys.argv:

    import argparse
    
    parser = argparse.ArgumentParser(description='Iteratively identify an optimal pixel aperture')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)
    parser.add_argument('--plotfile', default='None', help='name of output PNG plot file', type=str)
    parser.add_argument('--column', help='CCD column number of the target', dest='column', type=int)
    parser.add_argument('--row', help='CCD row number of the target', dest='row', type=int)
    parser.add_argument('--timescale', '-t', default=6.5, help='STDDEV timescale', dest='timescale', type=float)
    parser.add_argument('--nsig', default=3., help='Sigma clipping threshold for outliers', type=float)
    parser.add_argument('--stepsize', default=0.5, help='Stepsize on which to fit data [days]', type=float)
    parser.add_argument('--winsize', default=5.0, help='Window size of data to fit after each step (>= stepsize) [days]', type=float)
    parser.add_argument('--npoly', default=3, help='Polynomial order for each fit', type=int)
    parser.add_argument('--niter', default=1, help='Maximum number of clipping iterations', type=int)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()    
    cmdLine=True
    kepitermask(args.infile,args.outfile,args.plotfile,args.column,args.row,args.timescale,
                args.nsig,args.stepsize,args.winsize,args.npoly,args.niter,args.clobber,
                args.verbose,args.logfile,args.status,cmdLine)
else:

    from pyraf import iraf

    parfile = iraf.osfn("kepler$kepitermask.par")
    t = iraf.IrafTaskFactory(taskname="kepitermask", value=parfile, function=kepitermask)
