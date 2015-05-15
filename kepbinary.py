from pyraf import iraf
import pylab, numpy
from pylab import *
from matplotlib import *
import scipy
from scipy import optimize
from scipy.optimize import fmin
import kepio, kepmsg, kepkey, kepsim
import re

niter = 0

def kepbinary(infile,outfile,datacol,m1,m2,r1,r2,period,bjd0,eccn,omega,inclination,
              c1,c2,c3,c4,albedo,depth,contamination,gamma,fitparams,eclipses,dopboost,
              tides,job,clobber,verbose,logfile,status): 

# startup parameters

    status = 0
    labelsize = 24; ticksize = 16; xsize = 17; ysize = 7
    lcolor = '#0000ff'; lwidth = 1.0; fcolor = '#ffff00'; falpha = 0.2

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPBINARY -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+datacol+' '
    call += 'm1='+str(m1)+' '
    call += 'm2='+str(m2)+' '
    call += 'r1='+str(r1)+' '
    call += 'r2='+str(r2)+' '
    call += 'period='+str(period)+' '
    call += 'bjd0='+str(bjd0)+' '
    call += 'eccn='+str(eccn)+' '
    call += 'omega='+str(omega)+' '
    call += 'inclination='+str(inclination)+' '
    call += 'c1='+str(c1)+' '
    call += 'c2='+str(c2)+' '
    call += 'c3='+str(c3)+' '
    call += 'c4='+str(c4)+' '
    call += 'albedo='+str(albedo)+' '
    call += 'depth='+str(depth)+' '
    call += 'contamination='+str(contamination)+' '
    call += 'gamma='+str(gamma)+' '
    call += 'fitparams='+str(fitparams)+' '
    eclp = 'n'
    if (eclipses): eclp = 'y'
    call += 'eclipses='+eclp+ ' '
    boost = 'n'
    if (dopboost): boost = 'y'
    call += 'dopboost='+boost+ ' '
    distort = 'n'
    if (tides): distort = 'y'
    call += 'tides='+distort+ ' '
    call += 'job='+str(job)+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPBINARY started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# check and format the list of fit parameters

    if status == 0 and job == 'fit':
        allParams = [m1,m2,r1,r2,period,bjd0,eccn,omega,inclination]
        allNames = ['m1','m2','r1','r2','period','bjd0','eccn','omega','inclination']
        fitparams = re.sub('\|',',',fitparams.strip())
        fitparams = re.sub('\.',',',fitparams.strip())
        fitparams = re.sub(';',',',fitparams.strip())
        fitparams = re.sub(':',',',fitparams.strip())
        fitparams = re.sub('\s+',',',fitparams.strip())
        fitparams, status = kepio.parselist(fitparams,logfile,verbose)
        for fitparam in fitparams:
            if fitparam.strip() not in allNames:
                message = 'ERROR -- KEPBINARY: unknown field in list of fit parameters'
                status = kepmsg.err(logfile,message,verbose)

# clobber output file

    if status == 0:
        if clobber: status = kepio.clobber(outfile,logfile,verbose)
        if kepio.fileexists(outfile): 
            message = 'ERROR -- KEPBINARY: ' + outfile + ' exists. Use --clobber'
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

# check the data column exists

    if status == 0:
        try:
            instr[1].data.field(datacol)
        except:
            message = 'ERROR -- KEPBINARY: ' + datacol + ' column does not exist in ' + infile + '[1]'
            status = kepmsg.err(logfile,message,verbose)

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
                    if numpy.isfinite(table.field('barytime')[i]) and \
                            numpy.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
            except:
                for i in range(len(table.field(0))):
                    if numpy.isfinite(table.field('time')[i]) and \
                            numpy.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
            comment = 'NaN cadences removed from data'
            status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)
 
# read table columns

    if status == 0:
	try:
            time = instr[1].data.field('barytime')
	except:
            time, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
	indata, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)
    if status == 0:
        time = time + bjdref
        indata = indata / cadenom

# limb-darkening cofficients

    if status == 0:
        limbdark = numpy.array([c1,c2,c3,c4],dtype='float32')

# time details for model

    if status == 0:
        npt = len(time)
        exptime = numpy.zeros((npt),dtype='float64')
        dtype = numpy.zeros((npt),dtype='int')
        for i in range(npt):
            try:
                exptime[i] = time[i+1] - time[i]
            except:
                exptime[i] = time[i] - time[i-1]

# calculate binary model

    if status == 0:
        tmodel = kepsim.transitModel(1.0,m1,m2,r1,r2,period,inclination,bjd0,eccn,omega,depth,
                                     albedo,c1,c2,c3,c4,gamma,contamination,npt,time,exptime,
                                     dtype,eclipses,dopboost,tides)

# re-normalize binary model to data

    if status == 0 and (job == 'overlay' or job == 'fit'):
        dmedian = numpy.median(indata)
        tmodel = tmodel / numpy.median(tmodel) * dmedian

# define arrays of floating and frozen parameters

    if status == 0 and job =='fit':
        params = []; paramNames = []; arguments = []; argNames = []
        for i in range(len(allNames)):
            if allNames[i] in fitparams:
                params.append(allParams[i])
                paramNames.append(allNames[i])
            else:
                arguments.append(allParams[i])
                argNames.append(allNames[i])
        params.append(dmedian)
        params = numpy.array(params,dtype='float32')

# subtract model from data

    if status == 0 and job == 'fit':
        deltam = numpy.abs(indata - tmodel)

# fit statistics

    if status == 0 and job == 'fit':
        aveDelta = numpy.sum(deltam) / npt
        chi2 = math.sqrt(numpy.sum((indata - tmodel) * (indata - tmodel) / (npt - len(params))))

# fit model to data using downhill simplex

    if status == 0 and job == 'fit':
        print ''
        print '%4s %11s %11s' % ('iter', 'delta', 'chi^2')
        print '----------------------------'
        print '%4d %.5E %.5E' % (0,aveDelta,chi2)
        bestFit = scipy.optimize.fmin(fitModel,params,args=(paramNames,dmedian,m1,m2,r1,r2,period,bjd0,eccn,
                                                            omega,inclination,depth,albedo,c1,c2,c3,c4,
                                                            gamma,contamination,npt,time,exptime,indata,
                                                            dtype,eclipses,dopboost,tides),maxiter=1e4)

# calculate best fit binary model

    if status == 0 and job == 'fit':
        print ''
        for i in range(len(paramNames)):
            if 'm1' in paramNames[i].lower():
                m1 = bestFit[i]
                print '  M1 = %.3f Msun' % bestFit[i]
            elif 'm2' in paramNames[i].lower():
                m2 = bestFit[i]
                print '  M2 = %.3f Msun' % bestFit[i]
            elif 'r1' in paramNames[i].lower():
                r1 = bestFit[i]
                print '  R1 = %.4f Rsun' % bestFit[i]
            elif 'r2' in paramNames[i].lower():
                r2 = bestFit[i]
                print '  R2 = %.4f Rsun' % bestFit[i]
            elif 'period' in paramNames[i].lower():
                period = bestFit[i]
            elif 'bjd0' in paramNames[i].lower():
                bjd0 = bestFit[i]
                print 'BJD0 = %.8f' % bestFit[i]
            elif 'eccn' in paramNames[i].lower():
                eccn = bestFit[i]
                print '   e = %.3f' % bestFit[i]
            elif 'omega' in paramNames[i].lower():
                omega = bestFit[i]
                print '   w = %.3f deg' % bestFit[i]
            elif 'inclination' in paramNames[i].lower():
                inclination = bestFit[i]
                print '   i = %.3f deg' % bestFit[i]
        flux = bestFit[-1]
        print ''
        tmodel = kepsim.transitModel(flux,m1,m2,r1,r2,period,inclination,bjd0,eccn,omega,depth,
                                     albedo,c1,c2,c3,c4,gamma,contamination,npt,time,exptime,
                                     dtype,eclipses,dopboost,tides)

# subtract model from data

    if status == 0:
        deltaMod = indata - tmodel

# standard deviation of model

    if status == 0:
        stdDev = math.sqrt(numpy.sum((indata - tmodel) * (indata - tmodel)) / npt)

# clean up x-axis unit

    if status == 0:
	time0 = float(int(tstart / 100) * 100.0)
	ptime = time - time0
	xlab = 'BJD $-$ %d' % time0

# clean up y-axis units

    if status == 0:
	nrm = len(str(int(indata.max())))-1
	pout = indata / 10**nrm
	pmod = tmodel / 10**nrm
        pres = deltaMod / stdDev
        if job == 'fit' or job == 'overlay':
            try:
                ylab1 = 'Flux (10$^%d$ e$^-$ s$^{-1}$)' % nrm
                ylab2 = 'Residual ($\sigma$)'
            except:
                ylab1 = 'Flux (10**%d e-/s)' % nrm
                ylab2 = 'Residual (sigma)'
        else:
            ylab1 = 'Normalized Flux'

# dynamic range of model plot

    if status == 0 and job == 'model':
        xmin = ptime.min()
        xmax = ptime.max()
        ymin = tmodel.min()
        ymax = tmodel.max()

# dynamic range of model/data overlay or fit

    if status == 0 and (job == 'overlay' or job == 'fit'):
        xmin = ptime.min()
        xmax = ptime.max()
        ymin = pout.min()
        ymax = pout.max()
        tmin = pmod.min()
        tmax = pmod.max()
        ymin = numpy.array([ymin,tmin]).min()
        ymax = numpy.array([ymax,tmax]).max()
        rmin = pres.min()
        rmax = pres.max()
 
# pad the dynamic range

    if status == 0:
        xr = (xmax - xmin) / 80
        yr = (ymax - ymin) / 40
        if job == 'overlay' or job == 'fit':
            rr = (rmax - rmin) / 40

# set up plot style

    if status == 0:
        labelsize = 24; ticksize = 16; xsize = 17; ysize = 7
        lcolor = '#0000ff'; lwidth = 1.0; fcolor = '#ffff00'; falpha = 0.2
        params = {'backend': 'png',
                  'axes.linewidth': 2.5,
                  'axes.labelsize': 24,
                  'axes.font': 'sans-serif',
                  'axes.fontweight' : 'bold',
                  'text.fontsize': 12,
                  'legend.fontsize': 12,
                  'xtick.labelsize': 16,
                  'ytick.labelsize': 16}
        pylab.rcParams.update(params)
        pylab.figure(figsize=[14,10])
        pylab.clf()

# main plot window

        ax = pylab.axes([0.05,0.3,0.94,0.68])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        labels = ax.get_yticklabels()
        setp(labels, 'rotation', 90, fontsize=12)

# plot model time series 

    if status == 0 and job == 'model':
        pylab.plot(ptime,tmodel,color='#0000ff',linestyle='-',linewidth=1.0)
        ptime = numpy.insert(ptime,[0.0],ptime[0])
        ptime = numpy.append(ptime,ptime[-1])
        tmodel = numpy.insert(tmodel,[0.0],0.0)
        tmodel = numpy.append(tmodel,0.0)
        pylab.fill(ptime,tmodel,fc='#ffff00',linewidth=0.0,alpha=0.2)

# plot data time series and best fit

    if status == 0 and (job == 'overlay' or job == 'fit'):
        pylab.plot(ptime,pout,color='#0000ff',linestyle='-',linewidth=1.0)
        ptime = numpy.insert(ptime,[0.0],ptime[0])
        ptime = numpy.append(ptime,ptime[-1])
        pout = numpy.insert(pout,[0],0.0)
        pout = numpy.append(pout,0.0)
        pylab.fill(ptime,pout,fc='#ffff00',linewidth=0.0,alpha=0.2)
        pylab.plot(ptime[1:-1],pmod,color='r',linestyle='-',linewidth=2.0)

# ranges and labels

    if status == 0:
        pylab.xlim(xmin-xr,xmax+xr)
        pylab.ylim(ymin-yr,ymax+yr)
        pylab.xlabel(xlab, {'color' : 'k'})
        pylab.ylabel(ylab1, {'color' : 'k'})

# residual plot window

    if status == 0 and (job == 'overlay' or job == 'fit'):
        ax = pylab.axes([0.05,0.07,0.94,0.23])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        labels = ax.get_yticklabels()
        setp(labels, 'rotation', 90, fontsize=12)

# plot residual time series 

    if status == 0 and (job == 'overlay' or job == 'fit'):
        pylab.plot([ptime[0],ptime[-1]],[0.0,0.0],color='r',linestyle='--',linewidth=1.0)
        pylab.plot([ptime[0],ptime[-1]],[-1.0,-1.0],color='r',linestyle='--',linewidth=1.0)
        pylab.plot([ptime[0],ptime[-1]],[1.0,1.0],color='r',linestyle='--',linewidth=1.0)
        pylab.plot(ptime[1:-1],pres,color='#0000ff',linestyle='-',linewidth=1.0)
        pres = numpy.insert(pres,[0],rmin)
        pres = numpy.append(pres,rmin)
        pylab.fill(ptime,pres,fc='#ffff00',linewidth=0.0,alpha=0.2)

# ranges and labels of residual time series

    if status == 0 and (job == 'overlay' or job == 'fit'):
        pylab.xlim(xmin-xr,xmax+xr)
        pylab.ylim(rmin-rr,rmax+rr)
        pylab.xlabel(xlab, {'color' : 'k'})
        pylab.ylabel(ylab2, {'color' : 'k'})

# display the plot

    if status == 0:
        pylab.draw()

# -----------------------------------------------------------
# fit function

def fitModel(params,paramNames,flux,m1,m2,r1,r2,period,bjd0,eccn,omega,inclination,depth,albedo,
             c1,c2,c3,c4,gamma,contamination,npt,time,exptime,indata,dtype,eclipses,dopboost,tides):

# iteration number

    global niter
    niter += 1

# unpack first-guess model variables

    for i in range(len(paramNames)):
        if 'm1' in paramNames[i].lower():
            m1 = params[i]
        elif 'm2' in paramNames[i].lower():
            m2 = params[i]
        elif 'r1' in paramNames[i].lower():
            r1 = params[i]
        elif 'r2' in paramNames[i].lower():
            r2 = params[i]
        elif 'period' in paramNames[i].lower():
            period = params[i]
        elif 'bjd0' in paramNames[i].lower():
            bjd0 = params[i]
        elif 'eccn' in paramNames[i].lower():
            eccn = params[i]
        elif 'omega' in paramNames[i].lower():
            omega = params[i]
        elif 'inclination' in paramNames[i].lower():
            inclination = params[i]
    flux = params[-1]

# constuct time series model

    tmodel = kepsim.transitModel(flux,m1,m2,r1,r2,period,inclination,bjd0,eccn,omega,depth,
                                        albedo,c1,c2,c3,c4,gamma,contamination,npt,time,exptime,
                                        dtype,eclipses,dopboost,tides)

# subtract model from data

    deltam = numpy.abs(indata - tmodel)

# fit statistics

    aveDelta = numpy.sum(deltam) / npt
    chi2 = math.sqrt(numpy.sum((indata - tmodel) * (indata - tmodel) / (npt - len(params))))

# print fit

    print '%4d %.5E %.5E' % (niter,aveDelta,chi2)

    return chi2

# -----------------------------------------------------------
# main

parfile = iraf.osfn("kepler$kepbinary.par")
t = iraf.IrafTaskFactory(taskname="kepbinary", value=parfile, function=kepbinary)
