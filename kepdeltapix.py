import pylab, numpy
from pylab import *
from matplotlib import *
from numpy import *
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey, kepplot, kepfit, keparray
import sys, time, re, math, glob, random

# -----------------------------------------------------------
# core code

def kepdeltapix(infile,nexp,columns,rows,fluxes,prfdir,interpolation,tolerance,fittype,imscale,
           colmap,verbose,logfile,status,cmdLine=False): 

# input arguments

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDELTAPIX -- '
    call += 'infile='+infile+' '
    call += 'nexp='+str(nexp)+' '
    call += 'columns='+columns+' '
    call += 'rows='+rows+' '
    call += 'fluxes='+fluxes+' '
    call += 'prfdir='+prfdir+' '
    call += 'interpolation='+interpolation+' '
    call += 'tolerance='+str(tolerance)+' '
    call += 'fittype='+str(fittype)+' '
    call += 'imscale='+imscale+' '
    call += 'colmap='+colmap+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# test log file

    logfile = kepmsg.test(logfile)

# start time

    kepmsg.clock('KEPDELTAPIX started at',logfile,verbose)

# reference color map

    if colmap == 'browse':
        status = cmap_plot(cmdLine)

# open TPF FITS file

    if status == 0:
        try:
            kepid, channel, skygroup, module, output, quarter, season, \
                ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
                kepio.readTPF(infile,'TIME',logfile,verbose)
        except:
            message = 'ERROR -- KEPDELTAPIX: is %s a Target Pixel File? ' % infile
            status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, tcorr, status = \
            kepio.readTPF(infile,'TIMECORR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, cadno, status = \
            kepio.readTPF(infile,'CADENCENO',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, fluxpixels, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, errpixels, status = \
            kepio.readTPF(infile,'FLUX_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)

# print target data

    if status == 0 and verbose:
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

# determine suitable PRF calibration file

    if status == 0:
        if int(module) < 10:
            prefix = 'kplr0'
        else:
            prefix = 'kplr'
        prfglob = prfdir + '/' + prefix + str(module) + '.' + str(output) + '*' + '_prf.fits'
        try:
            prffile = glob.glob(prfglob)[0]
        except:
            message = 'ERROR -- KEPDELTAPIX: No PRF file found in ' + prfdir
            status = kepmsg.err(logfile,message,verbose)

# read PRF images

    if status == 0:
        prfn = [0,0,0,0,0]
        crpix1p = numpy.zeros((5),dtype='float32')
        crpix2p = numpy.zeros((5),dtype='float32')
        crval1p = numpy.zeros((5),dtype='float32')
        crval2p = numpy.zeros((5),dtype='float32')
        cdelt1p = numpy.zeros((5),dtype='float32')
        cdelt2p = numpy.zeros((5),dtype='float32')
        for i in range(5):
            prfn[i], crpix1p[i], crpix2p[i], crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i], status \
                = kepio.readPRFimage(prffile,i+1,logfile,verbose)    

# choose rows in the TPF table at random

    if status == 0:
        i = 0
        rownum = []
        while i < nexp:
            work = int(random.random() * len(barytime))
            if numpy.isfinite(barytime[work]) and numpy.isfinite(fluxpixels[work,ydim*xdim/2]):
                rownum.append(work)
                i += 1

# construct input pixel image

    if status == 0:
        fscat = numpy.empty((len(fluxes),nexp),dtype='float32')
        xscat = numpy.empty((len(columns),nexp),dtype='float32')
        yscat = numpy.empty((len(rows),nexp),dtype='float32')
        for irow in range(nexp):
            flux = fluxpixels[rownum[irow],:]

# image scale and intensity limits of pixel data

            if status == 0:
                flux_pl, zminfl, zmaxfl = kepplot.intScale1D(flux,imscale)
                n = 0
                imgflux_pl = empty((ydim,xdim))
                for i in range(ydim):
                    for j in range(xdim):
                        imgflux_pl[i,j] = flux_pl[n]
                        n += 1
        
# fit PRF model to pixel data

            if status == 0:
                start = time.time()
                f,y,x,prfMod,prfFit,prfRes = kepfit.fitMultiPRF(flux,ydim,xdim,column,row,prfn,crval1p,
                     crval2p,cdelt1p,cdelt2p,interpolation,tolerance,fluxes,columns,rows,fittype,
                     verbose,logfile)
                if verbose:
                    print '\nConvergence time = %.1fs' % (time.time() - start)

# best fit parameters

            if status == 0:
                for i in range(len(f)):
                    fscat[i,irow] = f[i]
                    xscat[i,irow] = x[i]
                    yscat[i,irow] = y[i]

# replace starting guess with previous fit parameters

            if status == 0:
                fluxes = copy(f)
                columns = copy(x)
                rows = copy(y)

# mean and rms results

    if status == 0:
        fmean = []; fsig = []
        xmean = []; xsig = []
        ymean = []; ysig = []
        for i in range(len(f)):
            fmean.append(numpy.mean(fscat[i,:]))
            xmean.append(numpy.mean(xscat[i,:]))
            ymean.append(numpy.mean(yscat[i,:]))
            fsig.append(numpy.std(fscat[i,:]))
            xsig.append(numpy.std(xscat[i,:]))
            ysig.append(numpy.std(yscat[i,:]))
            txt = 'Flux = %10.2f e-/s ' % fmean[-1]
            txt += 'X = %7.4f +/- %6.4f pix ' % (xmean[-1], xsig[i])
            txt += 'Y = %7.4f +/- %6.4f pix' % (ymean[-1], ysig[i])
            kepmsg.log(logfile,txt,True)

# output results for kepprfphot

    if status == 0:
        txt1 = 'columns=0.0'
        txt2 = '   rows=0.0'
        for i in range(1,len(f)):
            txt1 += ',%.4f' % (xmean[i] - xmean[0])
            txt2 += ',%.4f' % (ymean[i] - ymean[0])
        kepmsg.log(logfile,'\nkepprfphot input fields:',True)
        kepmsg.log(logfile,txt1,True)
        kepmsg.log(logfile,txt2,True)
        
# image scale and intensity limits for PRF model image

    if status == 0:
        imgprf_pl, zminpr, zmaxpr = kepplot.intScale2D(prfMod,imscale)
        
# image scale and intensity limits for PRF fit image

    if status == 0:
        imgfit_pl, zminfi, zmaxfi = kepplot.intScale2D(prfFit,imscale)
        
# image scale and intensity limits for data - fit residual

    if status == 0:
        imgres_pl, zminre, zmaxre = kepplot.intScale2D(prfRes,imscale)
        
# plot style

    if status == 0:
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': 24,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': 10,
                      'ytick.labelsize': 10}
            pylab.rcParams.update(params)
        except:
            pass
        pylab.figure(figsize=[10,10])
        pylab.clf()
        plotimage(imgflux_pl,zminfl,zmaxfl,1,row,column,xdim,ydim,0.06,0.52,'flux',colmap)
        plotimage(imgfit_pl,zminfl,zmaxfl,3,row,column,xdim,ydim,0.06,0.06,'fit',colmap)
        plotimage(imgres_pl,zminfl,zmaxfl,4,row,column,xdim,ydim,0.52,0.06,'residual',colmap)
        plotimage(imgprf_pl,zminpr,zmaxpr*0.9,2,row,column,xdim,ydim,0.52,0.52,'model',colmap)
        for i in range(len(f)):        
            pylab.plot(xscat[i,:],yscat[i,:],'o',color='k')            
            
# Plot creep of target position over time, relative to the central source

#	barytime0 = float(int(barytime[0] / 100) * 100.0)
#	barytime -= barytime0
#        xlab = 'BJD $-$ %d' % barytime0
#	xmin = numpy.nanmin(barytime)
#	xmax = numpy.nanmax(barytime)
#	y1min = numpy.nanmin(data)
#	y1max = numpy.nanmax(data)
#	xr = xmax - xmin
#	yr = ymax - ymin
#        barytime = insert(barytime,[0],[barytime[0]]) 
#        barytime = append(barytime,[barytime[-1]])
#        data = insert(data,[0],[0.0]) 
#        data = append(data,0.0)
#
#        pylab.figure(2,figsize=[10,10])
#        pylab.clf()
#	ax = pylab.subplot(211)
#	pylab.subplots_adjust(0.1,0.5,0.88,0.42)
#        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
#        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
#        labels = ax.get_yticklabels()
#        setp(labels, 'rotation', 90, fontsize=ticksize)
#        for i in range(1,len(f)):
#            pylab.plot(rownum,xscat[i,:]-xscat[0,:],'o')
#        pylab.ylabel('$\Delta$Columns', {'color' : 'k'})
#	ax = pylab.subplot(211)
#	pylab.subplots_adjust(0.1,0.1,0.88,0.42)
#        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
#        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
#        labels = ax.get_yticklabels()
#        setp(labels, 'rotation', 90, fontsize=ticksize)
#        for i in range(1,len(f)):
#            pylab.plot(rownum,yscat[i,:]-yscat[0,:],'o')
#	pylab.xlim(xmin-xr*0.01,xmax+xr*0.01)
#	if ymin-yr*0.01 <= 0.0 or fullrange:
#            pylab.ylim(1.0e-10,ymax+yr*0.01)
#	else:
#            pylab.ylim(ymin-yr*0.01,ymax+yr*0.01)
#        pylab.ylabel('$\Delta$Rows', {'color' : 'k'})
#        pylab.xlabel(xlab, {'color' : 'k'})

# render plot

    if status == 0:
        if cmdLine: 
            pylab.show(block=True)
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
# stop time

    kepmsg.clock('\nKEPDELTAPIX ended at',logfile,verbose)

    return

# -----------------------------------------------------------
# plot channel image

def plotimage(imgflux_pl,zminfl,zmaxfl,plmode,row,column,xdim,ydim,winx,winy,tlabel,colmap):

# pixel limits of the subimage

    ymin = row
    ymax = ymin + ydim
    xmin = column
    xmax = xmin + xdim

# plot limits for flux image

    ymin = float(ymin) - 0.5
    ymax = float(ymax) - 0.5
    xmin = float(xmin) - 0.5
    xmax = float(xmax) - 0.5

# plot the image window

    ax = pylab.axes([winx,winy,0.46,0.46])
    pylab.imshow(imgflux_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zminfl,vmax=zmaxfl,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    pylab.gca().set_autoscale_on(False)
    labels = ax.get_yticklabels()
    setp(labels, 'rotation', 90)
    pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    if plmode == 1:
        pylab.setp(pylab.gca(),xticklabels=[])
    if plmode == 2:
        pylab.setp(pylab.gca(),xticklabels=[],yticklabels=[])
    if plmode == 4:
        pylab.setp(pylab.gca(),yticklabels=[])
    if plmode == 3 or plmode == 4:
        pylab.xlabel('Pixel Column Number', {'color' : 'k'})
    if plmode == 1 or plmode == 3:
        pylab.ylabel('Pixel Row Number', {'color' : 'k'})
    pylab.text(0.05, 0.93,tlabel,horizontalalignment='left',verticalalignment='center',
               fontsize=28,fontweight=500,transform=ax.transAxes)

    return

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot(cmdLine):

    pylab.figure(figsize=[5,10])
    a=outer(ones(10),arange(0,1,0.01))
    subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        print m
        subplot(l,1,i+1)
        pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        pylab.text(100.85,0.5,m,fontsize=10)

# render plot

    if cmdLine: 
        pylab.show(block=True)
    else: 
        pylab.ion()
        pylab.plot([])
        pylab.ioff()
	
    status = 1
    return status

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Fitting PRF model to Target Pixel image')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('--nexp', '-n', default=10, help='Number of random pixel images to sample', dest='nexp', type=int)
    parser.add_argument('columns', help='Column number of each source to be fit', type=str)
    parser.add_argument('rows', help='Row number of each source to be fit', type=str)
    parser.add_argument('fluxes', help='Relative flux of each source to be fit', type=str)
    parser.add_argument('prfdir', help='Folder containing Point Response Function FITS files', type=str)
    parser.add_argument('--interpolation', '-n', help='Pixel Response Function interpolation scheme', default='linear', dest='interpolation', type=str,choices=['linear','spline'])
    parser.add_argument('--tolerance', '-t', default=1.0e-4, help='Fit tolerance', dest='tolerance', type=float)
    parser.add_argument('--fittype', '-f', help='Fitting scheme', default='2D', dest='fittype', type=str,choices=['1D','2D'])
    parser.add_argument('--imscale', '-i', help='Type of image intensity scale', default='linear', dest='imscale', type=str,choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', '-c', help='Image colormap', default='YlOrBr', dest='cmap', type=str,choices=['Accent','Blues','BrBG','BuGn','BuPu','Dark2','GnBu','Greens','Greys','OrRd','Oranges','PRGn','Paired','Pastel1','Pastel2','PiYG','PuBu','PuBuGn','PuOr','PuRd','Purples','RdBu','RdGy','RdPu','RdYlBu','RdYlGn','Reds','Set1','Set2','Set3','Spectral','YlGn','YlGnBu','YlOrBr','YlOrRd','afmhot','autumn','binary','bone','brg','bwr','cool','copper','flag','gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_yarg','gnuplot','gnuplot2','gray','hot','hsv','jet','ocean','pink','prism','rainbow','seismic','spectral','spring','summer','terrain','winter','browse'])
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepdeltapix.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepdeltapix(args.infile,args.rownum,args.columns,args.rows,args.fluxes,args.prfdir,args.interpolation,
           args.tolerance,args.fittype,args.imscale,args.cmap,args.verbose,args.logfile,args.status,cmdLine)
    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepdeltapix.par")
    t = iraf.IrafTaskFactory(taskname="kepdeltapix", value=parfile, function=kepdeltapix)
