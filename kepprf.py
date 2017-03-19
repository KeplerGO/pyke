import pylab, numpy, scipy
from astropy.io import fits as pyfits
from pylab import *
from matplotlib import *
from matplotlib import pyplot
from numpy import *
import kepio, kepmsg, kepkey, kepplot, kepfit, keparray, kepfunc, kepstat
import sys, time, re, math, glob
from scipy import interpolate, optimize, ndimage, stats
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline, interp2d
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate

# -----------------------------------------------------------
# core code

def kepprf(infile,plotfile,rownum,columns,rows,fluxes,border,background,focus,prfdir,xtol,ftol,
           imscale,colmap,labcol,apercol,plt,verbose,logfile,status,cmdLine=False):

# input arguments

    status = 0
    seterr(all="ignore")

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPPRF -- '
    call += 'infile='+infile+' '
    call += 'plotfile='+plotfile+' '
    call += 'rownum='+str(rownum)+' '
    call += 'columns='+columns+' '
    call += 'rows='+rows+' '
    call += 'fluxes='+fluxes+' '
    call += 'border='+str(border)+' '
    bground = 'n'
    if (background): bground = 'y'
    call += 'background='+bground+' '
    focs = 'n'
    if (focus): focs = 'y'
    call += 'focus='+focs+' '
    call += 'prfdir='+prfdir+' '
    call += 'xtol='+str(xtol)+' '
    call += 'ftol='+str(xtol)+' '
    call += 'imscale='+imscale+' '
    call += 'colmap='+colmap+' '
    call += 'labcol='+labcol+' '
    call += 'apercol='+apercol+' '
    plotit = 'n'
    if (plt): plotit = 'y'
    call += 'plot='+plotit+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# test log file

    logfile = kepmsg.test(logfile)

# start time

    kepmsg.clock('KEPPRF started at',logfile,verbose)

# reference color map

    if colmap == 'browse':
        status = cmap_plot(cmdLine)

# construct inital guess vector for fit

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            x = columns.strip().split(',')
            y = rows.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i])
        except:
            f = fluxes
            x = columns
            y = rows
        nsrc = len(f)
        for i in xrange(nsrc):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            if len(x) != nsrc or len(y) != nsrc:
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(nsrc):
                try:
                    guess.append(float(x[i]))
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(nsrc):
                try:
                    guess.append(float(y[i]))
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0 and background:
            if border == 0:
                guess.append(0.0)
            else:
                for i in range((border+1)*2):
                    guess.append(0.0)
        if status == 0 and focus:
            guess.append(1.0); guess.append(1.0); guess.append(0.0)

# open TPF FITS file

    if status == 0:
        try:
            kepid, channel, skygroup, module, output, quarter, season, \
                ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
                kepio.readTPF(infile,'TIME',logfile,verbose)
        except:
            message = 'ERROR -- KEPPRF: is %s a Target Pixel File? ' % infile
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

# read mask defintion data from TPF file

    if status == 0:
        maskimg, pixcoord1, pixcoord2, status = kepio.readMaskDefinition(infile,logfile,verbose)
        npix = numpy.size(numpy.nonzero(maskimg)[0])

# print target data

    if status == 0 and verbose:
        print ''
        print '      KepID: %s' % kepid
        print '        BJD: %.2f' % (barytime[rownum-1] + 2454833.0)
        print ' RA (J2000): %s' % ra
        print 'Dec (J2000):  %s' % dec
        print '     KepMag:  %s' % kepmag
        print '   SkyGroup:   %2s' % skygroup
        print '     Season:   %2s' % str(season)
        print '    Channel:   %2s' % channel
        print '     Module:   %2s' % module
        print '     Output:    %1s' % output
        print ''

# is this a good row with finite timestamp and pixels?

    if status == 0:
        if not numpy.isfinite(barytime[rownum-1]) or numpy.nansum(fluxpixels[rownum-1,:]) == numpy.nan:
            message = 'ERROR -- KEPFIELD: Row ' + str(rownum) + ' is a bad quality timestamp'
            status = kepmsg.err(logfile,message,verbose)

# construct input pixel image

    if status == 0:
        flux = fluxpixels[rownum-1,:]
        ferr = errpixels[rownum-1,:]
        DATx = arange(column,column+xdim)
        DATy = arange(row,row+ydim)
#        if numpy.nanmin > 420000.0: flux -= 420000.0

# image scale and intensity limits of pixel data

    if status == 0:
        n = 0
        DATimg = empty((ydim,xdim))
        ERRimg = empty((ydim,xdim))
        for i in range(ydim):
            for j in range(xdim):
                DATimg[i,j] = flux[n]
                ERRimg[i,j] = ferr[n]
                n += 1

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
            message = 'ERROR -- KEPPRF: No PRF file found in ' + prfdir
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
        prfn = array(prfn)
        PRFx = arange(0.5,shape(prfn[0])[1]+0.5)
        PRFy = arange(0.5,shape(prfn[0])[0]+0.5)
        PRFx = (PRFx - size(PRFx) / 2) * cdelt1p[0]
        PRFy = (PRFy - size(PRFy) / 2) * cdelt2p[0]

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e-6
            prf = prf + prfn[i] / prfWeight[i]
        prf = prf / nansum(prf) / cdelt1p[0] / cdelt2p[0]

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfDimY = int(ydim / cdelt1p[0])
        prfDimX = int(xdim / cdelt2p[0])
        PRFy0 = (shape(prf)[0] - prfDimY) / 2
        PRFx0 = (shape(prf)[1] - prfDimX) / 2

# interpolation function over the PRF

    if status == 0:
        splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)

# construct mesh for background model

    if status == 0 and background:
        bx = numpy.arange(1.,float(xdim+1))
        by = numpy.arange(1.,float(ydim+1))
        xx, yy = numpy.meshgrid(numpy.linspace(bx.min(), bx.max(), xdim),
                                numpy.linspace(by.min(), by.max(), ydim))

# fit PRF model to pixel data

    if status == 0:
        start = time.time()
        if focus and background:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,border,xx,yy,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithFocusAndBackground,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        elif focus and not background:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithFocus,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        elif background and not focus:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,border,xx,yy,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithBackground,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        else:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRF,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        print 'Convergence time = %.2fs\n' % (time.time() - start)

# pad the PRF data if the PRF array is smaller than the data array

    if status == 0:
        flux = []; OBJx = []; OBJy = []
        PRFmod = numpy.zeros((prfDimY,prfDimX))
        if PRFy0 < 0 or PRFx0 < 0.0:
            PRFmod = numpy.zeros((prfDimY,prfDimX))
            superPRF = zeros((prfDimY+1,prfDimX+1))
            superPRF[abs(PRFy0):abs(PRFy0)+shape(prf)[0],abs(PRFx0):abs(PRFx0)+shape(prf)[1]] = prf
            prf = superPRF * 1.0
            PRFy0 = 0
            PRFx0 = 0

# rotate the PRF model around its center

        if focus:
            angle = ans[-1]
            prf = rotate(prf,-angle,reshape=False,mode='nearest')

# iterate through the sources in the best fit PSF model

        for i in range(nsrc):
            flux.append(ans[i])
            OBJx.append(ans[nsrc+i])
            OBJy.append(ans[nsrc*2+i])

# calculate best-fit model

            y = (OBJy[i]-mean(DATy)) / cdelt1p[0]
            x = (OBJx[i]-mean(DATx)) / cdelt2p[0]
            prfTmp = shift(prf,[y,x],order=3,mode='constant')
            prfTmp = prfTmp[PRFy0:PRFy0+prfDimY,PRFx0:PRFx0+prfDimX]
            PRFmod = PRFmod + prfTmp * flux[i]
            wx = 1.0
            wy = 1.0
            angle = 0
            b = 0.0

# write out best fit parameters

            if verbose:
                txt = 'Flux = %10.2f e-/s ' % flux[i]
                txt += 'X = %9.4f pix ' % OBJx[i]
                txt += 'Y = %9.4f pix ' % OBJy[i]
                kepmsg.log(logfile,txt,True)

        if verbose and background:
            bterms = border + 1
            if bterms == 1:
                b = ans[nsrc*3]
            else:
                bcoeff = array([ans[nsrc*3:nsrc*3+bterms],ans[nsrc*3+bterms:nsrc*3+bterms*2]])
                bkg = kepfunc.polyval2d(xx,yy,bcoeff)
                b = nanmean(bkg.reshape(bkg.size))
            txt = '\n   Mean background = %.2f e-/s' % b
            kepmsg.log(logfile,txt,True)
        if focus:
            wx = ans[-3]
            wy = ans[-2]
            angle = ans[-1]
        if verbose and focus:
            if not background: kepmsg.log(logfile,'',True)
            kepmsg.log(logfile,' X/Y focus factors = %.3f/%.3f' % (wx,wy),True)
            kepmsg.log(logfile,'PRF rotation angle = %.2f deg' % angle,True)

# measure flux fraction and contamination

    if status == 0:
        PRFall = kepfunc.PRF2DET(flux,OBJx,OBJy,DATx,DATy,wx,wy,angle,splineInterpolation)
        PRFone = kepfunc.PRF2DET([flux[0]],[OBJx[0]],[OBJy[0]],DATx,DATy,wx,wy,angle,splineInterpolation)
        FluxInMaskAll = numpy.nansum(PRFall)
        FluxInMaskOne = numpy.nansum(PRFone)
        FluxInAperAll = 0.0
        FluxInAperOne = 0.0
        for i in range(1,ydim):
            for j in range(1,xdim):
                if kepstat.bitInBitmap(maskimg[i,j],2):
                    FluxInAperAll += PRFall[i,j]
                    FluxInAperOne += PRFone[i,j]
        FluxFraction = FluxInAperOne / flux[0]
        try:
            Contamination = (FluxInAperAll - FluxInAperOne) / FluxInAperAll
        except:
            Contamination = 0.0

        kepmsg.log(logfile,'\n                Total flux in mask = %.2f e-/s' % FluxInMaskAll,True)
        kepmsg.log(logfile,'               Target flux in mask = %.2f e-/s' % FluxInMaskOne,True)
        kepmsg.log(logfile,'            Total flux in aperture = %.2f e-/s' % FluxInAperAll,True)
        kepmsg.log(logfile,'           Target flux in aperture = %.2f e-/s' % FluxInAperOne,True)
        kepmsg.log(logfile,'  Target flux fraction in aperture = %.2f%%' % (FluxFraction * 100.0),True)
        kepmsg.log(logfile,'Contamination fraction in aperture = %.2f%%' % (Contamination * 100.0),True)


# constuct model PRF in detector coordinates

    if status == 0:
        PRFfit = PRFall + 0.0
        if background and bterms == 1:
            PRFfit = PRFall + b
        if background and bterms > 1:
            PRFfit = PRFall + bkg

# calculate residual of DATA - FIT

    if status == 0:
        PRFres = DATimg - PRFfit
        FLUXres = numpy.nansum(PRFres) / npix

# calculate the sum squared difference between data and model

    if status == 0:
        Pearson = abs(numpy.nansum(numpy.square(DATimg - PRFfit) / PRFfit))
        Chi2 = numpy.nansum(numpy.square(DATimg - PRFfit) / numpy.square(ERRimg))
        DegOfFreedom = npix - len(guess) - 1
        try:
            kepmsg.log(logfile,'\n       Residual flux = %.2f e-/s' % FLUXres,True)
            kepmsg.log(logfile,'Pearson\'s chi^2 test = %d for %d dof' % (Pearson,DegOfFreedom),True)
        except:
            pass
        kepmsg.log(logfile,'          Chi^2 test = %d for %d dof' % (Chi2,DegOfFreedom),True)

# image scale and intensity limits for plotting images

    if status == 0:
        imgdat_pl, zminfl, zmaxfl = kepplot.intScale2D(DATimg,imscale)
        imgprf_pl, zminpr, zmaxpr = kepplot.intScale2D(PRFmod,imscale)
        imgfit_pl, zminfi, zmaxfi = kepplot.intScale2D(PRFfit,imscale)
        imgres_pl, zminre, zmaxre = kepplot.intScale2D(PRFres,'linear')
        if imscale == 'linear':
            zmaxpr *= 0.9
        elif imscale == 'logarithmic':
            zmaxpr = numpy.max(zmaxpr)
            zminpr = zmaxpr / 2

# plot style

    if status == 0:
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': 28,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': 20,
                      'ytick.labelsize': 20,
                      'xtick.major.pad': 6,
                      'ytick.major.pad': 6}
            pylab.rcParams.update(params)
        except:
            pass
        pylab.figure(figsize=[12,10])
        pylab.clf()
        plotimage(imgdat_pl,zminfl,zmaxfl,1,row,column,xdim,ydim,0.07,0.53,'observation',colmap,labcol)
        plotimage(imgprf_pl,zminpr,zmaxpr,2,row,column,xdim,ydim,0.44,0.53,'model',colmap,labcol)
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,1,apercol,'--',0.5)
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,2,apercol,'-',3.0)
        plotimage(imgfit_pl,zminfl,zmaxfl,3,row,column,xdim,ydim,0.07,0.08,'fit',colmap,labcol)
        plotimage(imgres_pl,zminfl,zmaxfl,4,row,column,xdim,ydim,0.44,0.08,'residual',colmap,labcol)

# plot data color bar

    barwin = pylab.axes([0.84,0.08,0.06,0.9])
    if imscale == 'linear':
        brange = numpy.arange(zminfl,zmaxfl,(zmaxfl-zminfl)/1000)
    elif imscale == 'logarithmic':
        brange = numpy.arange(10.0**zminfl,10.0**zmaxfl,(10.0**zmaxfl-10.0**zminfl)/1000)
    elif imscale == 'squareroot':
        brange = numpy.arange(zminfl**2,zmaxfl**2,(zmaxfl**2-zminfl**2)/1000)
    if imscale == 'linear':
        barimg = numpy.resize(brange,(1000,1))
    elif imscale == 'logarithmic':
        barimg = numpy.log10(numpy.resize(brange,(1000,1)))
    elif imscale == 'squareroot':
        barimg = numpy.sqrt(numpy.resize(brange,(1000,1)))
    try:
        nrm = len(str(int(numpy.nanmax(brange))))-1
    except:
        nrm = 0
    brange = brange / 10**nrm
    pylab.imshow(barimg,aspect='auto',interpolation='nearest',origin='lower',
                 vmin=numpy.nanmin(barimg),vmax=numpy.nanmax(barimg),
                 extent=(0.0,1.0,brange[0],brange[-1]),cmap=colmap)
    barwin.yaxis.tick_right()
    barwin.yaxis.set_label_position('right')
    barwin.yaxis.set_major_locator(MaxNLocator(7))
    pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    pylab.gca().set_autoscale_on(False)
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[])
    pylab.ylabel('Flux (10$^%d$ e$^-$ s$^{-1}$)' % nrm)
    setp(barwin.get_yticklabels(), 'rotation', 90)
    barwin.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

# render plot

    if status == 0 and len(plotfile) > 0 and plotfile.lower() != 'none':
        pylab.savefig(plotfile)
    if status == 0 and plt:
        if cmdLine:
            pylab.show(block=True)
        else:
            pylab.ion()
            pylab.plot([])
            pylab.ioff()

# stop time

    kepmsg.clock('\nKEPPRF ended at',logfile,verbose)

    return

# -----------------------------------------------------------
# plot channel image

def plotimage(imgflux_pl,zminfl,zmaxfl,plmode,row,column,xdim,ydim,winx,winy,tlabel,colmap,labcol):

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

    ax = pylab.axes([winx,winy,0.37,0.45])
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
               fontsize=36,fontweight=500,color=labcol,transform=ax.transAxes)

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
    parser.add_argument('plotfile', help='Name of output PNG plot file', type=str)
    parser.add_argument('--rownum', '-r', default=2200, help='Row number of image stored in infile', dest='rownum', type=int)
    parser.add_argument('--columns', help='Column number of each source to be fit', type=str)
    parser.add_argument('--rows', help='Row number of each source to be fit', type=str)
    parser.add_argument('--fluxes', help='Relative flux of each source to be fit', type=str)
    parser.add_argument('--border', '-b', help='Order of background polynmial fit', default=1, dest='border', type=int)
    parser.add_argument('--background', action='store_true', help='Fit background?', default=False)
    parser.add_argument('--focus', action='store_true', help='Fit focus changes?', default=False)
    parser.add_argument('--prfdir', help='Folder containing Point Response Function FITS files', type=str)
    parser.add_argument('--xtol', '-x', default=1.0e-4, help='Fit parameter tolerance', dest='xtol', type=float)
    parser.add_argument('--ftol', '-f', default=1.0, help='Fit minimization tolerance', dest='ftol', type=float)
    parser.add_argument('--imscale', '-i', help='Type of image intensity scale', default='linear', dest='imscale', type=str,choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--colmap', '-c', help='Image colormap', default='YlOrBr', dest='cmap', type=str,choices=['Accent','Blues','BrBG','BuGn','BuPu','Dark2','GnBu','Greens','Greys','OrRd','Oranges','PRGn','Paired','Pastel1','Pastel2','PiYG','PuBu','PuBuGn','PuOr','PuRd','Purples','RdBu','RdGy','RdPu','RdYlBu','RdYlGn','Reds','Set1','Set2','Set3','Spectral','YlGn','YlGnBu','YlOrBr','YlOrRd','afmhot','autumn','binary','bone','brg','bwr','cool','copper','flag','gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_yarg','gnuplot','gnuplot2','gray','hot','hsv','jet','ocean','pink','prism','rainbow','seismic','spectral','spring','summer','terrain','winter','browse'])
    parser.add_argument('--labcol', help='Label color', default='#ffffff', type=str)
    parser.add_argument('--apercol', help='Aperture color', default='#ffffff', type=str)
    parser.add_argument('--plot', action='store_true', help='Plot fit results?', default=False)
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepprfphot.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepprf(args.infile,args.plotfile,args.rownum,args.columns,args.rows,args.fluxes,args.border,
           args.background,args.focus,args.prfdir,args.xtol,args.ftol,args.imscale,args.cmap,
           args.labcol,args.apercol,args.plot,args.verbose,args.logfile,args.status,cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepprf.par")
    t = iraf.IrafTaskFactory(taskname="kepprf", value=parfile, function=kepprf)
