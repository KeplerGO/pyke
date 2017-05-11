import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy.io import fits as pyfits
import kepio, kepmsg, kepplot, kepfunc, kepstat
import sys, time, re, math, glob
from scipy import interpolate, optimize, ndimage, stats
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import interpolation

def kepprf(infile, plotfile, rownum, columns, rows, fluxes, border,
           background, focus, prfdir, xtol, ftol, imscale, labcol,
           apercol, plot, verbose, logfile):
    """
    Fit a PSF model to a specific image within a Target Pixel File.

    Fit a PSF model, combined with spacecraft jitter and pixel scale
    drift (the Pixel Response Function; PRF) to a single observation of Kepler
    target pixels.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler Target
        Pixel data within the first data extension.
    plotfile : str
        Name of an optional output plot file containing the results of kepprf.
        An example is provided in Figure 1. Typically this is a PNG format
        file. If no file is required, plotfile can be 'None' or blank, in
        which case the plot will be generated but the plot will not be saved
        to a file. Any existing file with this name will be automatically
        overwritten.
    rownum : int
        The row number in the input file data table containing the pixels to
        plot. If the chosen observation has a non-zero quality flag set or the
        pixel set contains only NULLs then the task will halt with an error
        message.
    columns : float
        A starting guess for the CCD column position(s) of the source(s) that
        are to be fit. The model is unlikely to converge if the guess is too
        far away from the correct location. A rule of thumb is to provide a
        guess within 1 CCD pixel of the true position. If more than one source
        is being modeled then the column positions of each are separated by a
        comma. The same number of sources in the columns, rows and fluxes
        field is a requirement of this task.
    rows : float
        A starting guess for the CCD row position(s) of the source(s) that are
        to be fit. The model is unlikely to converge if the guess is too far
        away from the correct location. A rule of thumb is to provide a guess
        within 1 CCD pixel of the true position. If more than one source is
        being modeled then the row positions of each are separated by a comma.
        The same number of sources in the columns, rows and fluxes field is a
        requirement of this task.
    fluxes : float
        A starting guess for the flux(es) of the source(s) that are to be fit.
        Fit convergence is not particularly reliant on the accuracy of these
        guesses, but the fit will converge faster the more accurate the guess.
        If more than one source is being modeled then the row positions of
        each are separated by a comma. The same number of sources in the
        columns, rows and fluxes field is a requirement of this task.
    border : int (optional)
        If a background is included in the fit then it is modeled as a
        two-dimensional polynomial. This parameter is the polynomial order. A
        zero-order polynomial is generally recommended.
    background : boolean (optional)
        Whether to include a background component in the model. If `True` then
        the background will be represented by a two-dimensional polynomial of
        order `border`. This functionality is somewhat experimental, with one
        eye upon potential background gradients across large masks or on those
        detectors more prone to pattern noise. Generally it is recommended to
        set background as `False`.
    focus : boolean (optional)
        Whether to incude pixel scale and focus rotation with the fit
        parameters of the model. This is also an experimental function. This
        approach does not attempt to deal with inter- or intra-pixel
        variations. The recommended use is currently to set focus as `False`.
    xtol : float
        The dimensionless, relative model parameter convergence criterion for
        the fit algorithm.
    ftol : float
        The dimensionless, relative model residual convergence criterion for
        the fit algorithm.
    imscale : str
        kepprf can plot images with three choices of image scales. The choice
        is made using this argument.
        The options are:
            * linear
            * logarithmic
            * squareroot
    plot : boolean (optional)
        Plot fit results to the screen?
    verbose : boolean (optional)
        Print informative messages and warnings to the shell and logfile?
    logfile : string (optional)
        Name of the logfile containing error and warning messages.
    """

# input arguments

    status = 0
    np.seterr(all="ignore")

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
    if (plot): plotit = 'y'
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
        npix = np.size(np.nonzero(maskimg)[0])

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
        if not np.isfinite(barytime[rownum-1]) or np.nansum(fluxpixels[rownum-1,:]) == np.nan:
            message = 'ERROR -- KEPFIELD: Row ' + str(rownum) + ' is a bad quality timestamp'
            status = kepmsg.err(logfile,message,verbose)

# construct input pixel image

    if status == 0:
        flux = fluxpixels[rownum-1,:]
        ferr = errpixels[rownum-1,:]
        DATx = np.arange(column,column+xdim)
        DATy = np.arange(row,row+ydim)
#        if np.nanmin > 420000.0: flux -= 420000.0

# image scale and intensity limits of pixel data

    if status == 0:
        n = 0
        DATimg = np.empty((ydim,xdim))
        ERRimg = np.empty((ydim,xdim))
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
        crpix1p = np.zeros((5),dtype='float32')
        crpix2p = np.zeros((5),dtype='float32')
        crval1p = np.zeros((5),dtype='float32')
        crval2p = np.zeros((5),dtype='float32')
        cdelt1p = np.zeros((5),dtype='float32')
        cdelt2p = np.zeros((5),dtype='float32')
        for i in range(5):
            prfn[i], crpix1p[i], crpix2p[i], crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i], status \
                = kepio.readPRFimage(prffile,i+1,logfile,verbose)
        prfn = np.array(prfn)
        PRFx = np.arange(0.5,np.shape(prfn[0])[1]+0.5)
        PRFy = np.arange(0.5,np.shape(prfn[0])[0]+0.5)
        PRFx = (PRFx - np.size(PRFx) / 2) * cdelt1p[0]
        PRFy = (PRFy - np.size(PRFy) / 2) * cdelt2p[0]

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = np.zeros(np.shape(prfn[0]),dtype='float32')
        prfWeight = np.zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = math.sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e-6
            prf = prf + prfn[i] / prfWeight[i]
        prf = prf / np.nansum(prf) / cdelt1p[0] / cdelt2p[0]

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfDimY = int(ydim / cdelt1p[0])
        prfDimX = int(xdim / cdelt2p[0])
        PRFy0 = int(np.round((np.shape(prf)[0] - prfDimY) / 2))
        PRFx0 = int(np.round((np.shape(prf)[1] - prfDimX) / 2))

# interpolation function over the PRF

    if status == 0:
        splineInterpolation = RectBivariateSpline(PRFx,PRFy,prf)

# construct mesh for background model

    if status == 0 and background:
        bx = np.arange(1.,float(xdim+1))
        by = np.arange(1.,float(ydim+1))
        xx, yy = np.meshgrid(np.linspace(bx.min(), bx.max(), xdim),
                                np.linspace(by.min(), by.max(), ydim))

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
        PRFmod = np.zeros((prfDimY,prfDimX))
        if PRFy0 < 0 or PRFx0 < 0.0:
            PRFmod = np.zeros((prfDimY,prfDimX))
            superPRF = np.zeros((prfDimY+1,prfDimX+1))
            superPRF[abs(PRFy0):abs(PRFy0)+np.shape(prf)[0],abs(PRFx0):abs(PRFx0)+np.shape(prf)[1]] = prf
            prf = superPRF * 1.0
            PRFy0 = 0
            PRFx0 = 0

# rotate the PRF model around its center

        if focus:
            angle = ans[-1]
            prf = interpolation.rotate(prf,-angle,reshape=False,mode='nearest')

# iterate through the sources in the best fit PSF model

        for i in range(nsrc):
            flux.append(ans[i])
            OBJx.append(ans[nsrc+i])
            OBJy.append(ans[nsrc*2+i])

# calculate best-fit model

            y = (OBJy[i]-np.mean(DATy)) / cdelt1p[0]
            x = (OBJx[i]-np.mean(DATx)) / cdelt2p[0]
            prfTmp = interpolation.shift(prf,[y,x],order=3,mode='constant')
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
                bcoeff = np.array([ans[nsrc*3:nsrc*3+bterms],ans[nsrc*3+bterms:nsrc*3+bterms*2]])
                bkg = kepfunc.polyval2d(xx,yy,bcoeff)
                b = np.nanmean(bkg.reshape(bkg.size))
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
        FluxInMaskAll = np.nansum(PRFall)
        FluxInMaskOne = np.nansum(PRFone)
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
        FLUXres = np.nansum(PRFres) / npix

# calculate the sum squared difference between data and model

    if status == 0:
        Pearson = abs(np.nansum(np.square(DATimg - PRFfit) / PRFfit))
        Chi2 = np.nansum(np.square(DATimg - PRFfit) / np.square(ERRimg))
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
            zmaxpr = np.max(zmaxpr)
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
            plt.rcParams.update(params)
        except:
            pass
        plt.figure(figsize=[12,10])
        plt.clf()
        plotimage(imgdat_pl,zminfl,zmaxfl,1,row,column,xdim,ydim,0.07,0.53,'observation',colmap,labcol)
        plotimage(imgprf_pl,zminpr,zmaxpr,2,row,column,xdim,ydim,0.44,0.53,'model',colmap,labcol)
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,1,apercol,'--',0.5)
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,2,apercol,'-',3.0)
        plotimage(imgfit_pl,zminfl,zmaxfl,3,row,column,xdim,ydim,0.07,0.08,'fit',colmap,labcol)
        plotimage(imgres_pl,zminfl,zmaxfl,4,row,column,xdim,ydim,0.44,0.08,'residual',colmap,labcol)

# plot data color bar

    barwin = plt.axes([0.84,0.08,0.06,0.9])
    if imscale == 'linear':
        brange = np.arange(zminfl,zmaxfl,(zmaxfl-zminfl)/1000)
    elif imscale == 'logarithmic':
        brange = np.arange(10.0**zminfl,10.0**zmaxfl,(10.0**zmaxfl-10.0**zminfl)/1000)
    elif imscale == 'squareroot':
        brange = np.arange(zminfl**2,zmaxfl**2,(zmaxfl**2-zminfl**2)/1000)
    if imscale == 'linear':
        barimg = np.resize(brange,(1000,1))
    elif imscale == 'logarithmic':
        barimg = np.log10(np.resize(brange,(1000,1)))
    elif imscale == 'squareroot':
        barimg = np.sqrt(np.resize(brange,(1000,1)))
    try:
        nrm = len(str(int(np.nanmax(brange))))-1
    except:
        nrm = 0
    brange = brange / 10**nrm
    plt.imshow(barimg,aspect='auto',interpolation='nearest',origin='lower',
                 vmin=np.nanmin(barimg),vmax=np.nanmax(barimg),
                 extent=(0.0,1.0,brange[0],brange[-1]),cmap=colmap)
    barwin.yaxis.tick_right()
    barwin.yaxis.set_label_position('right')
    barwin.yaxis.set_major_locator(plt.MaxNLocator(7))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().set_autoscale_on(False)
    plt.setp(plt.gca(),xticklabels=[],xticks=[])
    plt.ylabel('Flux (10$^%d$ e$^-$ s$^{-1}$)' % nrm)
    plt.setp(barwin.get_yticklabels(), 'rotation', 90)
    barwin.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

# render plot

    if status == 0 and len(plotfile) > 0 and plotfile.lower() != 'none':
        plt.savefig(plotfile)
    if status == 0 and plot:
        plt.ion()
        plt.show()

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

    ax = plt.axes([winx,winy,0.37,0.45])
    plt.imshow(imgflux_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zminfl,vmax=zmaxfl,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    plt.gca().set_autoscale_on(False)
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    if plmode == 1:
        plt.setp(plt.gca(),xticklabels=[])
    if plmode == 2:
        plt.setp(plt.gca(),xticklabels=[],yticklabels=[])
    if plmode == 4:
        plt.setp(plt.gca(),yticklabels=[])
    if plmode == 3 or plmode == 4:
        plt.xlabel('Pixel Column Number', {'color' : 'k'})
    if plmode == 1 or plmode == 3:
        plt.ylabel('Pixel Row Number', {'color' : 'k'})
    plt.text(0.05, 0.93,tlabel,horizontalalignment='left',verticalalignment='center',
               fontsize=36,fontweight=500,color=labcol,transform=ax.transAxes)

    return
