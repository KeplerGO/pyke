import pylab, numpy, pyfits, scipy, multiprocessing, itertools
from pylab import *
from matplotlib import *
from numpy import *
from pyfits import *
import kepio, kepmsg, kepkey, kepplot, kepfit, keparray, kepfunc
import sys, time, re, math, glob
from scipy import stats, interpolate, optimize, ndimage
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline

# -----------------------------------------------------------
# core code

def kepprfphot(infile,outroot,columns,rows,fluxes,border,background,focus,prfdir,ranges,
               tolerance,ftolerance,qualflags,plt,clobber,verbose,logfile,status,cmdLine=False):

# input arguments

    status = 0
    seterr(all="ignore")

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPPRFPHOT -- '
    call += 'infile='+infile+' '
    call += 'outroot='+outroot+' '
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
    call += 'ranges='+ranges+' '
    call += 'xtol='+str(tolerance)+' '
    call += 'ftol='+str(ftolerance)+' '
    quality = 'n'
    if (qualflags): quality = 'y'
    call += 'qualflags='+quality+' '
    plotit = 'n'
    if (plt): plotit = 'y'
    call += 'plot='+plotit+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# test log file

    logfile = kepmsg.test(logfile)

# start time

    kepmsg.clock('KEPPRFPHOT started at',logfile,verbose)

# number of sources

    if status == 0:
        work = fluxes.strip()
        work = re.sub(' ',',',work)
        work = re.sub(';',',',work)
        nsrc = len(work.split(','))

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

# clobber output file

    for i in range(nsrc):
        outfile = '%s_%d.fits' % (outroot, i)
        if clobber: status = kepio.clobber(outfile,logfile,verbose)
        if kepio.fileexists(outfile):
            message = 'ERROR -- KEPPRFPHOT: ' + outfile + ' exists. Use --clobber'
            status = kepmsg.err(logfile,message,verbose)

# open TPF FITS file

    if status == 0:
        try:
            kepid, channel, skygroup, module, output, quarter, season, \
                ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
                kepio.readTPF(infile,'TIME',logfile,verbose)
        except:
            message = 'ERROR -- KEPPRFPHOT: is %s a Target Pixel File? ' % infile
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
            ra, dec, column, row, kepmag, xdim, ydim, poscorr1, status = \
            kepio.readTPF(infile,'POS_CORR1',logfile,verbose)
        if status != 0:
            poscorr1 = numpy.zeros((len(barytime)),dtype='float32')
            poscorr1[:] = numpy.nan
            status = 0
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, poscorr2, status = \
            kepio.readTPF(infile,'POS_CORR2',logfile,verbose)
        if status != 0:
            poscorr2 = numpy.zeros((len(barytime)),dtype='float32')
            poscorr2[:] = numpy.nan
            status = 0
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)
    if status == 0:
        struct, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(struct,infile,logfile,verbose,status)

# input file keywords and mask map

    if status == 0:
        cards0 = struct[0].header.cards
        cards1 = struct[1].header.cards
        cards2 = struct[2].header.cards
        maskmap = copy(struct[2].data)
        npix = numpy.size(numpy.nonzero(maskmap)[0])

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
            message = 'ERROR -- KEPPRFPHOT: No PRF file found in ' + prfdir
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
                prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
        prf = prf / nansum(prf)
        prf = prf / cdelt1p[0] / cdelt2p[0]

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfDimY = ydim / cdelt1p[0]
        prfDimX = xdim / cdelt2p[0]
        PRFy0 = (shape(prf)[0] - prfDimY) / 2
        PRFx0 = (shape(prf)[1] - prfDimX) / 2

# construct input pixel image

    if status == 0:
        DATx = arange(column,column+xdim)
        DATy = arange(row,row+ydim)

# interpolation function over the PRF

    if status == 0:
        splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf,kx=3,ky=3)

# construct mesh for background model

    if status == 0:
        bx = numpy.arange(1.,float(xdim+1))
        by = numpy.arange(1.,float(ydim+1))
        xx, yy = numpy.meshgrid(numpy.linspace(bx.min(), bx.max(), xdim),
                                numpy.linspace(by.min(), by.max(), ydim))

# Get time ranges for new photometry, flag good data

    if status == 0:
        barytime += bjdref
        tstart,tstop,status = kepio.timeranges(ranges,logfile,verbose)
        incl = numpy.zeros((len(barytime)),dtype='int')
        for rownum in xrange(len(barytime)):
            for winnum in xrange(len(tstart)):
                if barytime[rownum] >= tstart[winnum] and \
                        barytime[rownum] <= tstop[winnum] and \
                        (qual[rownum] == 0 or qualflags) and \
                        numpy.isfinite(barytime[rownum]) and \
                        numpy.isfinite(numpy.nansum(fluxpixels[rownum,:])):
                    incl[rownum] = 1
        if not numpy.in1d(1,incl):
            message = 'ERROR -- KEPPRFPHOT: No legal data within the range ' + ranges
            status = kepmsg.err(logfile,message,verbose)

# filter out bad data

    if status == 0:
        n = 0
        nincl = (incl == 1).sum()
        tim = zeros((nincl),'float64')
        tco = zeros((nincl),'float32')
        cad = zeros((nincl),'float32')
        flu = zeros((nincl,len(fluxpixels[0])),'float32')
        fer = zeros((nincl,len(fluxpixels[0])),'float32')
        pc1 = zeros((nincl),'float32')
        pc2 = zeros((nincl),'float32')
        qua = zeros((nincl),'float32')
        for rownum in xrange(len(barytime)):
            if incl[rownum] == 1:
                tim[n] = barytime[rownum]
                tco[n] = tcorr[rownum]
                cad[n] = cadno[rownum]
                flu[n,:] = fluxpixels[rownum]
                fer[n,:] = errpixels[rownum]
                pc1[n] = poscorr1[rownum]
                pc2[n] = poscorr2[rownum]
                qua[n] = qual[rownum]
                n += 1
        barytime = tim * 1.0
        tcorr = tco * 1.0
        cadno = cad * 1.0
        fluxpixels = flu * 1.0
        errpixels = fer * 1.0
        poscorr1 = pc1 * 1.0
        poscorr2 = pc2 * 1.0
        qual = qua * 1.0

# initialize plot arrays

    if status == 0:
        t = numpy.array([],dtype='float64')
        fl = []; dx = []; dy = []; bg = []; fx = []; fy = []; fa = []; rs = []; ch = []
        for i in range(nsrc):
            fl.append(numpy.array([],dtype='float32'))
            dx.append(numpy.array([],dtype='float32'))
            dy.append(numpy.array([],dtype='float32'))

# Preparing fit data message

    if status == 0:
        progress = numpy.arange(nincl)
        if verbose:
            txt  = 'Preparing...'
            sys.stdout.write(txt)
            sys.stdout.flush()

# single processor version

    if status == 0:# and not cmdLine:
        oldtime = 0.0
        for rownum in xrange(numpy.min([80,len(barytime)])):
            try:
                if barytime[rownum] - oldtime > 0.5:
                    ftol = 1.0e-10; xtol = 1.0e-10
            except:
                pass
            args = (fluxpixels[rownum,:],errpixels[rownum,:],DATx,DATy,nsrc,border,xx,yy,PRFx,PRFy,splineInterpolation,
                    guess,ftol,xtol,focus,background,rownum,80,float(x[i]),float(y[i]),False)
            guess = PRFfits(args)
            ftol = ftolerance; xtol = tolerance; oldtime = barytime[rownum]

# Fit the time series: multi-processing

    if status == 0 and cmdLine:
        anslist = []
        cad1 = 0; cad2 = 50
        for i in range(int(nincl/50) + 1):
            try:
                fluxp = fluxpixels[cad1:cad2,:]
                errp = errpixels[cad1:cad2,:]
                progress = numpy.arange(cad1,cad2)
            except:
                fluxp = fluxpixels[cad1:nincl,:]
                errp = errpixels[cad1:nincl,:]
                progress = numpy.arange(cad1,nincl)
            try:
                args = itertools.izip(fluxp,errp,itertools.repeat(DATx),itertools.repeat(DATy),
                                      itertools.repeat(nsrc),itertools.repeat(border),itertools.repeat(xx),
                                      itertools.repeat(yy),itertools.repeat(PRFx),itertools.repeat(PRFy),
                                      itertools.repeat(splineInterpolation),itertools.repeat(guess),
                                      itertools.repeat(ftolerance),itertools.repeat(tolerance),
                                      itertools.repeat(focus),itertools.repeat(background),progress,
                                      itertools.repeat(numpy.arange(cad1,nincl)[-1]),
                                      itertools.repeat(float(x[0])),
                                      itertools.repeat(float(y[0])),itertools.repeat(True))
                p = multiprocessing.Pool()
                model = [0.0]
                model = p.imap(PRFfits,args,chunksize=1)
                p.close()
                p.join()
                cad1 += 50; cad2 += 50
                ans = array([array(item) for item in zip(*model)])
                try:
                    anslist = numpy.concatenate((anslist,ans.transpose()),axis=0)
                except:
                    anslist = ans.transpose()
                guess = anslist[-1]
                ans = anslist.transpose()
            except:
                pass

# single processor version

    if status == 0 and not cmdLine:
        oldtime = 0.0; ans = []
#        for rownum in xrange(1,10):
        for rownum in xrange(nincl):
            proctime = time.time()
            try:
                if barytime[rownum] - oldtime > 0.5:
                    ftol = 1.0e-10; xtol = 1.0e-10
            except:
                pass
            args = (fluxpixels[rownum,:],errpixels[rownum,:],DATx,DATy,nsrc,border,xx,yy,PRFx,PRFy,splineInterpolation,
                    guess,ftol,xtol,focus,background,rownum,nincl,float(x[0]),float(y[0]),True)
            guess = PRFfits(args)
            ans.append(guess)
            ftol = ftolerance; xtol = tolerance; oldtime = barytime[rownum]
        ans = array(ans).transpose()

# unpack the best fit parameters

    if status == 0:
        flux = []; OBJx = []; OBJy = []
        na = shape(ans)[1]
        for i in range(nsrc):
            flux.append(ans[i,:])
            OBJx.append(ans[nsrc+i,:])
            OBJy.append(ans[nsrc*2+i,:])
        try:
            bterms = border + 1
            if bterms == 1:
                b = ans[nsrc*3,:]
            else:
                b = array([])
                bkg = []
                for i in range(na):
                    bcoeff = array([ans[nsrc*3:nsrc*3+bterms,i],ans[nsrc*3+bterms:nsrc*3+bterms*2,i]])
                    bkg.append(kepfunc.polyval2d(xx,yy,bcoeff))
                    b = numpy.append(b,nanmean(bkg[-1].reshape(bkg[-1].size)))
        except:
            b = zeros((na))
        if focus:
            wx = ans[-3,:]; wy = ans[-2,:]; angle = ans[-1,:]
        else:
            wx = ones((na)); wy = ones((na)); angle = zeros((na))

# constuct model PRF in detector coordinates

    if status == 0:
        residual = []; chi2 = []
        for i in range(na):
            f = empty((nsrc))
            x = empty((nsrc))
            y = empty((nsrc))
            for j in range(nsrc):
                f[j] = flux[j][i]
                x[j] = OBJx[j][i]
                y[j] = OBJy[j][i]
            PRFfit = kepfunc.PRF2DET(f,x,y,DATx,DATy,wx[i],wy[i],angle[i],splineInterpolation)
            if background and bterms == 1:
                PRFfit = PRFfit + b[i]
            if background and bterms > 1:
                PRFfit = PRFfit + bkg[i]

# calculate residual of DATA - FIT

            xdim = shape(xx)[1]
            ydim = shape(yy)[0]
            DATimg = numpy.empty((ydim,xdim))
            n = 0
            for k in range(ydim):
                for j in range(xdim):
                    DATimg[k,j] = fluxpixels[i,n]
                    n += 1
            PRFres = DATimg - PRFfit
            residual.append(numpy.nansum(PRFres) / npix)

# calculate the sum squared difference between data and model

            chi2.append(abs(numpy.nansum(numpy.square(DATimg - PRFfit) / PRFfit)))

# load the output arrays

    if status == 0:
        otime = barytime - bjdref
        otimecorr = tcorr
        ocadenceno = cadno
        opos_corr1 = poscorr1
        opos_corr2 = poscorr2
        oquality = qual
        opsf_bkg = b
        opsf_focus1 = wx
        opsf_focus2 = wy
        opsf_rotation = angle
        opsf_residual = residual
        opsf_chi2 = chi2
        opsf_flux_err = numpy.empty((na)); opsf_flux_err.fill(numpy.nan)
        opsf_centr1_err = numpy.empty((na)); opsf_centr1_err.fill(numpy.nan)
        opsf_centr2_err = numpy.empty((na)); opsf_centr2_err.fill(numpy.nan)
        opsf_bkg_err = numpy.empty((na)); opsf_bkg_err.fill(numpy.nan)
        opsf_flux = []
        opsf_centr1 = []
        opsf_centr2 = []
        for i in range(nsrc):
            opsf_flux.append(flux[i])
            opsf_centr1.append(OBJx[i])
            opsf_centr2.append(OBJy[i])

# load the plot arrays

    if status == 0:
        t = barytime
        for i in range(nsrc):
            fl[i] = flux[i]
            dx[i] = OBJx[i]
            dy[i] = OBJy[i]
        bg = b
        fx = wx
        fy = wy
        fa = angle
        rs = residual
        ch = chi2

# construct output primary extension

    if status == 0:
        for j in range(nsrc):
            hdu0 = pyfits.PrimaryHDU()
            for i in range(len(cards0)):
                if cards0[i].key not in hdu0.header.keys():
                    hdu0.header.update(cards0[i].key, cards0[i].value, cards0[i].comment)
                else:
                    hdu0.header.cards[cards0[i].key].comment = cards0[i].comment
            status = kepkey.history(call,hdu0,outfile,logfile,verbose)
            outstr = HDUList(hdu0)

# construct output light curve extension

            col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=otime)
            col2 = Column(name='TIMECORR',format='E',unit='d',array=otimecorr)
            col3 = Column(name='CADENCENO',format='J',array=ocadenceno)
            col4 = Column(name='PSF_FLUX',format='E',unit='e-/s',array=opsf_flux[j])
            col5 = Column(name='PSF_FLUX_ERR',format='E',unit='e-/s',array=opsf_flux_err)
            col6 = Column(name='PSF_BKG',format='E',unit='e-/s/pix',array=opsf_bkg)
            col7 = Column(name='PSF_BKG_ERR',format='E',unit='e-/s',array=opsf_bkg_err)
            col8 = Column(name='PSF_CENTR1',format='E',unit='pixel',array=opsf_centr1[j])
            col9 = Column(name='PSF_CENTR1_ERR',format='E',unit='pixel',array=opsf_centr1_err)
            col10 = Column(name='PSF_CENTR2',format='E',unit='pixel',array=opsf_centr2[j])
            col11 = Column(name='PSF_CENTR2_ERR',format='E',unit='pixel',array=opsf_centr2_err)
            col12 = Column(name='PSF_FOCUS1',format='E',array=opsf_focus1)
            col13 = Column(name='PSF_FOCUS2',format='E',array=opsf_focus2)
            col14 = Column(name='PSF_ROTATION',format='E',unit='deg',array=opsf_rotation)
            col15 = Column(name='PSF_RESIDUAL',format='E',unit='e-/s',array=opsf_residual)
            col16 = Column(name='PSF_CHI2',format='E',array=opsf_chi2)
            col17 = Column(name='POS_CORR1',format='E',unit='pixel',array=opos_corr1)
            col18 = Column(name='POS_CORR2',format='E',unit='pixel',array=opos_corr2)
            col19 = Column(name='SAP_QUALITY',format='J',array=oquality)
            cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,
                            col12,col13,col14,col15,col16,col17,col18,col19])
            hdu1 = new_table(cols)
            for i in range(len(cards1)):
                if (cards1[i].key not in hdu1.header.keys() and
                    cards1[i].key[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                              '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                              '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                              '12PC','21PC','22PC']):
                    hdu1.header.update(cards1[i].key, cards1[i].value, cards1[i].comment)
            outstr.append(hdu1)

# construct output mask bitmap extension

            hdu2 = ImageHDU(maskmap)
            for i in range(len(cards2)):
                if cards2[i].key not in hdu2.header.keys():
                    hdu2.header.update(cards2[i].key, cards2[i].value, cards2[i].comment)
                else:
                    hdu2.header.cards[cards2[i].key].comment = cards2[i].comment
            outstr.append(hdu2)

# write output file

            outstr.writeto(outroot + '_' + str(j) + '.fits',checksum=True)

# close input structure

            status = kepio.closefits(struct,logfile,verbose)

# clean up x-axis unit

    if status == 0:
	barytime0 = float(int(t[0] / 100) * 100.0)
	t -= barytime0
        t = numpy.insert(t,[0],[t[0]])
        t = numpy.append(t,[t[-1]])
        xlab = 'BJD $-$ %d' % barytime0

# plot the light curves

    if status == 0:
        bg = numpy.insert(bg,[0],[-1.0e10])
        bg = numpy.append(bg,-1.0e10)
        fx = numpy.insert(fx,[0],[fx[0]])
        fx = numpy.append(fx,fx[-1])
        fy = numpy.insert(fy,[0],[fy[0]])
        fy = numpy.append(fy,fy[-1])
        fa = numpy.insert(fa,[0],[fa[0]])
        fa = numpy.append(fa,fa[-1])
        rs = numpy.insert(rs,[0],[-1.0e10])
        rs = numpy.append(rs,-1.0e10)
        ch = numpy.insert(ch,[0],[-1.0e10])
        ch = numpy.append(ch,-1.0e10)
        for i in range(nsrc):

# clean up y-axis units

            nrm = math.ceil(math.log10(numpy.nanmax(fl[i]))) - 1.0
            fl[i] /= 10**nrm
            if nrm == 0:
                ylab1 = 'e$^-$ s$^{-1}$'
            else:
                ylab1 = '10$^{%d}$ e$^-$ s$^{-1}$' % nrm
            xx = copy(dx[i])
            yy = copy(dy[i])
            ylab2 = 'offset (pixels)'

# data limits

            xmin = numpy.nanmin(t)
            xmax = numpy.nanmax(t)
            ymin1 = numpy.nanmin(fl[i])
            ymax1 = numpy.nanmax(fl[i])
            ymin2 = numpy.nanmin(xx)
            ymax2 = numpy.nanmax(xx)
            ymin3 = numpy.nanmin(yy)
            ymax3 = numpy.nanmax(yy)
            ymin4 = numpy.nanmin(bg[1:-1])
            ymax4 = numpy.nanmax(bg[1:-1])
            ymin5 = numpy.nanmin([numpy.nanmin(fx),numpy.nanmin(fy)])
            ymax5 = numpy.nanmax([numpy.nanmax(fx),numpy.nanmax(fy)])
            ymin6 = numpy.nanmin(fa[1:-1])
            ymax6 = numpy.nanmax(fa[1:-1])
            ymin7 = numpy.nanmin(rs[1:-1])
            ymax7 = numpy.nanmax(rs[1:-1])
            ymin8 = numpy.nanmin(ch[1:-1])
            ymax8 = numpy.nanmax(ch[1:-1])
            xr = xmax - xmin
            yr1 = ymax1 - ymin1
            yr2 = ymax2 - ymin2
            yr3 = ymax3 - ymin3
            yr4 = ymax4 - ymin4
            yr5 = ymax5 - ymin5
            yr6 = ymax6 - ymin6
            yr7 = ymax7 - ymin7
            yr8 = ymax8 - ymin8
            fl[i] = numpy.insert(fl[i],[0],[0.0])
            fl[i] = numpy.append(fl[i],0.0)

# plot style

            try:
                params = {'backend': 'png',
                          'axes.linewidth': 2.5,
                          'axes.labelsize': 24,
                          'axes.font': 'sans-serif',
                          'axes.fontweight' : 'bold',
                          'text.fontsize': 12,
                          'legend.fontsize': 12,
                          'xtick.labelsize': 12,
                          'ytick.labelsize': 12}
                pylab.rcParams.update(params)
            except:
                pass

# define size of plot on monitor screen

            pylab.figure(str(i+1) + ' ' + str(time.asctime(time.localtime())),figsize=[12,16])

# delete any fossil plots in the matplotlib window

            pylab.clf()

# position first axes inside the plotting window

            ax = pylab.axes([0.11,0.523,0.78,0.45])

# force tick labels to be absolute rather than relative

            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# no x-label

            pylab.setp(pylab.gca(),xticklabels=[])

# plot flux vs time

            ltime = numpy.array([],dtype='float64')
            ldata = numpy.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = numpy.append(ltime,t[j])
                    ldata = numpy.append(ldata,fl[i][j])
                else:
                    pylab.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)
                    ltime = numpy.array([],dtype='float64')
                    ldata = numpy.array([],dtype='float32')
            pylab.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)

# plot the fill color below data time series, with no data gaps

            pylab.fill(t,fl[i],fc='#ffff00',linewidth=0.0,alpha=0.2)

# define plot x and y limits

            pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            if ymin1 - yr1 * 0.01 <= 0.0:
                pylab.ylim(1.0e-10, ymax1 + yr1 * 0.01)
            else:
                pylab.ylim(ymin1 - yr1 * 0.01, ymax1 + yr1 * 0.01)

# plot labels

#            pylab.xlabel(xlab, {'color' : 'k'})
            try:
                pylab.ylabel('Source (' + ylab1 + ')', {'color' : 'k'})
            except:
                ylab1 = '10**%d e-/s' % nrm
                pylab.ylabel('Source (' + ylab1 + ')', {'color' : 'k'})

# make grid on plot

            pylab.grid()

# plot centroid tracks - position second axes inside the plotting window

            if focus and background:
                axs = [0.11,0.433,0.78,0.09]
            elif background or focus:
                axs = [0.11,0.388,0.78,0.135]
            else:
                axs = [0.11,0.253,0.78,0.27]
            ax1 = pylab.axes(axs)

# force tick labels to be absolute rather than relative

            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.setp(pylab.gca(),xticklabels=[])

# plot dx vs time

            ltime = numpy.array([],dtype='float64')
            ldata = numpy.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = numpy.append(ltime,t[j])
                    ldata = numpy.append(ldata,xx[j-1])
                else:
                    ax1.plot(ltime,ldata,color='r',linestyle='-',linewidth=1.0)
                    ltime = numpy.array([],dtype='float64')
                    ldata = numpy.array([],dtype='float32')
            ax1.plot(ltime,ldata,color='r',linestyle='-',linewidth=1.0)

# define plot x and y limits

            pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            pylab.ylim(ymin2 - yr2 * 0.03, ymax2 + yr2 * 0.03)

# plot labels

            ax1.set_ylabel('X-' + ylab2, color='k', fontsize=11)

# position second axes inside the plotting window

            ax2 = ax1.twinx()

# force tick labels to be absolute rather than relative

            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.setp(pylab.gca(),xticklabels=[])

# plot dy vs time

            ltime = numpy.array([],dtype='float64')
            ldata = numpy.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = numpy.append(ltime,t[j])
                    ldata = numpy.append(ldata,yy[j-1])
                else:
                    ax2.plot(ltime,ldata,color='g',linestyle='-',linewidth=1.0)
                    ltime = numpy.array([],dtype='float64')
                    ldata = numpy.array([],dtype='float32')
            ax2.plot(ltime,ldata,color='g',linestyle='-',linewidth=1.0)

# define plot y limits

            pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            pylab.ylim(ymin3 - yr3 * 0.03, ymax3 + yr3 * 0.03)

# plot labels

            ax2.set_ylabel('Y-' + ylab2, color='k',fontsize=11)

# background - position third axes inside the plotting window

            if background and focus:
                axs = [0.11,0.343,0.78,0.09]
            if background and not focus:
                axs = [0.11,0.253,0.78,0.135]
            if background:
                ax1 = pylab.axes(axs)

# force tick labels to be absolute rather than relative

                pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.setp(pylab.gca(),xticklabels=[])

# plot background vs time

                ltime = numpy.array([],dtype='float64')
                ldata = numpy.array([],dtype='float32')
                dt = 0
                work1 = 2.0 * cadence / 86400
                for j in range(1,len(t)-1):
                    dt = t[j] - t[j-1]
                    if dt < work1:
                        ltime = numpy.append(ltime,t[j])
                        ldata = numpy.append(ldata,bg[j])
                    else:
                        ax1.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)
                        ltime = numpy.array([],dtype='float64')
                        ldata = numpy.array([],dtype='float32')
                ax1.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)

# plot the fill color below data time series, with no data gaps

                pylab.fill(t,bg,fc='#ffff00',linewidth=0.0,alpha=0.2)

# define plot x and y limits

                pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
                pylab.ylim(ymin4 - yr4 * 0.03, ymax4 + yr4 * 0.03)

# plot labels

                ax1.set_ylabel('Background \n(e$^-$ s$^{-1}$ pix$^{-1}$)',
                               multialignment='center', color='k',fontsize=11)

# make grid on plot

                pylab.grid()

# position focus axes inside the plotting window

            if focus and background:
                axs = [0.11,0.253,0.78,0.09]
            if focus and not background:
                axs = [0.11,0.253,0.78,0.135]
            if focus:
                ax1 = pylab.axes(axs)

# force tick labels to be absolute rather than relative

                pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.setp(pylab.gca(),xticklabels=[])

# plot x-axis PSF width vs time

                ltime = numpy.array([],dtype='float64')
                ldata = numpy.array([],dtype='float32')
                dt = 0
                work1 = 2.0 * cadence / 86400
                for j in range(1,len(t)-1):
                    dt = t[j] - t[j-1]
                    if dt < work1:
                        ltime = numpy.append(ltime,t[j])
                        ldata = numpy.append(ldata,fx[j])
                    else:
                        ax1.plot(ltime,ldata,color='r',linestyle='-',linewidth=1.0)
                        ltime = numpy.array([],dtype='float64')
                        ldata = numpy.array([],dtype='float32')
                ax1.plot(ltime,ldata,color='r',linestyle='-',linewidth=1.0)

# plot y-axis PSF width vs time

                ltime = numpy.array([],dtype='float64')
                ldata = numpy.array([],dtype='float32')
                dt = 0
                work1 = 2.0 * cadence / 86400
                for j in range(1,len(t)-1):
                    dt = t[j] - t[j-1]
                    if dt < work1:
                        ltime = numpy.append(ltime,t[j])
                        ldata = numpy.append(ldata,fy[j])
                    else:
                        ax1.plot(ltime,ldata,color='g',linestyle='-',linewidth=1.0)
                        ltime = numpy.array([],dtype='float64')
                        ldata = numpy.array([],dtype='float32')
                ax1.plot(ltime,ldata,color='g',linestyle='-',linewidth=1.0)

# define plot x and y limits

                pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
                pylab.ylim(ymin5 - yr5 * 0.03, ymax5 + yr5 * 0.03)

# plot labels

                ax1.set_ylabel('Pixel Scale\nFactor',
                               multialignment='center', color='k',fontsize=11)

# Focus rotation - position second axes inside the plotting window

                ax2 = ax1.twinx()

# force tick labels to be absolute rather than relative

                pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
                pylab.setp(pylab.gca(),xticklabels=[])

# plot dy vs time

                ltime = numpy.array([],dtype='float64')
                ldata = numpy.array([],dtype='float32')
                dt = 0
                work1 = 2.0 * cadence / 86400
                for j in range(1,len(t)-1):
                    dt = t[j] - t[j-1]
                    if dt < work1:
                        ltime = numpy.append(ltime,t[j])
                        ldata = numpy.append(ldata,fa[j])
                    else:
                        ax2.plot(ltime,ldata,color='#000080',linestyle='-',linewidth=1.0)
                        ltime = numpy.array([],dtype='float64')
                        ldata = numpy.array([],dtype='float32')
                ax2.plot(ltime,ldata,color='#000080',linestyle='-',linewidth=1.0)

# define plot y limits

                pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
                pylab.ylim(ymin6 - yr6 * 0.03, ymax6 + yr6 * 0.03)

# plot labels

                ax2.set_ylabel('Rotation (deg)', color='k',fontsize=11)

# fit residuals - position fifth axes inside the plotting window

            axs = [0.11,0.163,0.78,0.09]
            ax1 = pylab.axes(axs)

# force tick labels to be absolute rather than relative

            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.setp(pylab.gca(),xticklabels=[])

# plot residual vs time

            ltime = numpy.array([],dtype='float64')
            ldata = numpy.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = numpy.append(ltime,t[j])
                    ldata = numpy.append(ldata,rs[j])
                else:
                    ax1.plot(ltime,ldata,color='b',linestyle='-',linewidth=1.0)
                    ltime = numpy.array([],dtype='float64')
                    ldata = numpy.array([],dtype='float32')
            ax1.plot(ltime,ldata,color='b',linestyle='-',linewidth=1.0)

# plot the fill color below data time series, with no data gaps

            pylab.fill(t,rs,fc='#ffff00',linewidth=0.0,alpha=0.2)

# define plot x and y limits

            pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            pylab.ylim(ymin7 - yr7 * 0.03, ymax7 + yr7 * 0.03)

# plot labels

            ax1.set_ylabel('Residual \n(e$^-$ s$^{-1}$)',
                           multialignment='center', color='k',fontsize=11)

# make grid on plot

            pylab.grid()

# fit chi square - position sixth axes inside the plotting window

            axs = [0.11,0.073,0.78,0.09]
            ax1 = pylab.axes(axs)

# force tick labels to be absolute rather than relative

            pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# plot background vs time

            ltime = numpy.array([],dtype='float64')
            ldata = numpy.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for j in range(1,len(t)-1):
                dt = t[j] - t[j-1]
                if dt < work1:
                    ltime = numpy.append(ltime,t[j])
                    ldata = numpy.append(ldata,ch[j])
                else:
                    ax1.plot(ltime,ldata,color='b',linestyle='-',linewidth=1.0)
                    ltime = numpy.array([],dtype='float64')
                    ldata = numpy.array([],dtype='float32')
            ax1.plot(ltime,ldata,color='b',linestyle='-',linewidth=1.0)

# plot the fill color below data time series, with no data gaps

            pylab.fill(t,ch,fc='#ffff00',linewidth=0.0,alpha=0.2)

# define plot x and y limits

            pylab.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            pylab.ylim(ymin8 - yr8 * 0.03, ymax8 + yr8 * 0.03)

# plot labels

            ax1.set_ylabel('$\chi^2$ (%d dof)' % (npix-len(guess)-1),color='k',fontsize=11)
            pylab.xlabel(xlab, {'color' : 'k'})

# make grid on plot

            pylab.grid()

# render plot

            if status == 0:
                pylab.savefig(outroot + '_' + str(i) + '.png')
            if status == 0 and plt:
                if cmdLine:
                    pylab.show(block=True)
                else:
                    pylab.ion()
                    pylab.plot([])
                    pylab.ioff()

# stop time

    kepmsg.clock('\n\nKEPPRFPHOT ended at',logfile,verbose)

    return



def PRFfits(args):

# start time

    proctime = time.time()

# extract image from the time series

    xdim = shape(args[6])[1]
    ydim = shape(args[6])[0]
    DATimg = numpy.empty((ydim,xdim))
    DATerr = numpy.empty((ydim,xdim))
    n = 0
    for i in range(ydim):
        for j in range(xdim):
            DATimg[i,j] = args[0][n]
            DATerr[i,j] = args[1][n]
            n += 1

# minimize data and model

    if args[14] and args[15]:
        argm = (args[2],args[3],DATimg,DATerr,args[4],args[5],args[6],args[7],args[10],args[18],args[19])
        ans = fmin_powell(kepfunc.PRFwithFocusAndBackground,args[11],args=argm,xtol=args[12],
                          ftol=args[13],disp=False)
    elif args[14] and not args[15]:
        argm = (args[2],args[3],DATimg,DATerr,args[4],args[10],args[18],args[19])
        ans = fmin_powell(kepfunc.PRFwithFocus,args[11],args=argm,xtol=args[12],
                          ftol=args[13],disp=False)
    elif args[15] and not args[14]:
        argm = (args[2],args[3],DATimg,DATerr,args[4],args[5],args[6],args[7],args[10],args[18],args[19])
        ans = fmin_powell(kepfunc.PRFwithBackground,args[11],args=argm,xtol=args[12],
                          ftol=args[13],disp=False)
    else:
        argm = (args[2],args[3],DATimg,DATerr,args[4],args[10],args[18],args[19])
        ans = fmin_powell(kepfunc.PRF,args[11],args=argm,xtol=args[12],
                          ftol=args[13],disp=False)

# print progress

    if args[20]:
        txt  = '\r%3d%% ' % ((float(args[16]) + 1.0) / float(args[17]) * 100.0)
        txt += 'nrow = %d ' % (args[16]+1)
        txt += 't = %.1f sec' % (time.time() - proctime)
        txt += ' ' * 5
        sys.stdout.write(txt)
        sys.stdout.flush()

    return ans

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Fitting PRF model to Target Pixel time series')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('outroot', help='Root name of output light curve files', type=str)
    parser.add_argument('--columns', help='Column number of each source to be fit', type=str)
    parser.add_argument('--rows', help='Row number of each source to be fit', type=str)
    parser.add_argument('--fluxes', help='Relative flux of each source to be fit', type=str)
    parser.add_argument('--border', '-b', help='Order of background polynmial fit', default=0, dest='border', type=int)
    parser.add_argument('--background', action='store_true', help='Fit background?', default=False)
    parser.add_argument('--focus', action='store_true', help='Fit focus changes?', default=False)
    parser.add_argument('--prfdir', default='/Volumes/data/Kepler/PRF', help='Folder containing Point Response Function FITS files', dest='prfdir', type=str)
    parser.add_argument('--ranges', default='0,0', help='Time ranges to fit', dest='ranges', type=str)
    parser.add_argument('--xtol', default=1.0e-4, help='Fit parameter tolerance', dest='tolerance', type=float)
    parser.add_argument('--ftol', default=1.0e-2, help='Fit minimization tolerance', dest='ftolerance', type=float)
    parser.add_argument('--qualflags', action='store_true', help='Fit data that have quality flags?', default=False)
    parser.add_argument('--plot', action='store_true', help='Plot fit results?', default=False)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?', default=False)
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?', default=False)
    parser.add_argument('--logfile', '-l', default='kepprfphot.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepprfphot(args.infile,args.outroot,args.columns,args.rows,args.fluxes,args.border,
               args.background,args.focus,args.prfdir,args.ranges,args.tolerance,
               args.ftolerance,args.qualflags,args.plot,args.clobber,args.verbose,
               args.logfile,args.status,cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepprfphot.par")
    t = iraf.IrafTaskFactory(taskname="kepprfphot", value=parfile, function=kepprfphot)
