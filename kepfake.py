import pylab, numpy, pyfits, scipy
from pylab import *
from matplotlib import *
from numpy import *
from pyfits import *
import kepmsg, kepio, kepfunc
import sys, glob
#import kepio, kepmsg, kepkey, kepplot, kepfit, keparray, kepfunc
#import sys, time, re, math, glob
#from scipy import interpolate, optimize, ndimage, stats
#from scipy.optimize import fmin_powell
#from scipy.interpolate import RectBivariateSpline
#from scipy.ndimage import interpolation
#from scipy.ndimage.interpolation import shift, rotate
#from scipy.stats import nanmean

# -----------------------------------------------------------
# core code

def kepfake(outfile,prfdir,module,output,column,row,kepmag,background,xdim,ydim,
            clobber,verbose,logfile,status,cmdLine=False): 

# input arguments

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFAKE -- '
    call += 'outfile='+outfile+' '
    call += 'prfdir='+prfdir+' '
    call += 'module='+str(module)+' '
    call += 'output='+str(output)+' '
    call += 'column='+str(column)+' '
    call += 'row='+str(row)+' '
    call += 'kepmag='+str(kepmag)+' '
    call += 'background='+str(background)+' '
    call += 'xdim='+str(xdim)+' '
    call += 'ydim='+str(ydim)+' '
    clob = 'n'
    if (clobber): clob = 'y'
    call += 'clobber='+clob+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# test log file

    logfile = kepmsg.test(logfile)

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
        pylab.figure(figsize=[12,12])
        pylab.clf()

# start time

    kepmsg.clock('KEPFAKE started at',logfile,verbose)

# clobber output file

    if status == 0:
        if clobber: status = kepio.clobber(outfile,logfile,verbose)
        if kepio.fileexists(outfile): 
            message = 'ERROR -- KEPPRFPHOT: ' + outfile + ' exists. Use --clobber'
            status = kepmsg.err(logfile,message,verbose)

# Create time-tagged arrays

    if status == 0:
        nobs = 192
        time = zeros((nobs))
        timecorr = zeros((nobs))
        cadenceno = zeros((nobs))
        raw_cnts = zeros((nobs,ydim,xdim))
        flux = zeros((nobs,ydim,xdim))
        flux_err = zeros((nobs,ydim,xdim))
        flux_bkg = zeros((nobs,ydim,xdim))
        flux_bkg_err = zeros((nobs,ydim,xdim))
        cosmic_rays = zeros((nobs,ydim,xdim))
        quality = zeros((nobs))
        pos_corr1 = zeros((nobs))
        pos_corr2 = zeros((nobs))

# Create aperture bit mqp

    if status == 0:
        aperture = ones((ydim,xdim)) * 3.0

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
            message = 'ERROR -- KEPFAKE: No PRF file found in ' + prfdir
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
            if prfWeight[i] == 0.0: prfWeight[i] = 1.0e6
            prf = prf + prfn[i] / prfWeight[i]
        prf = prf / nansum(prf) / cdelt1p[0] / cdelt2p[0]

# interpolation function over the PRF

    if status == 0:
        splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf,kx=3,ky=3)

# flux from target

    if status == 0:
        zeropoint = 33.227
        kepflux = 10.0**((zeropoint - kepmag) / 2.5)

# range of the output image in detector coordinates
        
    if status == 0:
        if xdim % 2  == 0: xdim += 1
        if ydim % 2  == 0: ydim += 1
        DATx = arange(round(column) - floor(xdim/2),round(column) + floor(xdim/2) + 1.0)
        DATy = arange(round(row) - floor(ydim/2),round(row) + floor(ydim/2) + 1.0)
    ax = pylab.axes([0.05,0.5,0.45,0.45])
    pylab.imshow(log10(prf+0.001),aspect='auto',interpolation='nearest',origin='lower',cmap='jet',
                 extent=(min(PRFx),max(PRFx),min(PRFy),max(PRFy)))
    pylab.plot([-100.,100.],[0.0,0.0],'k',ls='--')
    pylab.plot([0.0,0.0],[-100.,100.],'k',ls='--')
    pylab.xlim(min(PRFx),max(PRFx))
    pylab.ylim(min(PRFy),max(PRFy))
    ax = pylab.axes([0.5,0.5,0.45,0.45])
    TMPx = arange(min(PRFx),max(PRFx),0.1) + floor(column)
    TMPy = arange(min(PRFy),max(PRFy),0.1) + floor(row)
    PRFfit = kepfunc.PRF2DET([kepflux],[column],[row],TMPx,TMPy,1.0,1.0,0.0,splineInterpolation)
    PRFfit = PRFfit + background
    pylab.imshow(log10(PRFfit),aspect='auto',interpolation='nearest',origin='lower',cmap='jet',
                 extent=(min(TMPx),max(TMPx),min(TMPy),max(TMPy)))
    pylab.plot([column,column],[-100.,2000.],'k',ls='--')
    pylab.plot([-100.,2000.],[row,row],'k',ls='--')
    pylab.xlim(min(TMPx),max(TMPx))
    pylab.ylim(min(TMPy),max(TMPy))
    ax = pylab.axes([0.05,0.05,0.45,0.45])
    TMPx = arange(min(PRFx),max(PRFx),0.5) + floor(column)
    TMPy = arange(min(PRFy),max(PRFy),0.5) + floor(row)
    PRFfit = kepfunc.PRF2DET([kepflux],[column],[row],TMPx,TMPy,1.0,1.0,0.0,splineInterpolation)
    PRFfit = PRFfit + background
    pylab.imshow(log10(PRFfit),aspect='auto',interpolation='nearest',origin='lower',cmap='jet',
                 extent=(min(TMPx),max(TMPx),min(TMPy),max(TMPy)))
    pylab.plot([column,column],[-100.,2000.],'k',ls='--')
    pylab.plot([-100.,2000.],[row,row],'k',ls='--')
    pylab.xlim(min(TMPx),max(TMPx))
    pylab.ylim(min(TMPy),max(TMPy))
    ax = pylab.axes([0.5,0.05,0.45,0.45])
    TMPx = arange(min(PRFx),max(PRFx+1.0),1.0) + floor(column) - 0.5
    TMPy = arange(min(PRFy),max(PRFy+1.0),1.0) + floor(row) - 0.5
    PRFfit = kepfunc.PRF2DET([kepflux],[column],[row],TMPx,TMPy,1.0,1.0,0.0,splineInterpolation)
    PRFfit = PRFfit + background
    pylab.imshow(log10(PRFfit),aspect='auto',interpolation='nearest',origin='lower',cmap='jet',
                 extent=(min(TMPx)-0.5,max(TMPx)-0.5,min(TMPy)-0.5,max(TMPy)-0.5))
    pylab.plot([column,column],[-100.,2000.],'k',ls='--')
    pylab.plot([-100.,2000.],[row,row],'k',ls='--')
    pylab.xlim(min(TMPx),max(TMPx))
    pylab.ylim(min(TMPy),max(TMPy))
    pylab.ion()
    pylab.plot([])
    pylab.ioff()
#    sys.exit()

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfDimY = int(ydim / cdelt1p[0])
        prfDimX = int(xdim / cdelt2p[0])
        PRFy0 = (shape(prf)[0] - prfDimY) / 2
        PRFx0 = (shape(prf)[1] - prfDimX) / 2

# iterate over each exposure

    if status == 0:
        for i in range(nobs):
#        for i in range(1):

# constuct model PRF in detector coordinates

            if status == 0:
                PRFfit = kepfunc.PRF2DET([kepflux],[column],[row],DATx,DATy,1.0,1.0,0.0,splineInterpolation)
                PRFfit = PRFfit + background

# add noise to image

            if status == 0:
                for index,value in ndenumerate(PRFfit):
                    PRFfit[index] += sqrt(value) * kepfunc.inv_normal_cummulative_function(random.random())

# populate output array

            if status == 0:
                time[i] = 1500.2 + float(i) * 0.020416667
                cadenceno[i] = 10000 + i
                flux[i,:,:] = PRFfit
                flux_err[i,:,:] = sqrt(PRFfit)

    pylab.imshow(log10(PRFfit),aspect='auto',interpolation='nearest',origin='lower',cmap='jet',
                 extent=(min(DATx),max(DATx),min(DATy),max(DATy)))
    pylab.ion()
    pylab.plot([])
    pylab.ioff()
            
# Create the outfile primary extension

    if status == 0:
        hdu0 = PrimaryHDU()
        hdu0.header.update('EXTNAME','PRIMARY','name of extension')
        hdu0.header.update('EXTEND',True,'file may contain standard extensions')
        hdu0.header.update('EXTVER',1.0,'extension version number')
        hdu0.header.update('ORIGIN','NASA/Ames','organization that generated this file')
        hdu0.header.update('DATE','2014-04-08','file creation date')
        hdu0.header.update('CREATOR','FluxExporter-91415','SW version used to create this file')
        hdu0.header.update('PROCVER',11.0,'processing script version')
        hdu0.header.update('FILEVER','COA','file format version')
        hdu0.header.update('TIMVERSN','OGIP/93-003','OGIP memo number for file format')
        hdu0.header.update('TELESCOP','Kepler','telescope')
        hdu0.header.update('INSTRUME','Kepler photometer','detector type')
        hdu0.header.update('OBJECT','KIC 12345678','string version of kepID')
        hdu0.header.update('KEPLERID',1234567,'unique Kepler target identifier')
        hdu0.header.update('CHANNEL',58,'CCD channel')
        hdu0.header.update('SKYGROUP',32,'roll-independent location of channel')
        hdu0.header.update('MODULE',17,'CCD module')
        hdu0.header.update('OUTPUT',2,'CCD output')
        hdu0.header.update('QUARTER',4,'mission quarter during which data was collected')
        hdu0.header.update('SEASON',2,'mission season during which data was collected')
        hdu0.header.update('DATA_REL',25,'version of data release notes describing data')
        hdu0.header.update('OBSMODE','long cadence','observing mode')
        hdu0.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
        hdu0.header.update('RA_OBJ',0.0,'[deg] right ascension from KIC')
        hdu0.header.update('DEC_OBJ',0.0,'[deg] declination from KIC')
        hdu0.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
        hdu0.header.update('PMRA',0.0,'[arcsec/yr] RA proper motion')
        hdu0.header.update('PMDEC',0.0,'[arcsec/yr] Dec proper motion')
        hdu0.header.update('PMTOTAL',0.0,'[arcsec/yr] total proper motion')
        hdu0.header.update('PARALLAX',0.0,'[arcsec] parallax')
        hdu0.header.update('GLON',0.0,'[deg] galactic longitude')
        hdu0.header.update('GLAT',0.0,'[deg] galactic latitude')
        hdu0.header.update('GMAG',kepmag,'[mag] SDSS g band magnitude from KIC')
        hdu0.header.update('RMAG',kepmag,'[mag] SDSS r band magnitude from KIC')
        hdu0.header.update('IMAG',kepmag,'[mag] SDSS i band magnitude from KIC')
        hdu0.header.update('ZMAG',kepmag,'[mag] SDSS z band magnitude from KIC')
        hdu0.header.update('D51MAG',kepmag,'[mag] D51 magnitude, from KIC')
        hdu0.header.update('JMAG',kepmag,'[mag] J band magnitude from 2MASS')
        hdu0.header.update('HMAG',kepmag,'[mag] H band magnitude from 2MASS')
        hdu0.header.update('KMAG',kepmag,'[mag] K band magnitude from 2MASS')
        hdu0.header.update('KEPMAG',kepmag,'[mag] Kepler magnitude (Kp) from KIC')
        hdu0.header.update('GRCOLOR',0.0,'[mag] (g-r) color, SDSS bands')
        hdu0.header.update('JKCOLOR',0.0,'[mag] (J-K) color, 2MASS bands')
        hdu0.header.update('GKCOLOR',0.0,'[mag] (g-K) color, SDSS g - 2MASS K')
        hdu0.header.update('TEFF',5000.0,'[K] effective temperature from KIC')
        hdu0.header.update('LOGG',4.5,'[cm/s2] log10 surface gravity from KIC')
        hdu0.header.update('FEH',0.0,'[log10([Fe/H])] metallicity from KIC')
        hdu0.header.update('EBMINUSV',0.0,'[mag] E(B-V) redenning from KIC')
        hdu0.header.update('AV',0.0,'[mag] A_v extinction from KIC')
        hdu0.header.update('RADIUS',1.0,'[solar radii] stellar radius from KIC')
        hdu0.header.update('TMINDEX',1117146912,'unique 2MASS catalog ID from KIC')
        hdu0.header.update('SCPID',1117146912,'unique SCP processing ID from KIC') 
        hdulist = HDUList(hdu0)

# create the outfile table extension

    if status == 0:
        eformat = str(ydim*xdim) + 'E'
        jformat = str(ydim*xdim) + 'J'
        kformat = str(ydim*xdim) + 'K'
        coldim = '(' + str(ydim) + ',' + str(ydim) + ')'
        col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
        col2 = Column(name='TIMECORR',format='E',unit='d',array=timecorr)
        col3 = Column(name='CADENCENO',format='J',array=cadenceno)
        col4 = Column(name='RAW_CNTS',format=jformat,unit='count',dim=coldim,array=raw_cnts)
        col5 = Column(name='FLUX',format=eformat,unit='e-/s',dim=coldim,array=flux)
        col6 = Column(name='FLUX_ERR',format=eformat,unit='e-/s',dim=coldim,array=flux_err)
        col7 = Column(name='FLUX_BKG',format=eformat,unit='e-/s',dim=coldim,array=flux_bkg)
        col8 = Column(name='FLUX_BKG_ERR',format=eformat,unit='e-/s',dim=coldim,array=flux_bkg_err)
        col9 = Column(name='COSMIC_RAYS',format=eformat,unit='e-/s',dim=coldim,array=cosmic_rays)
        col10 = Column(name='QUALITY',format='J',array=quality)
        col11 = Column(name='POS_CORR1',format='E',unit='pixel',array=pos_corr1)
        col12 = Column(name='POS_CORR2',format='E',unit='pixel',array=pos_corr2)
        cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
        hdu1 = new_table(cols)
        hdu1.header.update('TTYPE1','TIME','column title: data time stamps')
        hdu1.header.update('TFORM1','D','data type: float64')
        hdu1.header.update('TUNIT1','BJD - 2454833','column units: barycenter corrected JD')
        hdu1.header.update('TDISP1','D13.7','column display format')
        hdu1.header.update('TTYPE2','TIMECORR','column title: barycentric correction')
        hdu1.header.update('TFORM2','E','column format: float32')
        hdu1.header.update('TUNIT2','d','column units: days')
        hdu1.header.update('TDISP2','D12.7','column display format')
        hdu1.header.update('TTYPE3','CADENCENO','column title: unique cadence number')
        hdu1.header.update('TFORM3','J','column format: signed integer32')
        hdu1.header.update('TTYPE4','RAW_CNTS','column title: raw pixel count')
        hdu1.header.update('TFORM4',jformat,'column format: signed integer32')
        hdu1.header.update('TDIM4',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT4','count','column units: count')
        hdu1.header.update('TTYPE5','FLUX','column title: calibrated pixel flux')
        hdu1.header.update('TFORM5',eformat,'column format: float32')
        hdu1.header.update('TDIM5',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT5','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE6','FLUX_ERR','column title: 1-sigma calibrated uncertainty')
        hdu1.header.update('TFORM6',eformat,'column format: float32')
        hdu1.header.update('TDIM6',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT6','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE7','FLUX_BKG','column title: calibrated background flux')
        hdu1.header.update('TFORM7',eformat,'column format: float32')
        hdu1.header.update('TDIM7',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT7','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE8','FLUX_BKG_ERR','column title: 1-sigma cal. backgrnd uncertainty')
        hdu1.header.update('TFORM8',eformat,'column format: float32')
        hdu1.header.update('TDIM8',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT8','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE9','COSMIC_RAYS','column title: cosmic ray detections')
        hdu1.header.update('TFORM9',eformat,'column format: float32')
        hdu1.header.update('TDIM9',coldim,'column dimensions: pixel apeture array')
        hdu1.header.update('TUNIT9','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE10','QUALITY','column title: pixel quality flags')
        hdu1.header.update('TFORM10','J','column format: signed integer32')
        hdu1.header.update('TTYPE11','POS_CORR1','column title: col correction for vel. abbern')
        hdu1.header.update('TFORM11','E','column format: float32')
        hdu1.header.update('TUNIT11','pixel','column units: pixel')
        hdu1.header.update('TTYPE12','POS_CORR2','column title: row correction for vel. abbern')
        hdu1.header.update('TFORM12','E','column format: float32')
        hdu1.header.update('TUNIT12','pixel','column units: pixel')
        hdu1.header.update('WCAX4',2,'number of WCS axes')
        hdu1.header.update('1CTY4P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY4P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP4P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP4P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV4P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV4P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN4P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN4P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE4P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE4P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP4','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP4','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX4',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX4',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL4',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL4',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI4','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI4','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT4',-0.001106815552144,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT4',0.001106815552144,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC4',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC4',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC4',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC4',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX5',2,'number of WCS axes')
        hdu1.header.update('1CTY5P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY5P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP5P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP5P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV5P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV5P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN5P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN5P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE5P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE5P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP5','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP5','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX5',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX5',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL5',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL5',43.80033,'[deg] declination at reference pixel [deg]')
        hdu1.header.update('1CUNI5','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI5','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT5',-0.001106815552144,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT5',0.001106815552144,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC5',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC5',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC5',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC5',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX6',2,'number of WCS axes')
        hdu1.header.update('1CTY6P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY6P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP6P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP6P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV6P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV6P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN6P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN6P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE6P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE6P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP6','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP6','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX6',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX6',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL6',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL6',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI6','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI6','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT6',-0.00110558987335788,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT6',0.00110558987335788,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC6',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC6',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC6',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC6',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX7',2,'number of WCS axes')
        hdu1.header.update('1CTY7P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY7P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP7P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP7P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV7P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV7P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN7P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN7P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE7P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE7P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP7','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP7','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX7',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX7',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL7',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL7',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI7','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI7','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT7',-0.00110558987335788,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT7',0.00110558987335788,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC7',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC7',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC7',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC7',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX8',2,'number of WCS axes')
        hdu1.header.update('1CTY8P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY8P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP8P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP8P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV8P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV8P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN8P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN8P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE8P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE8P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP8','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP8','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX8',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX8',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL8',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL8',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI8','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI8','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT8',-0.00110558987335788,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT8',0.00110558987335788,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC8',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC8',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC8',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC8',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX9',2,'number of WCS axes')
        hdu1.header.update('1CTY9P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY9P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP9P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP9P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV9P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV9P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN9P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN9P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE9P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE9P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP9','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP9','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX9',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX9',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL9',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL9',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI9','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI9','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT9',-0.00110558987335788,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT9',0.00110558987335788,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC9',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC9',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC9',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC9',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('WCAX10',2,'number of WCS axes')
        hdu1.header.update('1CTY10P','RAWX','right ascension coordinate type')
        hdu1.header.update('2CTY10P','RAWY','declination coordinate type')
        hdu1.header.update('1CRP10P',1,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRP10P',1,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV10P',DATx[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV10P',DATy[0],'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CUN10P','pixel','physical unit in column dimension')
        hdu1.header.update('2CUN10P','pixel','physical unit in row dimension')
        hdu1.header.update('1CDE10P',1.0,'[pixel] pixel scale in column dimension')
        hdu1.header.update('2CDE10P',1.0,'[pixel] pixel scale in row dimension')
        hdu1.header.update('1CTYP10','RA---TAN','right ascension coordinate type')
        hdu1.header.update('2CTYP10','DEC--TAN','declination coordinate type')
        hdu1.header.update('1CRPX10',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX10',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRVL10',294.94017,'[deg] right ascension at reference pixel')
        hdu1.header.update('2CRVL10',43.80033,'[deg] declination at reference pixel')
        hdu1.header.update('1CUNI10','deg','physical unit in column dimension')
        hdu1.header.update('2CUNI10','deg','physical unit in row dimension')
        hdu1.header.update('1CDLT10',-0.00110558987335788,'[deg] pixel scale in RA dimension')
        hdu1.header.update('2CDLT10',0.00110558987335788,'[deg] pixel scale in Dec dimension')
        hdu1.header.update('11PC10',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu1.header.update('12PC10',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu1.header.update('21PC10',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu1.header.update('22PC10',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdu1.header.update('INHERIT',True,'inherit primary keywords')
        hdu1.header.update('EXTNAME','TARGETTABLES','name of extension')
        hdu1.header.update('EXTVER',1,'extension version number')
        hdu1.header.update('TELESCOP','Kepler','telescope')
        hdu1.header.update('INSTRUME','Kepler photometer','detector type')
        hdu1.header.update('OBJECT','KIC 12345678','string version of kepID')
        hdu1.header.update('KEPLERID',12345678,'unique Kepler target identifier')
        hdu1.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
        hdu1.header.update('RA_OBJ',0.0,'[deg] right ascension from KIC')
        hdu1.header.update('DEC_OBJ',0.0,'[deg] declination from KIC')
        hdu1.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
        hdu1.header.update('TIMEREF','SOLARSYSTEM','barycentric correction applied to times')
        hdu1.header.update('TASSIGN','SPACECRAFT','where time is assigned')
        hdu1.header.update('TIMESYS','TDB','time system is barycentric JD')
        hdu1.header.update('BJDREFI',2454833,'integer part of BJD reference date')
        hdu1.header.update('BJDREFF',0.0,'fraction of day in BJD reference date')
        hdu1.header.update('TIMEUNIT','d','time unit for TIME, TSTART and TSTOP')
        hdu1.header.update('TSTART',1500.0,'observation start time in JD - BJDREF')
        hdu1.header.update('TSTOP',1504.0,'observation stop time in JD - BJDREF')
        hdu1.header.update('LC_START',1500.0+54833.5,'observation start time in MJD')
        hdu1.header.update('LC_END',1504.0+54833.5,'observation stop time in MJD')
        hdu1.header.update('TELAPSE',93.0,'[d] TSTOP - TSTART')
        hdu1.header.update('LIVETIME',82.7273,'[d] TELAPSE multiplied by DEADC')
        hdu1.header.update('EXPOSURE',82.7273,'[d] time on source')
        hdu1.header.update('DEADC',0.909091,'deadtime correction')
        hdu1.header.update('TIMEPIXR',0.5,'bin time beginning=0 middle=0.5 end=1')
        hdu1.header.update('TIERRELA',5.78e-7,'[d] relative time error')
        hdu1.header.update('TIERABSO',5.78e-6,'[d] absolute time error')
        hdu1.header.update('INT_TIME',6.0198,'[s] photon accumulation time per frame')
        hdu1.header.update('READTIME',0.518948526144,'[s] readout time per frame')
        hdu1.header.update('FRAMETIM',6.53875,'[s] frame time (INT_TIME + READTIME)')
        hdu1.header.update('NUM_FRM',270,'number of frames per time stamp')
        hdu1.header.update('TIMEDEL',30.0/1440,'[d] time resolution of data')
        hdu1.header.update('DATE-OBS','2014-09-30T12:28:48.0','TSTART as UT calendar date')
        hdu1.header.update('DATE-END','2014-10-02T11:02:24.0','TSTOP as UT calendar date')
        hdu1.header.update('BACKAPP',True,'background is subtracted')
        hdu1.header.update('DEADAPP',False,'deadtime applied')
        hdu1.header.update('VIGNAPP',False,'vignetting or collimator correction applied')
        hdu1.header.update('GAIN',104.88,'[electrons/count] channel gain')
        hdu1.header.update('READNOIS',0.7845*104.88,'[electrons] read noise')
        hdu1.header.update('NREADOUT',276,'number of reads per cadence')
        hdu1.header.update('TIMSLICE',3,'time-slice readout sequence section')
        hdu1.header.update('MEANBLCK',738,'[count] FSW mean black level')
        hdulist.append(hdu1)
        
# create the outfile image extension

    if status == 0:
        hdu2 = ImageHDU(aperture)
        hdu2.header.update('INHERIT',True,'inherit primary keywords')
        hdu2.header.update('EXTNAME','APERTURE','extension name')
        hdu2.header.update('EXTVER',1,'extension version number')
        hdu2.header.update('TELESCOP','Kepler','telescope')
        hdu2.header.update('INSTRUME','Kepler photometer','detector type')
        hdu2.header.update('OBJECT','KIC 12345678','string version of kepID')
        hdu2.header.update('KEPLERID',1234567,'unique Kepler target identifier')
        hdu2.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
        hdu2.header.update('RA_OBJ',294.94017,'[deg] right ascension from KIC')
        hdu2.header.update('DEC_OBJ',43.80033,'[deg] declination from KIC')
        hdu2.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
        hdu2.header.update('WCSAXES',2,'number of WCS axes')
        hdu2.header.update('CTYPE1P','RAWX','pixel coordinate type')
        hdu2.header.update('CTYPE2P','RAWY','pixel coordinate type')
        hdu2.header.update('CRPIX1P',1,'[pixel] reference pixel along image axis 1')
        hdu2.header.update('CRPIX2P',1,'[pixel] reference pixel along image axis 2')
        hdu2.header.update('CRVAL1P',DATx[0],'[pixel] column number at reference pixel')
        hdu2.header.update('CRVAL2P',DATy[0],'[pixel] row number at reference pixel')
        hdu2.header.update('CUNIT1P','pixel','physical unit in column dimension')
        hdu2.header.update('CUNIT2P','pixel','physical unit in row dimension')
        hdu2.header.update('CDELT1P',1.0,'[pixel] pixel scale along columns')
        hdu2.header.update('CDELT2P',1.0,'[pixel] pixel scale along rows')
        hdu2.header.update('CTYPE1','RA---TAN','right ascension coordinate type')
        hdu2.header.update('CTYPE2','DEC--TAN','declination coordinate type')
        hdu2.header.update('CRPIX1',int(xdim/2),'[pixel] reference pixel along image axis 1')
        hdu2.header.update('CRPIX2',int(ydim/2),'[pixel] reference pixel along image axis 2')
        hdu2.header.update('CRVAL1',294.94017,'right ascension at reference pixel [deg]')
        hdu2.header.update('CRVAL2',43.80033,'declination at reference pixel [deg]')
        hdu2.header.update('CUNIT1','deg','physical unit in column dimension')
        hdu2.header.update('CUNIT2','deg','physical unit in row dimension')
        hdu2.header.update('CDELT1',-0.001106815552144,'pixel scale in RA dimension')
        hdu2.header.update('CDELT2',0.001106815552144,'pixel scale in Dec dimension')
        hdu2.header.update('PC1_1',0.46086006096337634,'linear transformation matrix element cos(th)')
        hdu2.header.update('PC1_2',-0.8897046441865888,'linear transformation matrix element -sin(th)')
        hdu2.header.update('PC2_1',0.8864170184391076,'linear transformation matrix element sin(th)')
        hdu2.header.update('PC2_2',0.45860051653617395,'linear transformation matrix element cos(th)')
        hdulist.append(hdu2)

# write output file
        
    if status == 0:
        hdulist.writeto(outfile,checksum=True)

# stop time

    kepmsg.clock('\nKEPFAKE ended at',logfile,verbose)

    return

# -----------------------------------------------------------
# plot channel image

def plotimage(imgflux_pl,zminfl,zmaxfl,plmode,row,column,xdim,ydim,tlabel,colmap):

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

    ax = pylab.axes([0.1,0.1,0.8,0.8])
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
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Fitting PRF model to Target Pixel image')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('outfile', help='Name of output target pixel file', type=str)
    parser.add_argument('--modile', '-b', help='Order of background polynmial fit', default=1, dest='border', type=int)


    parser.add_argument('--rownum', '-r', default=2200, help='Row number of image stored in infile', dest='rownum', type=int)
    parser.add_argument('columns', help='Column number of each source to be fit', type=str)
    parser.add_argument('rows', help='Row number of each source to be fit', type=str)
    parser.add_argument('fluxes', help='Relative flux of each source to be fit', type=str)
    parser.add_argument('--border', '-b', help='Order of background polynmial fit', default=1, dest='border', type=int)
    parser.add_argument('--background', action='store_true', help='Fit background?', default=False)
    parser.add_argument('--focus', action='store_true', help='Fit focus changes?', default=False)
    parser.add_argument('prfdir', help='Folder containing Point Response Function FITS files', type=str)
    parser.add_argument('--xtol', '-x', default=1.0e-4, help='Fit parameter tolerance', dest='xtol', type=float)
    parser.add_argument('--ftol', '-f', default=1.0, help='Fit minimization tolerance', dest='ftol', type=float)
    parser.add_argument('--imscale', '-i', help='Type of image intensity scale', default='linear', dest='imscale', type=str,choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', '-c', help='Image colormap', default='YlOrBr', dest='cmap', type=str,choices=['Accent','Blues','BrBG','BuGn','BuPu','Dark2','GnBu','Greens','Greys','OrRd','Oranges','PRGn','Paired','Pastel1','Pastel2','PiYG','PuBu','PuBuGn','PuOr','PuRd','Purples','RdBu','RdGy','RdPu','RdYlBu','RdYlGn','Reds','Set1','Set2','Set3','Spectral','YlGn','YlGnBu','YlOrBr','YlOrRd','afmhot','autumn','binary','bone','brg','bwr','cool','copper','flag','gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_yarg','gnuplot','gnuplot2','gray','hot','hsv','jet','ocean','pink','prism','rainbow','seismic','spectral','spring','summer','terrain','winter','browse'])
    parser.add_argument('--plot', action='store_true', help='Plot fit results?', default=False)
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepfakephot.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepfake(args.infile,args.plotfile,args.rownum,args.columns,args.rows,args.fluxes,args.border,
           args.background,args.focus,args.prfdir,args.xtol,args.ftol,args.imscale,args.cmap,
           args.plot,args.verbose,args.logfile,args.status,cmdLine)
    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepfake.par")
    t = iraf.IrafTaskFactory(taskname="kepfake", value=parfile, function=kepfake)
