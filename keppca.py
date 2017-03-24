import sys, re
import numpy, scipy, pylab
from astropy.io import fits as pyfits
from numpy import *
from pylab import *
from matplotlib import *
import kepmsg, kepio, kepkey, kepplot
import random

# -----------------------------------------------------------
# core code

def keppca(infile,maskfile,outfile,components,plotpca,nreps,clobber,verbose,logfile,status,cmdLine=False): 

    try:
        import mdp
    except:
        msg = 'ERROR -- KEPPCA: this task has an external python dependency to MDP, a Modular toolkit for Data Processing (http://mdp-toolkit.sourceforge.net). In order to take advantage of this PCA task, the user must first install MDP with their current python distribution. Note carefully that you may have more than python installation on your machine, and ensure that MDP is installed with the same version of python that the PyKE tools employ. Installation instructions for MDP can be found at the URL provided above.'
        status = kepmsg.err(None,msg,True)
    
# startup parameters

    status = 0
    labelsize = 32
    ticksize = 18
    xsize = 16
    ysize = 10
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2
    seterr(all="ignore") 

# log the call 

    if status == 0:
        hashline = '----------------------------------------------------------------------------'
        kepmsg.log(logfile,hashline,verbose)
        call = 'KEPPCA -- '
        call += 'infile='+infile+' '
        call += 'maskfile='+maskfile+' '
        call += 'outfile='+outfile+' '
        call += 'components='+components+' '
        ppca = 'n'
        if (plotpca): ppca = 'y'
        call += 'plotpca='+ppca+ ' '
        call += 'nmaps='+str(nreps)+' '
        overwrite = 'n'
        if (clobber): overwrite = 'y'
        call += 'clobber='+overwrite+ ' '
        chatter = 'n'
        if (verbose): chatter = 'y'
        call += 'verbose='+chatter+' '
        call += 'logfile='+logfile
        kepmsg.log(logfile,call+'\n',verbose)
        
# start time

    if status == 0:
        kepmsg.clock('KEPPCA started at',logfile,verbose)

# test log file

    if status == 0:
        logfile = kepmsg.test(logfile)
    
# clobber output file

    if status == 0:
        if clobber: status = kepio.clobber(outfile,logfile,verbose)
        if kepio.fileexists(outfile): 
            message = 'ERROR -- KEPPCA: ' + outfile + ' exists. Use clobber=yes'
            status = kepmsg.err(logfile,message,verbose)

# Set output file names - text file with data and plot

    if status == 0:
        dataout = copy(outfile)
        repname = re.sub('.fits','.png',outfile)

# open input file

    if status == 0:    
        instr = pyfits.open(infile,mode='readonly',memmap=True)
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# open TPF FITS file

    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
            kepio.readTPF(infile,'TIME',logfile,verbose)
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
            ra, dec, column, row, kepmag, xdim, ydim, flux_bkg, status = \
            kepio.readTPF(infile,'FLUX_BKG',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, flux_bkg_err, status = \
            kepio.readTPF(infile,'FLUX_BKG_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, pcorr1, status = \
            kepio.readTPF(infile,'POS_CORR1',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, pcorr2, status = \
            kepio.readTPF(infile,'POS_CORR2',logfile,verbose)

# Save original data dimensions, in case of using maskfile

    if status == 0:
        xdimorig = xdim
        ydimorig = ydim
    
# read mask definition file if it has been supplied

    if status == 0 and 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        maskx = array([],'int')
        masky = array([],'int')
        lines, status = kepio.openascii(maskfile,'r',logfile,verbose)
        for line in lines:
            line = line.strip().split('|')
            if len(line) == 6:
                y0 = int(line[3])
                x0 = int(line[4])
                line = line[5].split(';')
                for items in line:
                    try:
                        masky = numpy.append(masky,y0 + int(items.split(',')[0]))
                        maskx = numpy.append(maskx,x0 + int(items.split(',')[1]))
                    except:
                        continue
        status = kepio.closeascii(lines,logfile,verbose)
        if len(maskx) == 0 or len(masky) == 0:
            message = 'ERROR -- KEPPCA: ' + maskfile + ' contains no pixels.'
            status = kepmsg.err(logfile,message,verbose)
        xdim = max(maskx) - min(maskx) + 1   # Find largest x dimension of mask
        ydim = max(masky) - min(masky) + 1   # Find largest y dimension of mask

# pad mask to ensure it is rectangular

        workx = array([],'int')
        worky = array([],'int')
        for ip in arange(min(maskx),max(maskx) + 1):
            for jp in arange(min(masky),max(masky) + 1):
                workx = append(workx,ip)
                worky = append(worky,jp)
        maskx = workx
        masky = worky

# define new subimage bitmap...

    if status == 0 and maskfile.lower() != 'all':
        aperx = numpy.array([],'int')
        apery = numpy.array([],'int')
        aperb = maskx - x0 + xdimorig * (masky - y0)   # aperb is an array that contains the pixel numbers in the mask
        npix = len(aperb)

# ...or use all pixels

    if status == 0 and maskfile.lower() == 'all':
        npix = xdimorig*ydimorig
        aperb = array([],'int')
        aperb = numpy.r_[0:npix]

# legal mask defined?

    if status == 0:
        if len(aperb) == 0:
            message = 'ERROR -- KEPPCA: no legal pixels within the subimage are defined.'
            status = kepmsg.err(logfile,message,verbose)

# Identify principal components desired

    if status == 0:
        pcaout = []
        txt = components.strip().split(',')
        for work1 in txt:
            try:
                pcaout.append(int(work1.strip()))
            except:
                work2 = work1.strip().split('-')
                try:
                    for work3 in range(int(work2[0]),int(work2[1]) + 1):
                        pcaout.append(work3)
                except:
                    message = 'ERROR -- KEPPCA: cannot understand principal component list requested'
                    status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        pcaout = set(sort(pcaout))
    pcarem = array(list(pcaout))-1    # The list of pca component numbers to be removed

# Initialize arrays and variables, and apply pixel mask to the data

    if status == 0:
        ntim = 0
        time = numpy.array([],dtype='float64')
        timecorr = numpy.array([],dtype='float32')
        cadenceno = numpy.array([],dtype='int')
        pixseries = numpy.array([],dtype='float32')
        errseries = numpy.array([],dtype='float32')
        bkgseries = numpy.array([],dtype='float32')
        berseries = numpy.array([],dtype='float32')
        quality = numpy.array([],dtype='float32')
        pos_corr1 = numpy.array([],dtype='float32')
        pos_corr2 = numpy.array([],dtype='float32')
        nrows = numpy.size(fluxpixels,0)
        
# Apply the pixel mask so we are left with only the desired pixels       

    if status == 0:
        pixseriesb = fluxpixels[:,aperb]
        errseriesb = errpixels[:,aperb]
        bkgseriesb = flux_bkg[:,aperb]
        berseriesb = flux_bkg_err[:,aperb]

# Read in the data to various arrays 
   
    if status == 0:
        for i in range(nrows):
            if qual[i] < 10000 and \
                    numpy.isfinite(barytime[i]) and \
                    numpy.isfinite(fluxpixels[i,int(ydim*xdim/2+0.5)]) and \
                    numpy.isfinite(fluxpixels[i,1+int(ydim*xdim/2+0.5)]):
                ntim += 1
                time = numpy.append(time,barytime[i])
                timecorr = numpy.append(timecorr,tcorr[i])
                cadenceno = numpy.append(cadenceno,cadno[i])
                pixseries = numpy.append(pixseries,pixseriesb[i])
                errseries = numpy.append(errseries,errseriesb[i])
                bkgseries = numpy.append(bkgseries,bkgseriesb[i])
                berseries = numpy.append(berseries,berseriesb[i])
                quality = numpy.append(quality,qual[i])
                pos_corr1 = numpy.append(pos_corr1,pcorr1[i])
                pos_corr2 = numpy.append(pos_corr2,pcorr2[i])
        pixseries = numpy.reshape(pixseries,(ntim,npix))
        errseries = numpy.reshape(errseries,(ntim,npix))
        bkgseries = numpy.reshape(bkgseries,(ntim,npix))
        berseries = numpy.reshape(berseries,(ntim,npix))        
        tmp =  numpy.median(pixseries,axis=1)     
        for i in range(len(tmp)):
             pixseries[i] = pixseries[i] - tmp[i]

# Figure out which pixels are undefined/nan and remove them. Keep track for adding back in later

    if status == 0:
        nanpixels = numpy.array([],dtype='int')
        i = 0
        while (i < npix):
            if numpy.isnan(pixseries[0,i]):
                nanpixels = numpy.append(nanpixels,i)
                npix = npix - 1
            i = i + 1
        pixseries = numpy.delete(pixseries,nanpixels,1)
        errseries = numpy.delete(errseries,nanpixels,1)
        pixseries[numpy.isnan(pixseries)] = random.gauss(100,10)
        errseries[numpy.isnan(errseries)] = 10
 
# Compute statistical weights, means, standard deviations

    if status == 0:
        weightseries = (pixseries/errseries)**2
        pixMean = numpy.average(pixseries,axis=0,weights=weightseries)
        pixStd  = numpy.std(pixseries,axis=0)

# Normalize the input by subtracting the mean and divising by the standard deviation. 
# This makes it a correlation-based PCA, which is what we want.

    if status == 0:
        pixseriesnorm = (pixseries - pixMean)/pixStd

# Number of principal components to compute. Setting it equal to the number of pixels

    if status == 0:
        nvecin = npix  

# Run PCA using the MDP Whitening PCA, which produces normalized PCA components (zero mean and unit variance)
    
    if status == 0:
        pcan = mdp.nodes.WhiteningNode(svd=True)
        pcar = pcan.execute(pixseriesnorm)
        eigvec = pcan.get_recmatrix()
        model = pcar
 
# Re-insert nan columns as zeros

    if status == 0:
        for i in range(0,len(nanpixels)):
            nanpixels[i] = nanpixels[i]-i
        eigvec = numpy.insert(eigvec,nanpixels,0,1)
        pixMean = numpy.insert(pixMean,nanpixels,0,0)

#  Make output eigenvectors (correlation images) into xpix by ypix images

    if status == 0:
        eigvec = eigvec.reshape(nvecin,ydim,xdim)

# Calculate sum of all pixels to display as raw lightcurve and other quantities

    if status == 0:
        pixseriessum = sum(pixseries,axis=1)
        nrem=len(pcarem)  # Number of components to remove
        nplot = npix      # Number of pcas to plot - currently set to plot all components, but could set 
                          # nplot = nrem to just plot as many components as is being removed

# Subtract components by fitting them to the summed light curve

    if status == 0:
        x0 = numpy.tile(-1.0,1)
        for k in range(0,nrem):
            def f(x):
                fluxcor = pixseriessum
                for k in range(0,len(x)):
                    fluxcor = fluxcor - x[k]*model[:,pcarem[k]]
                return mad(fluxcor)
            if k==0:
                x0 = array([-1.0])
            else:
                x0 = numpy.append(x0,1.0)
            myfit = scipy.optimize.fmin(f,x0,maxiter=50000,maxfun=50000,disp=False)
            x0 = myfit
    
# Now that coefficients for all components have been found, subtract them to produce a calibrated time-series, 
# and then divide by the robust mean to produce a normalized time series as well

    if status == 0:
        c = myfit
        fluxcor = pixseriessum
        for k in range(0,nrem):
            fluxcor = fluxcor - c[k]*model[:,pcarem[k]]
            normfluxcor = fluxcor/mean(reject_outliers(fluxcor,2))

# input file data

    if status == 0:
        cards0 = instr[0].header.cards
        cards1 = instr[1].header.cards
        cards2 = instr[2].header.cards
        table = instr[1].data[:]
        maskmap = copy(instr[2].data)

# subimage physical WCS data

    if status == 0:
        crpix1p = cards2['CRPIX1P'].value
        crpix2p = cards2['CRPIX2P'].value
        crval1p = cards2['CRVAL1P'].value
        crval2p = cards2['CRVAL2P'].value
        cdelt1p = cards2['CDELT1P'].value
        cdelt2p = cards2['CDELT2P'].value

# dummy columns for output file

    if status == 0:
        sap_flux_err = numpy.empty(len(time)); sap_flux_err[:] = numpy.nan
        sap_bkg = numpy.empty(len(time)); sap_bkg[:] = numpy.nan
        sap_bkg_err = numpy.empty(len(time)); sap_bkg_err[:] = numpy.nan
        pdc_flux = numpy.empty(len(time)); pdc_flux[:] = numpy.nan
        pdc_flux_err = numpy.empty(len(time)); pdc_flux_err[:] = numpy.nan
        psf_centr1 = numpy.empty(len(time)); psf_centr1[:] = numpy.nan
        psf_centr1_err = numpy.empty(len(time)); psf_centr1_err[:] = numpy.nan
        psf_centr2 = numpy.empty(len(time)); psf_centr2[:] = numpy.nan
        psf_centr2_err = numpy.empty(len(time)); psf_centr2_err[:] = numpy.nan
        mom_centr1 = numpy.empty(len(time)); mom_centr1[:] = numpy.nan
        mom_centr1_err = numpy.empty(len(time)); mom_centr1_err[:] = numpy.nan
        mom_centr2 = numpy.empty(len(time)); mom_centr2[:] = numpy.nan
        mom_centr2_err = numpy.empty(len(time)); mom_centr2_err[:] = numpy.nan

# mask bitmap

    if status == 0 and 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i,j] == 0:
                    pass
                else:
                    maskmap[i,j] = 1
                    for k in range(len(maskx)):
                        if aperx[-1] == maskx[k] and apery[-1] == masky[k]:
                            maskmap[i,j] = 3

# construct output primary extension

    if status == 0:
        hdu0 = pyfits.PrimaryHDU()
        for i in range(len(cards0)):
            if cards0[i].keyword not in hdu0.header.keys():
                hdu0.header[cards0[i].keyword] = (cards0[i].value, cards0[i].comment)
            else:
                hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
        status = kepkey.history(call,hdu0,outfile,logfile,verbose)
        outstr = HDUList(hdu0)

# construct output light curve extension

    if status == 0:
        col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
        col2 = Column(name='TIMECORR',format='E',unit='d',array=timecorr)
        col3 = Column(name='CADENCENO',format='J',array=cadenceno)
        col4 = Column(name='SAP_FLUX',format='E',unit='e-/s',array=pixseriessum)
        col5 = Column(name='SAP_FLUX_ERR',format='E',unit='e-/s',array=sap_flux_err)
        col6 = Column(name='SAP_BKG',format='E',unit='e-/s',array=sap_bkg)
        col7 = Column(name='SAP_BKG_ERR',format='E',unit='e-/s',array=sap_bkg_err)
        col8 = Column(name='PDCSAP_FLUX',format='E',unit='e-/s',array=pdc_flux)
        col9 = Column(name='PDCSAP_FLUX_ERR',format='E',unit='e-/s',array=pdc_flux_err)
        col10 = Column(name='SAP_QUALITY',format='J',array=quality)
        col11 = Column(name='PSF_CENTR1',format='E',unit='pixel',array=psf_centr1)
        col12 = Column(name='PSF_CENTR1_ERR',format='E',unit='pixel',array=psf_centr1_err)
        col13 = Column(name='PSF_CENTR2',format='E',unit='pixel',array=psf_centr2)
        col14 = Column(name='PSF_CENTR2_ERR',format='E',unit='pixel',array=psf_centr2_err)
        col15 = Column(name='MOM_CENTR1',format='E',unit='pixel',array=mom_centr1)
        col16 = Column(name='MOM_CENTR1_ERR',format='E',unit='pixel',array=mom_centr1_err)
        col17 = Column(name='MOM_CENTR2',format='E',unit='pixel',array=mom_centr2)
        col18 = Column(name='MOM_CENTR2_ERR',format='E',unit='pixel',array=mom_centr2_err)
        col19 = Column(name='POS_CORR1',format='E',unit='pixel',array=pos_corr1)
        col20 = Column(name='POS_CORR2',format='E',unit='pixel',array=pos_corr2)
        col21 = Column(name='PCA_FLUX',format='E',unit='e-/s',array=fluxcor)
        col22 = Column(name='PCA_FLUX_NRM',format='E',array=normfluxcor)
        cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11, \
                            col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22])
        hdu1 = new_table(cols)
        hdu1.header['TTYPE1'] = ('TIME','column title: data time stamps')
        hdu1.header['TFORM1'] = ('D','data type: float64')
        hdu1.header['TUNIT1'] = ('BJD - 2454833','column units: barycenter corrected JD')
        hdu1.header['TDISP1'] = ('D12.7','column display format')
        hdu1.header['TTYPE2'] = ('TIMECORR','column title: barycentric-timeslice correction')
        hdu1.header['TFORM2'] = ('E','data type: float32')
        hdu1.header['TUNIT2'] = ('d','column units: days')
        hdu1.header['TTYPE3'] = ('CADENCENO','column title: unique cadence number')
        hdu1.header['TFORM3'] = ('J','column format: signed integer32')
        hdu1.header['TTYPE4'] = ('SAP_FLUX','column title: aperture photometry flux')
        hdu1.header['TFORM4'] = ('E','column format: float32')
        hdu1.header['TUNIT4'] = ('e-/s','column units: electrons per second')
        hdu1.header['TTYPE5'] = ('SAP_FLUX_ERR','column title: aperture phot. flux error')
        hdu1.header['TFORM5'] = ('E','column format: float32')
        hdu1.header['TUNIT5'] = ('e-/s','column units: electrons per second (1-sigma)')
        hdu1.header['TTYPE6'] = ('SAP_BKG','column title: aperture phot. background flux')
        hdu1.header['TFORM6'] = ('E','column format: float32')
        hdu1.header['TUNIT6'] = ('e-/s','column units: electrons per second')
        hdu1.header['TTYPE7'] = ('SAP_BKG_ERR','column title: ap. phot. background flux error')
        hdu1.header['TFORM7'] = ('E','column format: float32')
        hdu1.header['TUNIT7'] = ('e-/s','column units: electrons per second (1-sigma)')
        hdu1.header['TTYPE8'] = ('PDCSAP_FLUX','column title: PDC photometry flux')
        hdu1.header['TFORM8'] = ('E','column format: float32')
        hdu1.header['TUNIT8'] = ('e-/s','column units: electrons per second')
        hdu1.header['TTYPE9'] = ('PDCSAP_FLUX_ERR','column title: PDC flux error')
        hdu1.header['TFORM9'] = ('E','column format: float32')
        hdu1.header['TUNIT9'] = ('e-/s','column units: electrons per second (1-sigma)')
        hdu1.header['TTYPE10'] = ('SAP_QUALITY','column title: aperture photometry quality flag')
        hdu1.header['TFORM10'] = ('J','column format: signed integer32')
        hdu1.header['TTYPE11'] = ('PSF_CENTR1','column title: PSF fitted column centroid')
        hdu1.header['TFORM11'] = ('E','column format: float32')
        hdu1.header['TUNIT11'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE12'] = ('PSF_CENTR1_ERR','column title: PSF fitted column error')
        hdu1.header['TFORM12'] = ('E','column format: float32')
        hdu1.header['TUNIT12'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE13'] = ('PSF_CENTR2','column title: PSF fitted row centroid')
        hdu1.header['TFORM13'] = ('E','column format: float32')
        hdu1.header['TUNIT13'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE14'] = ('PSF_CENTR2_ERR','column title: PSF fitted row error')
        hdu1.header['TFORM14'] = ('E','column format: float32')
        hdu1.header['TUNIT14'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE15'] = ('MOM_CENTR1','column title: moment-derived column centroid')
        hdu1.header['TFORM15'] = ('E','column format: float32')
        hdu1.header['TUNIT15'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE16'] = ('MOM_CENTR1_ERR','column title: moment-derived column error')
        hdu1.header['TFORM16'] = ('E','column format: float32')
        hdu1.header['TUNIT16'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE17'] = ('MOM_CENTR2','column title: moment-derived row centroid')
        hdu1.header['TFORM17'] = ('E','column format: float32')
        hdu1.header['TUNIT17'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE18'] = ('MOM_CENTR2_ERR','column title: moment-derived row error')
        hdu1.header['TFORM18'] = ('E','column format: float32')
        hdu1.header['TUNIT18'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE19'] = ('POS_CORR1','column title: col correction for vel. abbern')
        hdu1.header['TFORM19'] = ('E','column format: float32')
        hdu1.header['TUNIT19'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE20'] = ('POS_CORR2','column title: row correction for vel. abbern')
        hdu1.header['TFORM20'] = ('E','column format: float32')
        hdu1.header['TUNIT20'] = ('pixel','column units: pixel')
        hdu1.header['TTYPE21'] = ('PCA_FLUX','column title: PCA-corrected flux')
        hdu1.header['TFORM21'] = ('E','column format: float32')
        hdu1.header['TUNIT21'] = ('pixel','column units: e-/s')
        hdu1.header['TTYPE22'] = ('PCA_FLUX_NRM','column title: normalized PCA-corrected flux')
        hdu1.header['TFORM22'] = ('E','column format: float32')
        hdu1.header['EXTNAME'] = ('LIGHTCURVE','name of extension')
        for i in range(len(cards1)):
            if (cards1[i].keyword not in hdu1.header.keys() and
                cards1[i].keyword[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC']):
                hdu1.header[cards1[i].keyword] = (cards1[i].value, cards1[i].comment)
        outstr.append(hdu1)

# construct output mask bitmap extension

    if status == 0:
        hdu2 = ImageHDU(maskmap)
        for i in range(len(cards2)):
            if cards2[i].keyword not in hdu2.header.keys():
                hdu2.header[cards2[i].keyword] = (cards2[i].value, cards2[i].comment)
            else:
                hdu2.header.cards[cards2[i].keyword].comment = cards2[i].comment
        outstr.append(hdu2)

# construct principal component table

    if status == 0:
        cols = [Column(name='TIME',format='E',unit='BJD - 2454833',array=time)]
        for i in range(len(pcar[0,:])):
            colname = 'PC' + str(i + 1)
            col = Column(name=colname,format='E',array=pcar[:,i])
            cols.append(col)
        hdu3 = new_table(ColDefs(cols))
        hdu3.header['EXTNAME'] = ('PRINCIPAL_COMPONENTS','name of extension')
        hdu3.header['TTYPE1'] = ('TIME','column title: data time stamps')
        hdu3.header['TFORM1'] = ('D','data type: float64')
        hdu3.header['TUNIT1'] = ('BJD - 2454833','column units: barycenter corrected JD')
        hdu3.header['TDISP1'] = ('D12.7','column display format')
        for i in range(len(pcar[0,:])):
            hdu3.header['TTYPE' + str(i + 2)] = \
                ('PC' + str(i + 1), 'column title: principal component number' + str(i + 1))
            hdu3.header['TFORM' + str(i + 2)] = ('E','column format: float32')
        outstr.append(hdu3)

# write output file

    if status == 0:
        outstr.writeto(outfile)
    
# close input structure

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)
        
# Create PCA report 

    if status == 0 and plotpca:
        npp = 7 # Number of plots per page
        l = 1
        repcnt = 1
        for k in range(nreps):

# First plot of every pagewith flux image, flux and calibrated time series 

            status = kepplot.define(16,12,logfile,verbose)
            if (k % (npp - 1) == 0):     
                pylab.figure(figsize=[10,16])
                subplot2grid((npp,6),(0,0), colspan=2)
#                imshow(log10(pixMean.reshape(xdim,ydim).T-min(pixMean)+1),interpolation="nearest",cmap='RdYlBu')
                imshow(log10(flipud(pixMean.reshape(ydim,xdim))-min(pixMean)+1),interpolation="nearest",cmap='RdYlBu')
                xticks([])
                yticks([])
                ax1 = subplot2grid((npp,6),(0,2), colspan=4)
                px = copy(time) + bjdref
                py = copy(pixseriessum)
                px, xlab, status = kepplot.cleanx(px,logfile,verbose) 
                py, ylab, status = kepplot.cleany(py,1.0,logfile,verbose)
                kepplot.RangeOfPlot(px,py,0.01,False)
                kepplot.plot1d(px,py,cadence,lcolor,lwidth,fcolor,falpha,True)
                py = copy(fluxcor)
                py, ylab, status = kepplot.cleany(py,1.0,logfile,verbose)
                plot(px,py,marker='.',color='r',linestyle='',markersize=1.0)
                kepplot.labels('',re.sub('\)','',re.sub('Flux \(','',ylab)),'k',18)
                grid()
                setp(ax1.get_xticklabels(), visible=False)

# plot principal components

            subplot2grid((npp,6),(l,0), colspan=2)
            imshow(eigvec[k],interpolation="nearest",cmap='RdYlBu')
            xlim(-0.5,xdim-0.5)
            ylim(-0.5,ydim-0.5)
            xticks([])
            yticks([])

# The last plot on the page that should have the xlabel

            if ( k% (npp - 1) == npp - 2 or k == nvecin - 1):  
                subplot2grid((npp,6),(l,2), colspan=4)
                py = copy(model[:,k])
                kepplot.RangeOfPlot(px,py,0.01,False)
                kepplot.plot1d(px,py,cadence,'r',lwidth,'g',falpha,True)
                kepplot.labels(xlab,'PC ' + str(k+1),'k',18)
                pylab.grid()
                pylab.tight_layout()
                l = 1
                pylab.savefig(re.sub('.png','_%d.png' % repcnt,repname))
                if not cmdLine: kepplot.render(cmdLine)
                repcnt += 1

# The other plots on the page that should have no xlabel

            else:
                ax2 = subplot2grid((npp,6),(l,2), colspan=4)
                py = copy(model[:,k])
                kepplot.RangeOfPlot(px,py,0.01,False)
                kepplot.plot1d(px,py,cadence,'r',lwidth,'g',falpha,True)
                kepplot.labels('','PC ' + str(k+1),'k',18)
                grid()
                setp(ax2.get_xticklabels(), visible=False)
                pylab.tight_layout()
                l=l+1
        pylab.savefig(re.sub('.png','_%d.png' % repcnt,repname))
        if not cmdLine: kepplot.render(cmdLine)

# plot style and size

    if status == 0 and plotpca:
        status = kepplot.define(labelsize,ticksize,logfile,verbose)
        pylab.figure(figsize=[xsize,ysize])
        pylab.clf()

# plot aperture photometry and PCA corrected data

    if status == 0 and plotpca:
        ax = kepplot.location([0.06,0.54,0.93,0.43])
        px = copy(time) + bjdref
        py = copy(pixseriessum)
        px, xlab, status = kepplot.cleanx(px,logfile,verbose) 
        py, ylab, status = kepplot.cleany(py,1.0,logfile,verbose)
        kepplot.RangeOfPlot(px,py,0.01,False)
        kepplot.plot1d(px,py,cadence,lcolor,lwidth,fcolor,falpha,True)
        py = copy(fluxcor)
        py, ylab, status = kepplot.cleany(py,1.0,logfile,verbose)
        kepplot.plot1d(px,py,cadence,'r',2,fcolor,0.0,True)
        pylab.setp(pylab.gca(),xticklabels=[])
        kepplot.labels('',ylab,'k',24)
        pylab.grid()

# plot aperture photometry and PCA corrected data

    if status == 0 and plotpca:
        ax = kepplot.location([0.06,0.09,0.93,0.43])
        yr = array([],'float32')
        npc = min([6,nrem])
        for i in range(npc-1,-1,-1):
            py = pcar[:,i] * c[i]
            py, ylab, status = kepplot.cleany(py,1.0,logfile,verbose)
            cl = float(i) / (float(npc))
            kepplot.plot1d(px,py,cadence,[1.0-cl,0.0,cl],2,fcolor,0.0,True)
            yr = append(yr,py)
        y1 = max(yr)
        y2 = -min(yr)
        kepplot.RangeOfPlot(px,array([-y1,y1,-y2,y2]),0.01,False)
        kepplot.labels(xlab,'Principal Components','k',24)
        pylab.grid()

# save plot to file

    if status == 0 and plotpca:
        pylab.savefig(repname)

# render plot

    if status == 0 and plotpca:
        kepplot.render(cmdLine)

# stop time

    if status == 0:
        kepmsg.clock('KEPPCA ended at',logfile,verbose)

    return


# -----------------------------------------------------------
# Outlier rejection for computing robust mean later

def reject_outliers(data, m):
    return data[abs(data - numpy.mean(data)) < m * numpy.std(data)]

# -----------------------------------------------------------
# Mean absolute deviation function used for fitting the PCA components to the data and subtracting them out 
# Could replace with whateve minimization function you want

def mad(data):
    return mean(absolute(data - mean(data)))

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Correct aperture photmetry using target motion')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input target pixel FITS file', type=str)
    parser.add_argument('maskfile', help='Name of mask defintion ASCII file', type=str)
    parser.add_argument('outfile', help='Name of output FITS file', type=str)
    parser.add_argument('--components', default='1-3', help='Principal components to be removed', type=str)
    parser.add_argument('--plotpca', action='store_true', help='Create PCA plots?')
    parser.add_argument('--nmaps', default=10, help='Number of principal components to include in report', type=int)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='keppca.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine = True
    keppca(args.infile,args.maskfile,args.outfile,args.components,args.plotpca,
           args.nmaps,args.clobber,args.verbose,args.logfile,args.status,cmdLine)    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$keppca.par")
    t = iraf.IrafTaskFactory(taskname="keppca", value=parfile, function=keppca)
