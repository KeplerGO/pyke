
import numpy, sys, time, pyfits, pylab, math
from numpy import *
from pyfits import *
from pylab import *
from matplotlib import *
from scipy.optimize import leastsq
import kepio, kepmsg, kepkey, kepstat, kepfunc

# global variables

def kepextract(infile,maskfile,outfile,subback,clobber,verbose,logfile,status): 

# startup parameters

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPEXTRACT -- '
    call += 'infile='+infile+' '
    call += 'maskfile='+maskfile+' '
    call += 'outfile='+outfile+' '
    backgr = 'n'
    if (subback): backgr = 'y'
    call += 'background='+backgr+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPEXTRACT started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPEXTRACT: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    status = 0
    instr = pyfits.open(infile,mode='readonly',memmap=True)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# input file data

    if status == 0:
        cards0 = instr[0].header.cards
        cards1 = instr[1].header.cards
        cards2 = instr[2].header.cards
        table = instr[1].data[:]
        maskmap = copy(instr[2].data)

# input table data

    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, time, status = \
            kepio.readTPF(infile,'TIME',logfile,verbose)
        time = numpy.array(time,dtype='float64')
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, timecorr, status = \
            kepio.readTPF(infile,'TIMECORR',logfile,verbose)
        timecorr = numpy.array(timecorr,dtype='float32')
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, cadenceno, status = \
            kepio.readTPF(infile,'CADENCENO',logfile,verbose)
        cadenceno = numpy.array(cadenceno,dtype='int')
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, raw_cnts, status = \
            kepio.readTPF(infile,'RAW_CNTS',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, flux, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, flux_err, status = \
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
            ra, dec, column, row, kepmag, xdim, ydim, cosmic_rays, status = \
            kepio.readTPF(infile,'COSMIC_RAYS',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, quality, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)
        quality = numpy.array(quality,dtype='int')
    if status == 0:
        try:
            pos_corr1 = numpy.array(table.field('POS_CORR1'),dtype='float64')  #  ---for FITS wave #2
        except:
            pos_corr1 = empty(len(time)); pos_corr1[:] = numpy.nan   # ---temporary before FITS wave #2
        try:
            pos_corr2 = numpy.array(table.field('POS_CORR2'),dtype='float64')  #  ---for FITS wave #2
        except:
            pos_corr2 = empty(len(time)); pos_corr2[:] = numpy.nan   # ---temporary before FITS wave #2

# dummy columns for output file

        psf_centr1 = empty(len(time)); psf_centr1[:] = numpy.nan
        psf_centr1_err = empty(len(time)); psf_centr1_err[:] = numpy.nan
        psf_centr2 = empty(len(time)); psf_centr2[:] = numpy.nan
        psf_centr2_err = empty(len(time)); psf_centr2_err[:] = numpy.nan
#        mom_centr1 = empty(len(time)); mom_centr1[:] = numpy.nan
        mom_centr1_err = empty(len(time)); mom_centr1_err[:] = numpy.nan
#        mom_centr2 = empty(len(time)); mom_centr2[:] = numpy.nan
        mom_centr2_err = empty(len(time)); mom_centr2_err[:] = numpy.nan

# read mask definition file

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
                        masky = append(masky,y0 + int(items.split(',')[0]))
                        maskx = append(maskx,x0 + int(items.split(',')[1]))
                    except:
                        continue
        status = kepio.closeascii(lines,logfile,verbose)
        if len(maskx) == 0 or len(masky) == 0:
            message = 'ERROR -- KEPEXTRACT: ' + maskfile + ' contains no pixels.'
            status = kepmsg.err(logfile,message,verbose)

# subimage physical WCS data

    if status == 0:
        crpix1p = cards2['CRPIX1P'].value
        crpix2p = cards2['CRPIX2P'].value
        crval1p = cards2['CRVAL1P'].value
        crval2p = cards2['CRVAL2P'].value
        cdelt1p = cards2['CDELT1P'].value
        cdelt2p = cards2['CDELT2P'].value

# define new subimage bitmap...

    if status == 0 and 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        aperx = array([],'int')
        apery = array([],'int')
        aperb = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i,j] == 0:
                    aperb = append(aperb,0)
                else:
                    aperb = append(aperb,1)
                    maskmap[i,j] = 1
                    for k in range(len(maskx)):
                        if aperx[-1] == maskx[k] and apery[-1] == masky[k]:
                            aperb[-1] = 3
                            maskmap[i,j] = 3

# trap case where no aperture needs to be defined but pixel positions are still required for centroiding

    if status == 0 and maskfile.lower() == 'all':
        aperx = array([],'int')
        apery = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)

# ...or use old subimage bitmap

    if status == 0 and 'aper' in maskfile.lower():
        aperx = array([],'int')
        apery = array([],'int')
        aperb = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperb = append(aperb,maskmap[i,j])
                aperx = append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)

# ...or use all pixels

    if status == 0 and maskfile.lower() == 'all':
        aperb = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                if maskmap[i,j] == 0:
                    aperb = append(aperb,0)
                else:
                    aperb = append(aperb,3)
                    maskmap[i,j] = 3

# subtract median pixel value for background?

    if status == 0:
        sky = array([],'float32')
        for i in range(len(time)):
            sky = append(sky,median(flux[i,:]))
        if not subback:
            sky[:] = 0.0

# legal mask defined?

    if status == 0:
        if len(aperb) == 0:
            message = 'ERROR -- KEPEXTRACT: no legal pixels within the subimage are defined.'
            status = kepmsg.err(logfile,message,verbose)
        
# construct new table flux data

    if status == 0:
        naper = (aperb == 3).sum()
        ntime = len(time)
        sap_flux = array([],'float32')
        sap_flux_err = array([],'float32')
        sap_bkg = array([],'float32')
        sap_bkg_err = array([],'float32')
        raw_flux = array([],'float32')
        for i in range(len(time)):
            work1 = array([],'float64')
            work2 = array([],'float64')
            work3 = array([],'float64')
            work4 = array([],'float64')
            work5 = array([],'float64')
            for j in range(len(aperb)):
                if (aperb[j] == 3):
                    work1 = append(work1,flux[i,j]-sky[i])
                    work2 = append(work2,flux_err[i,j])
                    work3 = append(work3,flux_bkg[i,j])
                    work4 = append(work4,flux_bkg_err[i,j])
                    work5 = append(work5,raw_cnts[i,j])
            sap_flux = append(sap_flux,kepstat.sum(work1))
            sap_flux_err = append(sap_flux_err,kepstat.sumerr(work2))
            sap_bkg = append(sap_bkg,kepstat.sum(work3))
            sap_bkg_err = append(sap_bkg_err,kepstat.sumerr(work4))
            raw_flux = append(raw_flux,kepstat.sum(work5))

# construct new table moment data

    if status == 0:
        mom_centr1 = zeros(shape=(ntime))
        mom_centr2 = zeros(shape=(ntime))
        mom_centr1_err = zeros(shape=(ntime))
        mom_centr2_err = zeros(shape=(ntime))
        for i in range(ntime):
            xf = zeros(shape=(naper))
            yf = zeros(shape=(naper))
            f = zeros(shape=(naper))
            xfe = zeros(shape=(naper))
            yfe = zeros(shape=(naper))
            fe = zeros(shape=(naper))
            k = -1
            for j in range(len(aperb)):
                if (aperb[j] == 3):
                    k += 1
                    xf[k] = aperx[j] * flux[i,j]
                    xfe[k] = aperx[j] * flux_err[i,j]
                    yf[k] = apery[j] * flux[i,j]
                    yfe[k] = apery[j] * flux_err[i,j]
                    f[k] = flux[i,j]
                    fe[k] = flux_err[i,j]
            xfsum = kepstat.sum(xf)
            yfsum = kepstat.sum(yf)
            fsum = kepstat.sum(f)
            xfsume = sqrt(kepstat.sum(square(xfe)) / naper)
            yfsume = sqrt(kepstat.sum(square(yfe)) / naper)
            fsume = sqrt(kepstat.sum(square(fe)) / naper)
            mom_centr1[i] = xfsum / fsum
            mom_centr2[i] = yfsum / fsum
            mom_centr1_err[i] = sqrt((xfsume / xfsum)**2 + ((fsume / fsum)**2))
            mom_centr2_err[i] = sqrt((yfsume / yfsum)**2 + ((fsume / fsum)**2))
        mom_centr1_err = mom_centr1_err * mom_centr1
        mom_centr2_err = mom_centr2_err * mom_centr2

# construct new table PSF data

    if status == 0:
        psf_centr1 = zeros(shape=(ntime))
        psf_centr2 = zeros(shape=(ntime))
        psf_centr1_err = zeros(shape=(ntime))
        psf_centr2_err = zeros(shape=(ntime))
        modx = zeros(shape=(naper))
        mody = zeros(shape=(naper))
        k = -1
        for j in range(len(aperb)):
            if (aperb[j] == 3):
                k += 1
                modx[k] = aperx[j]
                mody[k] = apery[j]
        for i in range(ntime):
            modf = zeros(shape=(naper))
            k = -1
            guess = [mom_centr1[i], mom_centr2[i], nanmax(flux[i:]), 1.0, 1.0, 0.0, 0.0]
            for j in range(len(aperb)):
                if (aperb[j] == 3):
                    k += 1
                    modf[k] = flux[i,j]
                    args = (modx, mody, modf)
            try:
                ans = leastsq(kepfunc.PRFgauss2d,guess,args=args,xtol=1.0e-8,ftol=1.0e-4,full_output=True)
                s_sq = (ans[2]['fvec']**2).sum() / (ntime-len(guess))
                psf_centr1[i] = ans[0][0]
                psf_centr2[i] = ans[0][1]
            except:
                pass
            try:
                psf_centr1_err[i] = sqrt(diag(ans[1] * s_sq))[0]
            except:
                psf_centr1_err[i] = numpy.nan
            try:
                psf_centr2_err[i] = sqrt(diag(ans[1] * s_sq))[1]
            except:
                psf_centr2_err[i] = numpy.nan

# construct output primary extension

    if status == 0:
        hdu0 = pyfits.PrimaryHDU()
        for i in range(len(cards0)):
            if cards0[i].key not in hdu0.header.keys():
                hdu0.header.update(cards0[i].key, cards0[i].value, cards0[i].comment)
            else:
                hdu0.header.cards[cards0[i].key].comment = cards0[i].comment
        status = kepkey.history(call,hdu0,outfile,logfile,verbose)
        outstr = HDUList(hdu0)

# construct output light curve extension

    if status == 0:
        col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
        col2 = Column(name='TIMECORR',format='E',unit='d',array=timecorr)
        col3 = Column(name='CADENCENO',format='J',array=cadenceno)
        col4 = Column(name='SAP_FLUX',format='E',array=sap_flux)
        col5 = Column(name='SAP_FLUX_ERR',format='E',array=sap_flux_err)
        col6 = Column(name='SAP_BKG',format='E',array=sap_bkg)
        col7 = Column(name='SAP_BKG_ERR',format='E',array=sap_bkg_err)
        col8 = Column(name='PDCSAP_FLUX',format='E',array=sap_flux)
        col9 = Column(name='PDCSAP_FLUX_ERR',format='E',array=sap_flux_err)
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
        col21 = Column(name='RAW_FLUX',format='E',array=raw_flux)
        cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11, \
                            col12,col13,col14,col15,col16,col17,col18,col19,col20,col21])
        hdu1 = new_table(cols)
        hdu1.header.update('TTYPE1','TIME','column title: data time stamps')
        hdu1.header.update('TFORM1','D','data type: float64')
        hdu1.header.update('TUNIT1','BJD - 2454833','column units: barycenter corrected JD')
        hdu1.header.update('TDISP1','D12.7','column display format')
        hdu1.header.update('TTYPE2','TIMECORR','column title: barycentric-timeslice correction')
        hdu1.header.update('TFORM2','E','data type: float32')
        hdu1.header.update('TUNIT2','d','column units: days')
        hdu1.header.update('TTYPE3','CADENCENO','column title: unique cadence number')
        hdu1.header.update('TFORM3','J','column format: signed integer32')
        hdu1.header.update('TTYPE4','SAP_FLUX','column title: aperture photometry flux')
        hdu1.header.update('TFORM4','E','column format: float32')
        hdu1.header.update('TUNIT4','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE5','SAP_FLUX_ERR','column title: aperture phot. flux error')
        hdu1.header.update('TFORM5','E','column format: float32')
        hdu1.header.update('TUNIT5','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE6','SAP_BKG','column title: aperture phot. background flux')
        hdu1.header.update('TFORM6','E','column format: float32')
        hdu1.header.update('TUNIT6','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE7','SAP_BKG_ERR','column title: ap. phot. background flux error')
        hdu1.header.update('TFORM7','E','column format: float32')
        hdu1.header.update('TUNIT7','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE8','PDCSAP_FLUX','column title: PDC photometry flux')
        hdu1.header.update('TFORM8','E','column format: float32')
        hdu1.header.update('TUNIT8','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE9','PDCSAP_FLUX_ERR','column title: PDC flux error')
        hdu1.header.update('TFORM9','E','column format: float32')
        hdu1.header.update('TUNIT9','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE10','SAP_QUALITY','column title: aperture photometry quality flag')
        hdu1.header.update('TFORM10','J','column format: signed integer32')
        hdu1.header.update('TTYPE11','PSF_CENTR1','column title: PSF fitted column centroid')
        hdu1.header.update('TFORM11','E','column format: float32')
        hdu1.header.update('TUNIT11','pixel','column units: pixel')
        hdu1.header.update('TTYPE12','PSF_CENTR1_ERR','column title: PSF fitted column error')
        hdu1.header.update('TFORM12','E','column format: float32')
        hdu1.header.update('TUNIT12','pixel','column units: pixel')
        hdu1.header.update('TTYPE13','PSF_CENTR2','column title: PSF fitted row centroid')
        hdu1.header.update('TFORM13','E','column format: float32')
        hdu1.header.update('TUNIT13','pixel','column units: pixel')
        hdu1.header.update('TTYPE14','PSF_CENTR2_ERR','column title: PSF fitted row error')
        hdu1.header.update('TFORM14','E','column format: float32')
        hdu1.header.update('TUNIT14','pixel','column units: pixel')
        hdu1.header.update('TTYPE15','MOM_CENTR1','column title: moment-derived column centroid')
        hdu1.header.update('TFORM15','E','column format: float32')
        hdu1.header.update('TUNIT15','pixel','column units: pixel')
        hdu1.header.update('TTYPE16','MOM_CENTR1_ERR','column title: moment-derived column error')
        hdu1.header.update('TFORM16','E','column format: float32')
        hdu1.header.update('TUNIT16','pixel','column units: pixel')
        hdu1.header.update('TTYPE17','MOM_CENTR2','column title: moment-derived row centroid')
        hdu1.header.update('TFORM17','E','column format: float32')
        hdu1.header.update('TUNIT17','pixel','column units: pixel')
        hdu1.header.update('TTYPE18','MOM_CENTR2_ERR','column title: moment-derived row error')
        hdu1.header.update('TFORM18','E','column format: float32')
        hdu1.header.update('TUNIT18','pixel','column units: pixel')
        hdu1.header.update('TTYPE19','POS_CORR1','column title: col correction for vel. abbern')
        hdu1.header.update('TFORM19','E','column format: float32')
        hdu1.header.update('TUNIT19','pixel','column units: pixel')
        hdu1.header.update('TTYPE20','POS_CORR2','column title: row correction for vel. abbern')
        hdu1.header.update('TFORM20','E','column format: float32')
        hdu1.header.update('TUNIT20','pixel','column units: pixel')
        hdu1.header.update('TTYPE21','RAW_FLUX','column title: raw aperture photometry flux')
        hdu1.header.update('TFORM21','E','column format: float32')
        hdu1.header.update('TUNIT21','e-/s','column units: electrons per second')
        hdu1.header.update('EXTNAME','LIGHTCURVE','name of extension')
        for i in range(len(cards1)):
            if (cards1[i].key not in hdu1.header.keys() and
                cards1[i].key[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC']):
                hdu1.header.update(cards1[i].key, cards1[i].value, cards1[i].comment)
        outstr.append(hdu1)

# construct output mask bitmap extension

    if status == 0:
        hdu2 = ImageHDU(maskmap)
        for i in range(len(cards2)):
            if cards2[i].key not in hdu2.header.keys():
                hdu2.header.update(cards2[i].key, cards2[i].value, cards2[i].comment)
            else:
                hdu2.header.cards[cards2[i].key].comment = cards2[i].comment
        outstr.append(hdu2)

# write output file

    if status == 0:
        outstr.writeto(outfile,checksum=True)

# close input structure

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

# end time

    kepmsg.clock('KEPEXTRACT finished at',logfile,verbose)

# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description=
        'Derive a light curve from a target pixel file, with user-defined apertures')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('maskfile', help='Name of mask defintion ASCII file', type=str)
    parser.add_argument('outfile', help='Name of output light curve FITS file', type=str)
    parser.add_argument('--background', action='store_true', help='Subtract background from data?')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', 
        dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    kepextract(args.infile,args.maskfile,args.outfile, args.background, args.clobber, args.verbose,
               args.logfile, args.status)
    
    

else:
    from pyraf import iraf
    
    
    parfile = iraf.osfn("kepler$kepextract.par")
    t = iraf.IrafTaskFactory(taskname="kepextract", value=parfile, function=kepextract)

