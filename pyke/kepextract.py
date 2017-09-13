from .utils import PyKEArgumentHelpFormatter
import math
import numpy as np
from astropy.io import fits as pyfits
from scipy.optimize import leastsq
from copy import copy
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepstat, kepfunc


__all__ = ['kepextract']


def kepextract(infile, outfile=None, maskfile='ALL', bkg=False, psfcentroid=False,
               overwrite=False, verbose=False,
               logfile='kepextract.log'):
    """
    kepextract -- create a light curve from a target pixel file by summing
    user-selected pixels

    kepextract calculates simple aperture photometry, from a target pixel file,
    for a user-supplied set of pixels. The Kepler pipeline sums a specific set
    of pixels to produce the standard light curves delivered to users. Termed
    the optimal aperture, the default pixel set is designed to maximize the
    signal-to-noise ratio of the resulting light curve, optimizing for transit
    detection. This tool provides users with a straightforward capability to
    alter the summed pixel set. Applications include:

        * Use of all pixels in the aperture
        * The Kepler pipeline does not produce a light curve for sources observed with
          custom or dedicated pixel masks. The user can create a light curve for
          these sources using kepextract.
        * Construction of pixel light curves, in which the time series for a single
          pixel can be examined.
        * Light curves for extended sources which may be
          poorly sampled by the optimal aperture.

    Parameters
    ----------
    infile : str
        Filename for the target pixel file.
    outfile : str
        Filename for the output light curve. This product will be written to
        the same FITS format as archived light curves.
    maskfile : str
        This string can be one of three options:

            * 'ALL' tells the task to calculate principal components from all
              pixels within the pixel mask stored in the input file.
            * 'APER' tells the task to calculate principal components from only
              the pixels within the photometric aperture stored in the input file
              (e.g. only those pixels summed by the Kepler pipeline to produce
              the light curve archived at MAST.
            * A filename describing the desired photometric aperture. Such a
              file can be constructed using the kepmask or kepffi tools, or can
              be created manually using the format described in the documentation
              for those tools.
    bkg : bool
        Option to subtract an estimate of the background. Background is
        calculated by identifying the median pixel value for each exposure.
        This method requires an adequate number of pixels within the target
        mask that contain background and negligible source flux. Note that
        background has already been subtracted from calibrated Kepler Target
        Pixel Files, but not early campaign data from the K2 mission.
    psfcentroid : bool
        Measure the star's position by fitting a 2D Gaussian PSF to the pixels
        in the mask. This will populate values for PSF_CENTR1 (column position)
        and PSF_CENTR2 (row position) in the output file.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Option for verbose mode, in which informative messages and warnings
        print to the shell and a logfile.
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepextract kplr008256049-2010174085026_lpd-targ.fits --maskfile ALL

    One further can plot the resulted light curve by doing

    .. code-block:: python

        import matplotlib.pyplot as plt
        from astropy.io import fits

        f = fits.open('outlc.fits')
        plt.plot(f[1].data['TIME'], f[1].data['SAP_FLUX'])

    or

    .. code-block:: bash

        $ kepdraw outlc.fits

    .. image:: ../_static/images/api/kepextract.png
        :align: center
    """
    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPEXTRACT -- '
            + ' infile={}'.format(infile)
            + ' maskfile={}'.format(maskfile)
            + ' outfile={}'.format(outfile)
            + ' background={}'.format(bkg)
            + ' psfcentroid={}'.format(psfcentroid)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPEXTRACT started at',logfile,verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPEXTRACT: {} exists. Use --overwrite'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile, mode='readonly', memmap=True)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    # input file data
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards
    table = instr[1].data[:]
    maskmap = copy(instr[2].data)

    kepmsg.log(logfile,"Extracting information from Target Pixel File...",verbose)

    # input table data
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, time = \
        kepio.readTPF(infile, 'TIME', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, timecorr = \
        kepio.readTPF(infile, 'TIMECORR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, cadenceno = \
        kepio.readTPF(infile, 'CADENCENO', logfile, verbose)
    cadenceno = np.array(cadenceno, dtype='int')

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, raw_cnts = \
        kepio.readTPF(infile, 'RAW_CNTS', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux = \
        kepio.readTPF(infile, 'FLUX', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux_err = \
        kepio.readTPF(infile, 'FLUX_ERR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux_bkg = \
        kepio.readTPF(infile, 'FLUX_BKG', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux_bkg_err = \
        kepio.readTPF(infile, 'FLUX_BKG_ERR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, cosmic_rays = \
        kepio.readTPF(infile, 'COSMIC_RAYS', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, quality = \
        kepio.readTPF(infile, 'QUALITY', logfile, verbose)

    try:
        #  ---for FITS wave #2
        pos_corr1 = np.array(table.field('POS_CORR1'), dtype='float64')
    except:
        pos_corr1 = np.empty(len(time))
        # ---temporary before FITS wave #2
        pos_corr1[:] = np.nan
    try:
        #  ---for FITS wave #2
        pos_corr2 = np.array(table.field('POS_CORR2'), dtype='float64')
    except:
        pos_corr2 = np.empty(len(time))
        # ---temporary before FITS wave #2
        pos_corr2[:] = np.nan

    # dummy columns for output file
    psf_centr1 = np.empty(len(time))
    psf_centr1[:] = np.nan
    psf_centr2 = np.empty(len(time))
    psf_centr2[:] = np.nan
    mom_centr1 = np.empty(len(time))
    mom_centr1[:] = np.nan
    mom_centr2 = np.empty(len(time))
    mom_centr2[:] = np.nan
    psf_centr1_err = np.empty(len(time))
    psf_centr1_err[:] = np.nan
    psf_centr2_err = np.empty(len(time))
    psf_centr2_err[:] = np.nan
    mom_centr1_err = np.empty(len(time))
    mom_centr1_err[:] = np.nan
    mom_centr2_err = np.empty(len(time))
    mom_centr2_err[:] = np.nan

    # read mask definition file
    if 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        maskx = np.array([], 'int')
        masky = np.array([], 'int')
        lines = kepio.openascii(maskfile, 'r', logfile, verbose)
        for line in lines:
            line = line.strip().split('|')
            if len(line) == 6:
                y0 = int(line[3])
                x0 = int(line[4])
                line = line[5].split(';')
                for items in line:
                    try:
                        masky = np.append(masky,y0 + int(items.split(',')[0]))
                        maskx = np.append(maskx,x0 + int(items.split(',')[1]))
                    except:
                        continue
        kepio.closeascii(lines, logfile, verbose)
        if len(maskx) == 0 or len(masky) == 0:
            errmsg = 'ERROR -- KEPEXTRACT: {} contains no pixels.'.format(maskfile)
            kepmsg.err(logfile, errmsg, verbose)

    # subimage physical WCS data
    crpix1p = cards2['CRPIX1P'].value
    crpix2p = cards2['CRPIX2P'].value
    crval1p = cards2['CRVAL1P'].value
    crval2p = cards2['CRVAL2P'].value
    cdelt1p = cards2['CDELT1P'].value
    cdelt2p = cards2['CDELT2P'].value

    # define new subimage bitmap...
    if 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        aperx = np.array([], 'int')
        apery = np.array([], 'int')
        aperb = np.array([], 'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = np.append(aperx, crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = np.append(apery, crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i, j] == 0:
                    aperb = np.append(aperb, 0)
                else:
                    aperb = np.append(aperb, 1)
                    maskmap[i, j] = 1
                    for k in range(len(maskx)):
                        if aperx[-1] == maskx[k] and apery[-1] == masky[k]:
                            aperb[-1] = 3
                            maskmap[i, j] = 3

    # trap case where no aperture needs to be defined but pixel positions are
    # still required for centroiding
    if maskfile.lower() == 'all':
        aperx = np.array([], 'int')
        apery = np.array([], 'int')
        aperb = np.array([], 'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = np.append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = np.append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i, j] == 0:
                    aperb = np.append(aperb, 0)
                else:
                    aperb = np.append(aperb, 3)
                    maskmap[i, j] = 3

    # ...or use old subimage bitmap
    if 'aper' in maskfile.lower():
        aperx = np.array([], 'int')
        apery = np.array([], 'int')
        aperb = np.array([], 'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperb = np.append(aperb, maskmap[i, j])
                aperx = np.append(aperx, crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = np.append(apery, crval2p + (i + 1 - crpix2p) * cdelt2p)

    # subtract median pixel value for background?
    sky = np.zeros(len(time), 'float32')
    if bkg:
        for i in range(len(time)):
            sky[i] = np.nanmedian(flux[i, :])

    # legal mask defined?
    if len(aperb) == 0:
        errmsg = ('ERROR -- KEPEXTRACT: no legal pixels within the subimage'
                  ' are defined.')
        kepmsg.err(logfile, errmsg, verbose)

    # construct new table flux data
    naper = (aperb == 3).sum()
    ntime = len(time)
    sap_flux = np.array([], 'float32')
    sap_flux_err = np.array([], 'float32')
    sap_bkg = np.array([], 'float32')
    sap_bkg_err = np.array([], 'float32')
    raw_flux = np.array([],'float32')
    kepmsg.log(logfile,"Aperture photometry...",verbose)
    for i in tqdm(range(len(time)), desc="Aperture photometry"):
        work1 = np.array([], 'float64')
        work2 = np.array([], 'float64')
        work3 = np.array([], 'float64')
        work4 = np.array([], 'float64')
        work5 = np.array([], 'float64')
        for j in range(len(aperb)):
            if aperb[j] == 3:
                work1 = np.append(work1, flux[i, j] - sky[i])
                work2 = np.append(work2, flux_err[i, j])
                work3 = np.append(work3, flux_bkg[i, j])
                work4 = np.append(work4, flux_bkg_err[i, j])
                work5 = np.append(work5, raw_cnts[i, j])
        sap_flux = np.append(sap_flux, np.nansum(work1))
        sap_flux_err = np.append(sap_flux_err, math.sqrt(np.nansum(work2 * work2)))
        sap_bkg = np.append(sap_bkg, np.nansum(work3))
        sap_bkg_err = np.append(sap_bkg_err, math.sqrt(np.nansum(work4 * work4)))
        raw_flux = np.append(raw_flux, np.nansum(work5))

    kepmsg.log(logfile,"Sample moments...",verbose)
    # construct new table moment data
    for i in tqdm(range(ntime), desc="Computing moments"):
        xf = np.zeros(shape=(naper))
        yf = np.zeros(shape=(naper))
        f = np.zeros(shape=(naper))
        xfe = np.zeros(shape=(naper))
        yfe = np.zeros(shape=(naper))
        fe = np.zeros(shape=(naper))
        k = -1
        for j in range(len(aperb)):
            if aperb[j] == 3:
                k += 1
                xf[k] = aperx[j] * flux[i, j]
                xfe[k] = aperx[j] * flux_err[i, j]
                yf[k] = apery[j] * flux[i, j]
                yfe[k] = apery[j] * flux_err[i, j]
                f[k] = flux[i, j]
                fe[k] = flux_err[i, j]
        xfsum = np.nansum(xf)
        yfsum = np.nansum(yf)
        fsum = np.nansum(f)
        xfsume = math.sqrt(np.nansum(xfe * xfe) / naper)
        yfsume = math.sqrt(np.nansum(yfe * yfe) / naper)
        fsume = math.sqrt(np.nansum(fe * fe) / naper)
        # Ignore "RuntimeWarning: invalid value encountered in double_scalars"
        with np.errstate(divide='ignore', invalid='ignore'):
            mom_centr1[i] = xfsum / fsum
            mom_centr2[i] = yfsum / fsum
            mom_centr1_err[i] = math.sqrt((xfsume / xfsum) ** 2 + ((fsume / fsum) ** 2))
            mom_centr2_err[i] = math.sqrt((yfsume / yfsum) ** 2 + ((fsume / fsum) ** 2))
    mom_centr1_err = mom_centr1_err * mom_centr1
    mom_centr2_err = mom_centr2_err * mom_centr2

    if psfcentroid:
        kepmsg.log(logfile,"PSF Centroiding...",verbose)
        # construct new table PSF data
        psf_centr1 = np.zeros(shape=(ntime))
        psf_centr2 = np.zeros(shape=(ntime))
        psf_centr1_err = np.zeros(shape=(ntime))
        psf_centr2_err = np.zeros(shape=(ntime))
        modx = np.zeros(shape=(naper))
        mody = np.zeros(shape=(naper))
        k = -1
        for j in range(len(aperb)):
            if (aperb[j] == 3):
                k += 1
                modx[k] = aperx[j]
                mody[k] = apery[j]
        for i in tqdm(range(ntime), desc='PSF centroiding'):
            modf = np.zeros(shape=(naper))
            k = -1
            guess = [mom_centr1[i], mom_centr2[i], np.nanmax(flux[i:]), 1.0, 1.0, 0.0, 0.0]
            for j in range(len(aperb)):
                if (aperb[j] == 3):
                    k += 1
                    modf[k] = flux[i,j]
                    args = (modx, mody, modf)
            try:
                ans = leastsq(kepfunc.PRFgauss2d, guess, args=args, xtol=1.0e-8,
                              ftol=1.0e-4, full_output=True)
                s_sq = (ans[2]['fvec'] ** 2).sum() / (ntime - len(guess))
                psf_centr1[i] = ans[0][0]
                psf_centr2[i] = ans[0][1]
            except:
                pass
            try:
                psf_centr1_err[i] = sqrt(diag(ans[1] * s_sq))[0]
            except:
                psf_centr1_err[i] = np.nan
            try:
                psf_centr2_err[i] = sqrt(diag(ans[1] * s_sq))[1]
            except:
                psf_centr2_err[i] = np.nan

    # construct output primary extension
    hdu0 = pyfits.PrimaryHDU()
    for i in range(len(cards0)):
        if cards0[i].keyword not in hdu0.header.keys():
            hdu0.header[cards0[i].keyword] = (cards0[i].value, cards0[i].comment)
        else:
            hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
    kepkey.history(call, hdu0, outfile, logfile, verbose)
    outstr = pyfits.HDUList(hdu0)

    # construct output light curve extension
    col1 = pyfits.Column(name='TIME', format='D', unit='BJD - 2454833',
                         array=time)
    col2 = pyfits.Column(name='TIMECORR', format='E', unit='d', array=timecorr)
    col3 = pyfits.Column(name='CADENCENO', format='J', array=cadenceno)
    col4 = pyfits.Column(name='SAP_FLUX', format='E', array=sap_flux)
    col5 = pyfits.Column(name='SAP_FLUX_ERR', format='E', array=sap_flux_err)
    col6 = pyfits.Column(name='SAP_BKG', format='E', array=sap_bkg)
    col7 = pyfits.Column(name='SAP_BKG_ERR', format='E', array=sap_bkg_err)
    col8 = pyfits.Column(name='PDCSAP_FLUX', format='E', array=sap_flux)
    col9 = pyfits.Column(name='PDCSAP_FLUX_ERR', format='E',
                         array=sap_flux_err)
    col10 = pyfits.Column(name='SAP_QUALITY', format='J', array=quality)
    col11 = pyfits.Column(name='PSF_CENTR1', format='E', unit='pixel',
                          array=psf_centr1)
    col12 = pyfits.Column(name='PSF_CENTR1_ERR', format='E', unit='pixel',
                          array=psf_centr1_err)
    col13 = pyfits.Column(name='PSF_CENTR2', format='E', unit='pixel',
                          array=psf_centr2)
    col14 = pyfits.Column(name='PSF_CENTR2_ERR', format='E', unit='pixel',
                          array=psf_centr2_err)
    col15 = pyfits.Column(name='MOM_CENTR1', format='E', unit='pixel',
                          array=mom_centr1)
    col16 = pyfits.Column(name='MOM_CENTR1_ERR', format='E', unit='pixel',
                          array=mom_centr1_err)
    col17 = pyfits.Column(name='MOM_CENTR2', format='E', unit='pixel',
                          array=mom_centr2)
    col18 = pyfits.Column(name='MOM_CENTR2_ERR', format='E', unit='pixel',
                          array=mom_centr2_err)
    col19 = pyfits.Column(name='POS_CORR1', format='E', unit='pixel',
                          array=pos_corr1)
    col20 = pyfits.Column(name='POS_CORR2', format='E', unit='pixel',
                          array=pos_corr2)
    col21 = pyfits.Column(name='RAW_FLUX', format='E', array=raw_flux)
    cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8,
                           col9, col10, col11, col12, col13, col14, col15,
                           col16, col17, col18, col19, col20, col21])
    hdu1 = pyfits.BinTableHDU.from_columns(cols)
    hdu1.header['TTYPE1'] = ('TIME', 'column title: data time stamps')
    hdu1.header['TFORM1'] = ('D', 'data type: float64')
    hdu1.header['TUNIT1'] = ('BJD - 2454833',
                             'column units: barycenter corrected JD')
    hdu1.header['TDISP1'] = ('D12.7', 'column display format')
    hdu1.header['TTYPE2'] = ('TIMECORR',
                             'column title: barycentric-timeslice correction')
    hdu1.header['TFORM2'] = ('E', 'data type: float32')
    hdu1.header['TUNIT2'] = ('d', 'column units: days')
    hdu1.header['TTYPE3'] = ('CADENCENO',
                             'column title: unique cadence number')
    hdu1.header['TFORM3'] = ('J', 'column format: signed integer32')
    hdu1.header['TTYPE4'] = ('SAP_FLUX',
                             'column title: aperture photometry flux')
    hdu1.header['TFORM4'] = ('E', 'column format: float32')
    hdu1.header['TUNIT4'] = ('e-/s', 'column units: electrons per second')
    hdu1.header['TTYPE5'] = ('SAP_FLUX_ERR',
                             'column title: aperture phot. flux error')
    hdu1.header['TFORM5'] = ('E', 'column format: float32')
    hdu1.header['TUNIT5'] = ('e-/s',
                             'column units: electrons per second (1-sigma)')
    hdu1.header['TTYPE6'] = ('SAP_BKG',
                             'column title: aperture phot. background flux')
    hdu1.header['TFORM6'] = ('E','column format: float32')
    hdu1.header['TUNIT6'] = ('e-/s', 'column units: electrons per second')
    hdu1.header['TTYPE7'] = ('SAP_BKG_ERR',
                             'column title: ap. phot. background flux error')
    hdu1.header['TFORM7'] = ('E', 'column format: float32')
    hdu1.header['TUNIT7'] = ('e-/s',
                             'column units: electrons per second (1-sigma)')
    hdu1.header['TTYPE8'] = ('PDCSAP_FLUX',
                             'column title: PDC photometry flux')
    hdu1.header['TFORM8'] = ('E', 'column format: float32')
    hdu1.header['TUNIT8'] = ('e-/s', 'column units: electrons per second')
    hdu1.header['TTYPE9'] = ('PDCSAP_FLUX_ERR', 'column title: PDC flux error')
    hdu1.header['TFORM9'] = ('E', 'column format: float32')
    hdu1.header['TUNIT9'] = ('e-/s',
                             'column units: electrons per second (1-sigma)')
    hdu1.header['TTYPE10'] = ('SAP_QUALITY',
                              'column title: aperture photometry quality flag')
    hdu1.header['TFORM10'] = ('J', 'column format: signed integer32')
    hdu1.header['TTYPE11'] = ('PSF_CENTR1',
                              'column title: PSF fitted column centroid')
    hdu1.header['TFORM11'] = ('E', 'column format: float32')
    hdu1.header['TUNIT11'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE12'] = ('PSF_CENTR1_ERR',
                              'column title: PSF fitted column error')
    hdu1.header['TFORM12'] = ('E', 'column format: float32')
    hdu1.header['TUNIT12'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE13'] = ('PSF_CENTR2',
                              'column title: PSF fitted row centroid')
    hdu1.header['TFORM13'] = ('E', 'column format: float32')
    hdu1.header['TUNIT13'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE14'] = ('PSF_CENTR2_ERR',
                              'column title: PSF fitted row error')
    hdu1.header['TFORM14'] = ('E', 'column format: float32')
    hdu1.header['TUNIT14'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE15'] = ('MOM_CENTR1',
                              'column title: moment-derived column centroid')
    hdu1.header['TFORM15'] = ('E', 'column format: float32')
    hdu1.header['TUNIT15'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE16'] = ('MOM_CENTR1_ERR',
                              'column title: moment-derived column error')
    hdu1.header['TFORM16'] = ('E', 'column format: float32')
    hdu1.header['TUNIT16'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE17'] = ('MOM_CENTR2',
                              'column title: moment-derived row centroid')
    hdu1.header['TFORM17'] = ('E', 'column format: float32')
    hdu1.header['TUNIT17'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE18'] = ('MOM_CENTR2_ERR',
                              'column title: moment-derived row error')
    hdu1.header['TFORM18'] = ('E', 'column format: float32')
    hdu1.header['TUNIT18'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE19'] = ('POS_CORR1',
                              'column title: col correction for vel. abbern')
    hdu1.header['TFORM19'] = ('E', 'column format: float32')
    hdu1.header['TUNIT19'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE20'] = ('POS_CORR2',
                              'column title: row correction for vel. abbern')
    hdu1.header['TFORM20'] = ('E', 'column format: float32')
    hdu1.header['TUNIT20'] = ('pixel', 'column units: pixel')
    hdu1.header['TTYPE21'] = ('RAW_FLUX',
                              'column title: raw aperture photometry flux')
    hdu1.header['TFORM21'] = ('E', 'column format: float32')
    hdu1.header['TUNIT21'] = ('e-/s', 'column units: electrons per second')
    hdu1.header['EXTNAME'] = ('LIGHTCURVE', 'name of extension')

    for i in range(len(cards1)):
        if (cards1[i].keyword not in hdu1.header.keys() and
            cards1[i].keyword[:4] not in ['TTYP', 'TFOR', 'TUNI', 'TDIS',
                                          'TDIM', 'WCAX', '1CTY', '2CTY',
                                          '1CRP', '2CRP', '1CRV', '2CRV',
                                          '1CUN', '2CUN', '1CDE', '2CDE',
                                          '1CTY', '2CTY', '1CDL', '2CDL',
                                          '11PC', '12PC', '21PC', '22PC']):
            hdu1.header[cards1[i].keyword] = (cards1[i].value,
                                              cards1[i].comment)
    outstr.append(hdu1)

    # construct output mask bitmap extension
    hdu2 = pyfits.ImageHDU(maskmap)
    for i in range(len(cards2)):
        if cards2[i].keyword not in hdu2.header.keys():
            hdu2.header[cards2[i].keyword] = (cards2[i].value,
                                              cards2[i].comment)
        else:
            hdu2.header.cards[cards2[i].keyword].comment = cards2[i].comment
    outstr.append(hdu2)

    # write output file
    kepmsg.log(logfile,"Writing output file {}...".format(outfile),verbose)
    outstr.writeto(outfile, checksum=True)
    # close input structure
    instr.close()
    # end time
    kepmsg.clock('KEPEXTRACT finished at', logfile, verbose)

def kepextract_main():
    import argparse

    parser = argparse.ArgumentParser(
            description=('Derive a light curve from a target pixel file, with'
                         ' user-defined apertures'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input target pixel file',
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepextract.'),
                        default=None)
    parser.add_argument('--maskfile', default='ALL',
                        help='Name of mask defintion ASCII file',
                        type=str)
    parser.add_argument('--bkg', action='store_true',
                        help='Subtract background from data?')
    parser.add_argument('--psfcentroid', action='store_true',
                        help=("Measure the star's position by fitting a 2D"
                              " Gaussian PSF to the pixels in the mask. This"
                              " will populate values for PSF_CENTR1 (column"
                              " position) and PSF_CENTR2 (row position) in the"
                              " output file."))
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepextract.log', type=str)
    args = parser.parse_args()
    kepextract(args.infile, args.outfile, args.maskfile, args.bkg,
               args.psfcentroid, args.overwrite, args.verbose,
               args.logfile)
