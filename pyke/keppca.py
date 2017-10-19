from .utils import PyKEArgumentHelpFormatter
from . import kepmsg, kepio, kepkey, kepplot
import re
import numpy as np
from astropy.io import fits as pyfits
from scipy import optimize as opt
from matplotlib import pyplot as plt
from tqdm import tqdm
import random


__all__ = ['keppca']


def keppca(infile, outfile=None, maskfile='ALL', components='1-3', plotpca=False,
           nmaps=10, overwrite=False, verbose=False, logfile='keppca.log'):
    """
    keppca -- Perform principal component analysis upon a target pixel file

    keppca provides a method to mitigate for motion-derived systematic
    artifacts via Principle Component Analysis (PCA). This method was
    demonstrated on Kepler light curves by Harrison et al. (2012). It provides
    an alternative to cotrending data using basis vectors (kepcotrend) and
    correlating aperture photometry struture with time-series centroid
    measurements (kepsff). PCA will perhaps become a more widespread tool in
    the K2 era where the magnitde of target motion across the detector over a
    Kepler quarter is experienced by a K2 target over just 6-hours during its
    regular sequence of thruster firings that counteract boresight roll motion
    Pixel-level PCA employs only those pixels collected around a specific
    target and separates photometric trends common to all pixels from trends
    localized to individual targets or pixels in a series of principal
    component curves.

    The user has the option to choose the specific set of pixels to sample in
    this analysis. Principal components are plotted by the tool and written out
    to an output FITS file in an output extension called PRINCIPAL_COMPONENTS.
    The extension contains a 2D table with one row per timestamp recorded in
    the input file and one column for every principal component. Summing all
    principal components together will reconstruct a normalized version of the
    summed pixel within the chosen aperture. The user also has the choice of
    which principal components to optimally-subtract from the aperture-derived
    light curve in order to remove motion systematics from the time-series
    data. The aperture light curve and the corrected light curve are written to
    the LIGHTCURVE extension of the output file. The first populates the
    SAP_FLUX data column and the second is written to a column called PCA_FLUX.
    This output file can be used as input for other PyKE tasks and can be e.g.
    inspected using kepdraw.

    Parameters
    ----------
    infile : str
        The name of a standard format FITS file containing Kepler or K2 target
        pixels within the first data extension.
    outfile : str
        Filename for the output light curves and principal components. This
        product will be written to the same FITS format as archived light
        curves. Aperture photometry will be stored in the SAP_FLUX column of
        the first FITS extension called LIGHTCURVE. A version of this light
        curve with principal components subtracted is stored in column PCA_FLUX
        and a normalized version is stored in PCA_FLUX_NRM. The individual
        principal components are stored within a new FITS extension called
        PRINCIPAL_COMPONENTS.
    maskfile : str
        This string can be one of three options:

        * 'ALL' tells the task to calculate principal components from all
          pixels within the pixel mask stored in the input file.

        * 'APER' tells the task to calculate principal components from only the
          pixels within the photometric aperture stored in the input file (e.g.
          only those pixels summed by the Kepler pipeline to produce the light
          curve archived at MAST (note that no such light curves are currently
          being created for the K2 mission)

        * A filename describing the desired photometric aperture. Such a file
          can be constructed using the kepmask or kepffi tools, or can be created
          manually using the format described in the documentation for those
          tools. Note that if an aperture provided is not stricly rectangular,
          keppca will increase the size of the aperture so that it defines the
          smallest possible rectangle that contains all of the specified pixels.

    components : str
        A list of the principal components to subtract from the aperture light
        curve. The strings '1 2 3 4 5', 1,'2,3,4,5' and '1,2,3-5' yield the
        same result.
    plotpca : bool
        If True, keppca will produce plots containing individual principal
        components, correlation maps and light curves, both aperture and
        PCA-corrected versions. The will be stored as hardcopies in PNG format.
    nmaps : int
        The number of correlation maps and principal components to plot as
        output. This can be any positive integer up to the number of pixels
        within the mask, although note that many hundreds of plots will likely
        become prohibitive and is unlikely to be informative.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning message

    Examples
    --------
    .. code-block:: bash

        $ keppca ktwo202073445-c00_lpd-targ.fits.gz --plotpca

    .. image:: ../_static/images/api/keppca.png
        :align: center
    """

    import mdp

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPPCA -- '
            + ' infile={}'.format(infile)
            + ' maskfile={}'.format(maskfile)
            + ' outfile={}'.format(outfile)
            + ' components={}'.format(components)
            + ' plotpca={}'.format(plotpca)
            + ' nmaps={}'.format(nmaps)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call + '\n', verbose)

    kepmsg.clock('KEPPCA started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = ('ERROR -- KEPPCA: {} exists. Use overwrite=True'
                  .format(outfile))
        kepmsg.err(logfile, errmsg, verbose)

    # Set output file names - text file with data and plot
    dataout = np.copy(outfile)
    repname = re.sub('.fits', '.png', outfile)

    # open input file
    instr = pyfits.open(infile, mode='readonly', memmap=True)
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)

    # open TPF FITS file
    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, barytime = \
    kepio.readTPF(infile, 'TIME', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, tcorr = \
    kepio.readTPF(infile, 'TIMECORR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, cadno = \
    kepio.readTPF(infile, 'CADENCENO', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, fluxpixels = \
    kepio.readTPF(infile, 'FLUX', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, errpixels = \
    kepio.readTPF(infile, 'FLUX_ERR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux_bkg = \
    kepio.readTPF(infile, 'FLUX_BKG', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, flux_bkg_err = \
    kepio.readTPF(infile, 'FLUX_BKG_ERR', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, qual = \
    kepio.readTPF(infile, 'QUALITY', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, pcorr1 = \
    kepio.readTPF(infile, 'POS_CORR1', logfile, verbose)

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, pcorr2 = \
    kepio.readTPF(infile, 'POS_CORR2', logfile ,verbose)

    # Save original data dimensions, in case of using maskfile
    xdimorig = xdim
    ydimorig = ydim

    # read mask definition file if it has been supplied
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
                        masky = np.append(masky, y0 + int(items.split(',')[0]))
                        maskx = np.append(maskx, x0 + int(items.split(',')[1]))
                    except:
                        continue
        kepio.closeascii(lines, logfile, verbose)
        if len(maskx) == 0 or len(masky) == 0:
            errmsg = 'ERROR -- KEPPCA: {} contains no pixels.'.format(maskfile)
            kepmsg.err(logfile, errmsg, verbose)
        xdim = max(maskx) - min(maskx) + 1   # Find largest x dimension of mask
        ydim = max(masky) - min(masky) + 1   # Find largest y dimension of mask

        # pad mask to ensure it is rectangular
        workx = np.array([], 'int')
        worky = np.array([], 'int')
        for ip in np.arange(min(maskx), max(maskx) + 1):
            for jp in np.arange(min(masky), max(masky) + 1):
                workx = np.append(workx, ip)
                worky = np.append(worky, jp)
        maskx = workx
        masky = worky

    # define new subimage bitmap...
    if maskfile.lower() != 'all':
        aperx = np.array([], 'int')
        apery = np.array([], 'int')
        # aperb is an array that contains the pixel numbers in the mask
        aperb = maskx - x0 + xdimorig * (masky - y0)
        npix = len(aperb)

    # ...or use all pixels
    if maskfile.lower() == 'all':
        npix = xdimorig * ydimorig
        aperb = np.array([], 'int')
        aperb = np.r_[0: npix]

    # legal mask defined?
    if len(aperb) == 0:
        message = ('ERROR -- KEPPCA: no legal pixels within the subimage are'
                   ' defined.')
        kepmsg.err(logfile, message, verbose)

    # Identify principal components desired
    pcaout = []
    txt = components.strip().split(',')
    for work1 in txt:
        try:
            pcaout.append(int(work1.strip()))
        except:
            work2 = work1.strip().split('-')
            try:
                for work3 in range(int(work2[0]), int(work2[1]) + 1):
                    pcaout.append(work3)
            except:
                errmsg = ('ERROR -- KEPPCA: cannot understand principal'
                          ' component list requested')
                kepmsg.err(logfile, message, verbose)

    pcaout = set(np.sort(pcaout))
    # The list of pca component numbers to be removed
    pcarem = np.array(list(pcaout)) - 1

    # Initialize arrays and variables, and apply pixel mask to the data
    ntim = 0
    time = np.array([], dtype='float64')
    timecorr = np.array([], dtype='float32')
    cadenceno = np.array([], dtype='int')
    pixseries = np.array([], dtype='float32')
    errseries = np.array([], dtype='float32')
    bkgseries = np.array([], dtype='float32')
    berseries = np.array([], dtype='float32')
    quality = np.array([], dtype='float32')
    pos_corr1 = np.array([], dtype='float32')
    pos_corr2 = np.array([], dtype='float32')
    nrows = np.size(fluxpixels, 0)

    # Apply the pixel mask so we are left with only the desired pixels
    pixseriesb = fluxpixels[:, aperb]
    errseriesb = errpixels[:, aperb]
    bkgseriesb = flux_bkg[:, aperb]
    berseriesb = flux_bkg_err[:, aperb]

    # Read in the data to various arrays
    for i in range(nrows):
        if (qual[i] < 10000 and np.isfinite(barytime[i])
            and np.isfinite(fluxpixels[i, int(ydim * xdim / 2 + 0.5)])
            and np.isfinite(fluxpixels[i, 1 + int(ydim * xdim / 2 + 0.5)])):
            ntim += 1
            time = np.append(time, barytime[i])
            timecorr = np.append(timecorr, tcorr[i])
            cadenceno = np.append(cadenceno, cadno[i])
            pixseries = np.append(pixseries, pixseriesb[i])
            errseries = np.append(errseries, errseriesb[i])
            bkgseries = np.append(bkgseries, bkgseriesb[i])
            berseries = np.append(berseries, berseriesb[i])
            quality = np.append(quality, qual[i])
            pos_corr1 = np.append(pos_corr1, pcorr1[i])
            pos_corr2 = np.append(pos_corr2, pcorr2[i])
    pixseries = np.reshape(pixseries,(ntim, npix))
    errseries = np.reshape(errseries,(ntim, npix))
    bkgseries = np.reshape(bkgseries,(ntim, npix))
    berseries = np.reshape(berseries,(ntim, npix))
    tmp =  np.ma.median(np.ma.masked_invalid(pixseries), axis=1)

    for i in range(len(tmp)):
         pixseries[i] = pixseries[i] - tmp[i]
    pixseries = np.ma.masked_invalid(pixseries)

    # Figure out which pixels are undefined/nan and remove them.
    # Keep track for adding back in later
    nanpixels = np.array([], dtype='int')
    i = 0
    while i < npix:
        if np.isnan(pixseries[0, i]):
            nanpixels = np.append(nanpixels, i)
            npix = npix - 1
        i = i + 1
    pixseries = np.delete(pixseries, nanpixels, 1)
    errseries = np.delete(errseries, nanpixels, 1)
    pixseries[np.isnan(pixseries)] = random.gauss(100, 10)
    errseries[np.isnan(errseries)] = 10

    # Compute statistical weights, means, standard deviations
    weightseries = (pixseries / errseries) ** 2
    pixMean = np.average(pixseries, axis=0, weights=weightseries)
    pixStd  = np.std(pixseries, axis=0)

    # Normalize the input by subtracting the mean and divising by the standard
    # deviation.
    # This makes it a correlation-based PCA, which is what we want.
    pixseriesnorm = (pixseries - pixMean) / pixStd

    # Number of principal components to compute. Setting it equal to the number
    # of pixels
    nvecin = npix

    # Run PCA using the MDP Whitening PCA, which produces normalized PCA
    # components (zero mean and unit variance)
    pcan = mdp.nodes.WhiteningNode(svd=True)
    pcar = pcan.execute(pixseriesnorm)
    eigvec = pcan.get_recmatrix()
    model = pcar

    # Re-insert nan columns as zeros
    for i in range(len(nanpixels)):
        nanpixels[i] = nanpixels[i] - i
    eigvec = np.insert(eigvec, nanpixels, 0, 1)
    pixMean = np.insert(pixMean, nanpixels, 0, 0)

    #  Make output eigenvectors (correlation images) into xpix by ypix images
    eigvec = eigvec.reshape(nvecin, ydim, xdim)
    # Calculate sum of all pixels to display as raw lightcurve and other quantities
    pixseriessum = np.sum(pixseries, axis=1)
    # Number of components to remove
    nrem = len(pcarem)
    # Number of pcas to plot - currently set to plot all components, but could set
    # nplot = nrem to just plot as many components as is being removed
    nplot = npix

    # Subtract components by fitting them to the summed light curve
    x0 = np.tile(-1.0, 1)
    for k in tqdm(range(nrem)):
        def f(x):
            fluxcor = pixseriessum
            for k in range(len(x)):
                fluxcor = fluxcor - x[k]*model[:, pcarem[k]]
            return mad(fluxcor)
        if k == 0:
            x0 = np.array([-1.0])
        else:
            x0 = np.append(x0, 1.0)
        myfit = opt.fmin(f, x0, maxiter=50000, maxfun=50000, disp=False)
        x0 = myfit

    # Now that coefficients for all components have been found, subtract them
    # to produce a calibrated time-series,
    # and then divide by the robust mean to produce a normalized time series
    # as well
    c = myfit
    fluxcor = pixseriessum
    for k in range(0, nrem):
        fluxcor = fluxcor - c[k] * model[:, pcarem[k]]

    normfluxcor = fluxcor / np.nanmean(reject_outliers(fluxcor, 2))

    # input file data
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards
    table = instr[1].data[:]
    maskmap = np.copy(instr[2].data)

    # subimage physical WCS data
    crpix1p = cards2['CRPIX1P'].value
    crpix2p = cards2['CRPIX2P'].value
    crval1p = cards2['CRVAL1P'].value
    crval2p = cards2['CRVAL2P'].value
    cdelt1p = cards2['CDELT1P'].value
    cdelt2p = cards2['CDELT2P'].value

    # dummy columns for output file
    sap_flux_err = np.empty(len(time))
    sap_flux_err[:] = np.nan
    sap_bkg = np.empty(len(time))
    sap_bkg[:] = np.nan
    sap_bkg_err = np.empty(len(time))
    sap_bkg_err[:] = np.nan
    pdc_flux = np.empty(len(time))
    pdc_flux[:] = np.nan
    pdc_flux_err = np.empty(len(time))
    pdc_flux_err[:] = np.nan
    psf_centr1 = np.empty(len(time))
    psf_centr1[:] = np.nan
    psf_centr1_err = np.empty(len(time))
    psf_centr1_err[:] = np.nan
    psf_centr2 = np.empty(len(time))
    psf_centr2[:] = np.nan
    psf_centr2_err = np.empty(len(time))
    psf_centr2_err[:] = np.nan
    mom_centr1 = np.empty(len(time))
    mom_centr1[:] = np.nan
    mom_centr1_err = np.empty(len(time))
    mom_centr1_err[:] = np.nan
    mom_centr2 = np.empty(len(time))
    mom_centr2[:] = np.nan
    mom_centr2_err = np.empty(len(time))
    mom_centr2_err[:] = np.nan

    # mask bitmap
    if 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i, j] == 0:
                    pass
                else:
                    maskmap[i, j] = 1
                    for k in range(len(maskx)):
                        if aperx[-1] == maskx[k] and apery[-1] == masky[k]:
                            maskmap[i, j] = 3

    # construct output primary extension
    hdu0 = pyfits.PrimaryHDU()
    for i in range(len(cards0)):
        if cards0[i].keyword not in hdu0.header.keys():
            hdu0.header[cards0[i].keyword] = (cards0[i].value,
                                              cards0[i].comment)
        else:
            hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
    kepkey.history(call, hdu0, outfile, logfile, verbose)
    outstr = pyfits.HDUList(hdu0)

    # construct output light curve extension
    col1 = pyfits.Column(name='TIME', format='D', unit='BJD - 2454833',
                         array=time)
    col2 = pyfits.Column(name='TIMECORR', format='E', unit='d', array=timecorr)
    col3 = pyfits.Column(name='CADENCENO', format='J', array=cadenceno)
    col4 = pyfits.Column(name='SAP_FLUX', format='E', unit='e-/s',
                         array=pixseriessum)
    col5 = pyfits.Column(name='SAP_FLUX_ERR', format='E', unit='e-/s',
                         array=sap_flux_err)
    col6 = pyfits.Column(name='SAP_BKG', format='E', unit='e-/s',
                         array=sap_bkg)
    col7 = pyfits.Column(name='SAP_BKG_ERR', format='E', unit='e-/s',
                         array=sap_bkg_err)
    col8 = pyfits.Column(name='PDCSAP_FLUX', format='E', unit='e-/s',
                         array=pdc_flux)
    col9 = pyfits.Column(name='PDCSAP_FLUX_ERR', format='E', unit='e-/s',
                         array=pdc_flux_err)
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
    col21 = pyfits.Column(name='PCA_FLUX', format='E', unit='e-/s',
                          array=fluxcor)
    col22 = pyfits.Column(name='PCA_FLUX_NRM', format='E', array=normfluxcor)
    cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8,
                           col9, col10, col11, col12, col13, col14, col15,
                           col16, col17, col18, col19, col20, col21, col22])
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
    hdu1.header['TFORM6'] = ('E', 'column format: float32')
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
    hdu1.header['TTYPE21'] = ('PCA_FLUX', 'column title: PCA-corrected flux')
    hdu1.header['TFORM21'] = ('E', 'column format: float32')
    hdu1.header['TUNIT21'] = ('pixel', 'column units: e-/s')
    hdu1.header['TTYPE22'] = ('PCA_FLUX_NRM',
                              'column title: normalized PCA-corrected flux')
    hdu1.header['TFORM22'] = ('E', 'column format: float32')
    hdu1.header['EXTNAME'] = ('LIGHTCURVE', 'name of extension')
    for i in range(len(cards1)):
        if (cards1[i].keyword not in hdu1.header.keys() and
            cards1[i].keyword[:4] not in ['TTYP', 'TFOR', 'TUNI', 'TDIS',
                                          'TDIM', 'WCAX', '1CTY', '2CTY',
                                          '1CRP', '2CRP', '1CRV', '2CRV',
                                          '1CUN', '2CUN', '1CDE', '2CDE',
                                          '1CTY', '2CTY', '1CDL', '2CDL',
                                          '11PC', '12PC', '21PC', '22PC']):
            hdu1.header[cards1[i].keyword] = (cards1[i].value, cards1[i].comment)
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

    # construct principal component table
    cols = [pyfits.Column(name='TIME', format='E', unit='BJD - 2454833',
                          array=time)]
    for i in range(len(pcar[0, :])):
        colname = 'PC' + str(i + 1)
        col = pyfits.Column(name=colname, format='E', array=pcar[:, i])
        cols.append(col)
    hdu3 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))
    hdu3.header['EXTNAME'] = ('PRINCIPAL_COMPONENTS', 'name of extension')
    hdu3.header['TTYPE1'] = ('TIME', 'column title: data time stamps')
    hdu3.header['TFORM1'] = ('D', 'data type: float64')
    hdu3.header['TUNIT1'] = ('BJD - 2454833',
                             'column units: barycenter corrected JD')
    hdu3.header['TDISP1'] = ('D12.7', 'column display format')
    for i in range(len(pcar[0, :])):
        hdu3.header['TTYPE' + str(i + 2)] = ("PC" + str(i + 1),
                                             "column title: principal "
                                             "component number " + str(i + 1))
        hdu3.header['TFORM' + str(i + 2)] = ('E', 'column format: float32')
    outstr.append(hdu3)

    # write output file
    print("Writing output file {}...".format(outfile))
    outstr.writeto(outfile)
    # close input structure
    instr.close()
    # Create PCA report
    if plotpca:
        npp = 7 # Number of plots per page
        l = 1
        repcnt = 1
        for k in range(nmaps):
            # First plot of every pagewith flux image,
            # flux and calibrated time series
            if (k % (npp - 1) == 0):
                plt.figure(figsize=[10, 16])
                plt.subplot2grid((npp,6), (0, 0), colspan=2)
                plt.imshow(np.log10(np.flipud(pixMean.reshape(ydim,xdim)) - min(pixMean) + 1),
                           interpolation="nearest",cmap='RdYlBu')
                plt.xticks([])
                plt.yticks([])
                ax1 = plt.subplot2grid((npp, 6), (0, 2), colspan=4)
                px = np.copy(time) + bjdref
                py = np.copy(pixseriessum)
                px, xlab = kepplot.cleanx(px, logfile, verbose)
                py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
                kepplot.RangeOfPlot(px, py, 0.01, False)
                kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2,
                               True)
                py = np.copy(fluxcor)
                py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
                plt.plot(px, py, marker='.', color='r', linestyle='',
                         markersize=1.0)
                kepplot.labels('', re.sub('\)', '', re.sub('Flux \(','', ylab)),
                               'k', 14)
                plt.grid()
                plt.setp(ax1.get_xticklabels(), visible=False)

            # plot principal components
            plt.subplot2grid((npp, 6), (l, 0), colspan=2)
            plt.imshow(eigvec[k], interpolation="nearest", cmap='RdYlBu')
            plt.xlim(-0.5, xdim-0.5)
            plt.ylim(-0.5, ydim-0.5)
            plt.xticks([])
            plt.yticks([])

            # The last plot on the page that should have the xlabel
            if (k % (npp - 1) == npp - 2 or k == nvecin - 1):
                plt.subplot2grid((npp, 6), (l, 2), colspan=4)
                py = np.copy(model[:, k])
                kepplot.RangeOfPlot(px, py, 0.01, False)
                kepplot.plot1d(px, py, cadence, 'r', 1.0, 'g', 0.2, True)
                kepplot.labels(xlab, 'PC ' + str(k + 1), 'k', 14)
                plt.grid()
                plt.tight_layout()
                l = 1
                plt.savefig(re.sub('.png', '_%d.png' % repcnt,repname))
                repcnt += 1
            # The other plots on the page that should have no xlabel
            else:
                ax2 = plt.subplot2grid((npp, 6), (l, 2), colspan=4)
                py = np.copy(model[:, k])
                kepplot.RangeOfPlot(px, py, 0.01, False)
                kepplot.plot1d(px, py, cadence, 'r', 1.0, 'g', 0.2, True)
                kepplot.labels('', 'PC ' + str(k + 1), 'k', 14)
                plt.grid()
                plt.setp(ax2.get_xticklabels(), visible=False)
                plt.tight_layout()
                l=l+1
        plt.savefig(re.sub('.png', '_%d.png' % repcnt, repname))

    # plot style and size
    if plotpca:
        plt.figure()
        plt.clf()
        # plot aperture photometry and PCA corrected data
        ax = kepplot.location([0.06, 0.54, 0.93, 0.43])
        px = np.copy(time) + bjdref
        py = np.copy(pixseriessum)
        px, xlab = kepplot.cleanx(px, logfile, verbose)
        py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
        kepplot.RangeOfPlot(px, py, 0.01, False)
        kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2, True)
        py = np.copy(fluxcor)
        py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
        kepplot.plot1d(px, py, cadence, 'r', 2, '#ffff00', 0.0, True)
        plt.setp(plt.gca(), xticklabels=[])
        kepplot.labels('', ylab, 'k', 14)
        plt.grid()
        # plot aperture photometry and PCA corrected data
        ax = kepplot.location([0.06, 0.09, 0.93, 0.43])
        yr = np.array([], 'float32')
        npc = min([6, nrem])
        for i in range(npc - 1, -1, -1):
            py = pcar[:, i] * c[i]
            py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
            cl = float(i) / (float(npc))
            kepplot.plot1d(px, py, cadence, [1.0 - cl, 0.0, cl], 2, '#ffff00',
                           0.0, True)
            yr = np.append(yr, py)
        y1 = max(yr)
        y2 = -min(yr)
        kepplot.RangeOfPlot(px, np.array([-y1, y1, -y2, y2]), 0.01, False)
        kepplot.labels(xlab, 'Principal Components', 'k', 14)
        plt.grid()
        # save plot to file
        plt.savefig(repname)
        # render plot
        plt.show()

    # stop time
    kepmsg.clock('KEPPCA ended at', logfile, verbose)

def reject_outliers(data, m):
    """Outlier rejection for computing robust mean"""
    try:
        return data[np.abs(data - np.nanmean(data)) < m * np.std(data)]
    except:
        print("Warning: Could not reject outliers.")
        return data

def mad(data):
    """
    Mean absolute deviation function used for fitting the PCA components to
    the data and subtracting them out
    """
    return np.nanmean(np.absolute(data - np.nanmean(data)))

def keppca_main():
    import argparse
    parser = argparse.ArgumentParser(
       description='Pixel-level principal component analysis of time series',
       formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input target pixel FITS file',
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-keppca.'),
                        default=None)
    parser.add_argument('--maskfile', help='Name of mask defintion ASCII file',
                        default='ALL', type=str)
    parser.add_argument('--components', default='1-3',
                        help='Principal components to be removed', type=str)
    parser.add_argument('--plotpca', action='store_true',
                        help='Create PCA plots?')
    parser.add_argument('--nmaps', default=10,
                        help='Number of principal components to include in report',
                        type=int)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='keppca.log', dest='logfile', type=str)
    args = parser.parse_args()
    keppca(args.infile, args.outfile, args.maskfile, args.components,
           args.plotpca, args.nmaps, args.overwrite, args.verbose, args.logfile)
