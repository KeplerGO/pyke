from .utils import PyKEArgumentHelpFormatter
import time, urllib
import sys
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from . import kepio
from . import kepmsg
from . import kepkey


__all__ = ['keptrim']


def keptrim(infile, column, row, imsize, outfile=None, kepid=None,
            overwrite=False, verbose=False, logfile='keptrim.log'):
    """
    keptrim -- trim pixels from Target Pixel Files

    keptrim will extract a square-shaped series of sub-images from a Target
    Pixel File. The simple purpose of this task is to reduce the size of large
    data sets such as the superstamps or two-wheel engineering data for the
    sake of processing efficiency. Performing a keptrim step speeds up
    calculations such as kepprfphot considertably and provides manual
    convenience for tasks such as kepmask.

    Parameters
    ----------
    infile : str
        Filename for the input Target Pixel File.
    column : int
        The CCD column number on which to center the output subimage.
    row : int
        The CCD row number on which to center the output subimage.
    imsize : int
        The pixel size of the subimage along either the row or column
        dimension. The subimage will be square.
    outfile : str
        Filename for the output Target Pixel File. This product will be written
        to the same FITS format as archived light curves.
    kepid : None or int
        If the target is catalogued within the Kepler Input Catalog (KIC), then
        the pixel row and column location will be extracted from the KIC
        provided the Kepler ID is provided. The user must be online for this
        feature to execute. If provided kepid will override column and row.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Option for verbose mode, in which informative messages and warnings to
        the shell and a logfile.
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ keptrim ktwo251248961-c112_lpd-targ.fits 14 770 --imsize 3
        --overwrite --verbose

    .. image:: ../_static/images/api/keptrim.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPTRIM -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' column={}'.format(column)
            + ' row={}'.format(row)
            + ' imsize={}'.format(imsize)
            + ' kepid={}'.format(kepid)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPTRIM started at', logfile, verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPTRIM: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile, mode='readonly', memmap=True)
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)

    # identify the season of observation
    try:
        season = cards0['SEASON'].value
    except:
        season = 0
    # retrieve column and row from KIC
    try:
        kic = FOVKepID(str(kepid))
        column = int(kic[98 + season * 5])
        row = int(kic[97 + season * 5])
    except:
        pass

    # convert CCD column and row to image column and row
    if imsize % 2 == 0:
        imsize += 1
    crpix1p = cards2['CRPIX1P'].value
    crpix2p = cards2['CRPIX2P'].value
    crval1p = cards2['CRVAL1P'].value
    crval2p = cards2['CRVAL2P'].value
    cdelt1p = cards2['CDELT1P'].value
    cdelt2p = cards2['CDELT2P'].value
    imcol = (column - crval1p) * cdelt1p + crpix1p - 1
    imrow = (row - crval2p) * cdelt2p + crpix2p - 1
    crval1p = column - imsize / 2 + 0.5
    crval2p = row - imsize / 2 + 0.5

    # check subimage is contained inside the input image
    naxis1 = cards2['NAXIS1'].value
    naxis2 = cards2['NAXIS2'].value
    x1 = int(imcol - imsize // 2 + 0.5)
    x2 = x1 + imsize
    y1 = int(imrow - imsize // 2 + 0.5)
    y2 = y1 + imsize

    if x1 < 0 or y1 < 0 or x2 > naxis1 or y2 > naxis2:
        errmsg = ('ERROR -- KEPTRIM: Requested pixel area falls outside of '
                  'the pixel image in file {}. Make the pixel area smaller '
                  'or relocate it''s center.'.format(infile))
        kepmsg.err(logfile, errmsg, verbose)

    # time series data
    time = instr[1].data.field('TIME')[:]
    timecorr = instr[1].data.field('TIMECORR')[:]
    cadenceno = instr[1].data.field('CADENCENO')[:]
    raw_cnts = instr[1].data.field('RAW_CNTS')[:]
    flux = instr[1].data.field('FLUX')[:]
    flux_err = instr[1].data.field('FLUX_ERR')[:]
    flux_bkg = instr[1].data.field('FLUX_BKG')[:]
    flux_bkg_err = instr[1].data.field('FLUX_BKG_ERR')[:]
    cosmic_rays = instr[1].data.field('COSMIC_RAYS')[:]
    quality = instr[1].data.field('QUALITY')[:]
    pos_corr1 = instr[1].data.field('POS_CORR1')[:]
    pos_corr2 = instr[1].data.field('POS_CORR2')[:]

    # resize time series
    raw_cnts = raw_cnts[:, y1:y2, x1:x2]
    flux = flux[:, y1:y2, x1:x2]
    flux_err = flux_err[:, y1:y2, x1:x2]
    flux_bkg = flux_bkg[:, y1:y2, x1:x2]
    flux_bkg_err = flux_bkg_err[:, y1:y2, x1:x2]
    cosmic_rays = cosmic_rays[:, y1:y2, x1:x2]

    # reshape time series images
    isize = np.shape(flux)[0]
    jsize = np.shape(flux)[1]
    ksize = np.shape(flux)[2]
    raw_cnts = np.reshape(raw_cnts, (isize, jsize * ksize))
    flux = np.reshape(flux, (isize, jsize * ksize))
    flux_err = np.reshape(flux_err, (isize, jsize * ksize))
    flux_bkg = np.reshape(flux_bkg, (isize, jsize * ksize))
    flux_bkg_err = np.reshape(flux_bkg_err, (isize, jsize * ksize))
    cosmic_rays = np.reshape(cosmic_rays, (isize, jsize * ksize))

    # pixel map data
    maskmap = np.array(instr[2].data[y1:y2,x1:x2])

    # construct output primary extension
    hdu0 = pyfits.PrimaryHDU()
    for i in range(len(cards0)):
        try:
            if cards0[i].keyword not in hdu0.header.keys():
                hdu0.header[cards0[i].keyword] = (cards0[i].value, cards0[i].comment)
            else:
                hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
        except:
            pass
    kepkey.history(call, hdu0, outfile, logfile, verbose)
    outstr = pyfits.HDUList(hdu0)

    # construct output light curve extension
    coldim = '(' + str(imsize) + ',' + str(imsize) + ')'
    eformat = str(imsize*imsize) + 'E'
    jformat = str(imsize*imsize) + 'J'
    kformat = str(imsize*imsize) + 'K'
    col1 =  pyfits.Column(name='TIME', format='D', unit='BJD - 2454833',
                          array=time)
    col2 =  pyfits.Column(name='TIMECORR', format='E', unit='d',
                          array=timecorr)
    col3 =  pyfits.Column(name='CADENCENO', format='J', array=cadenceno)
    col4 =  pyfits.Column(name='RAW_CNTS', format=jformat, unit='count',
                          dim=coldim, array=raw_cnts)
    col5 =  pyfits.Column(name='FLUX', format=eformat, unit='e-/s', dim=coldim,
                          array=flux)
    col6 =  pyfits.Column(name='FLUX_ERR', format=eformat, unit='e-/s',
                          dim=coldim,array=flux_err)
    col7 =  pyfits.Column(name='FLUX_BKG', format=eformat, unit='e-/s',
                          dim=coldim,array=flux_bkg)
    col8 =  pyfits.Column(name='FLUX_BKG_ERR', format=eformat, unit='e-/s',
                          dim=coldim, array=flux_bkg_err)
    col9 =  pyfits.Column(name='COSMIC_RAYS', format=eformat,unit='e-/s',
                          dim=coldim, array=cosmic_rays)
    col10 = pyfits.Column(name='QUALITY', format='J', array=quality)
    col11 = pyfits.Column(name='POS_CORR1', format='E', unit='pixel',
                          array=pos_corr1)
    col12 = pyfits.Column(name='POS_CORR2', format='E', unit='pixel',
                          array=pos_corr2)
    cols =  pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8,
                            col9, col10, col11, col12])
    hdu1 =  pyfits.BinTableHDU.from_columns(cols)
    for i in range(len(cards1)):
        try:
            if cards1[i].keyword not in hdu1.header.keys():
                hdu1.header[cards1[i].keyword] = (cards1[i].value,
                                                  cards1[i].comment)
            else:
                hdu1.header.cards[cards1[i].keyword].comment = cards1[i].comment
        except:
            pass
    hdu1.header['1CRV4P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV4P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX4'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX4'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    hdu1.header['1CRV5P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV5P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX5'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX5'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    hdu1.header['1CRV6P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV6P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX6'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX6'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    hdu1.header['1CRV7P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV7P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX7'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX7'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    hdu1.header['1CRV8P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV8P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX8'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX8'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    hdu1.header['1CRV9P'] = (crval1p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['2CRV9P'] = (crval2p,
                             '[pixel] detector coordinate at reference pixel')
    hdu1.header['1CRPX9'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 1')
    hdu1.header['2CRPX9'] = ((imsize + 1) / 2,
                             '[pixel] reference pixel along image axis 2')
    outstr.append(hdu1)

    # construct output mask bitmap extension
    hdu2 = pyfits.ImageHDU(maskmap)
    for i in range(len(cards2)):
        try:
            if cards2[i].keyword not in hdu2.header.keys():
                hdu2.header[cards2[i].keyword] = (cards2[i].value,
                                                  cards2[i].comment)
            else:
                hdu2.header.cards[cards2[i].keyword].comment = cards2[i].comment
        except:
            pass
    hdu2.header['NAXIS1' ] = (imsize, '')
    hdu2.header['NAXIS2' ] = (imsize, '')
    hdu2.header['CRVAL1P'] = (crval1p,
                              '[pixel] detector coordinate at reference pixel')
    hdu2.header['CRVAL2P'] = (crval2p,
                              '[pixel] detector coordinate at reference pixel')
    hdu2.header['CRPIX1' ] = ((imsize + 1) / 2,
                              '[pixel] reference pixel along image axis 1')
    hdu2.header['CRPIX2' ] = ((imsize + 1) / 2,
                              '[pixel] reference pixel along image axis 2')
    outstr.append(hdu2)

    # write output file
    print("Writing output file {}...".format(outfile))
    outstr.writeto(outfile,checksum=True)
    # close input structure
    instr.close()
    # end time
    kepmsg.clock('KEPTRIM finished at', logfile, verbose)

def FOVKepID(id):
    """KIC retrieval based upon KepID"""

    # build mast query
    url  = ('http://archive.stsci.edu/kepler/kepler_fov/search.php?'
            'action=Search&kic_kepler_id={}'.format(id) + '&max_records=100'
            '&verb=3&outputformat=CSV')

    # retrieve results from MAST
    out = ''
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0
            and 'Kepler' not in line
            and 'integer' not in line
            and 'no rows found' not in line):
            out = line.split(',')
    return out

def keptrim_main():
    import argparse
    parser = argparse.ArgumentParser(
             description='Trim unwanted pixels from a Target Pixel File',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input target pixel file',
                        type=str)
    parser.add_argument('column', help='CCD column number of the target',
                        type=int)
    parser.add_argument('row', help='CCD row number of the target', type=int)
    parser.add_argument('imsize',
                        help=('Number of pixels to extract in both row and'
                              ' column dimensions'), type=int)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-keptrim.'),
                        default=None)
    parser.add_argument('--kepid', type=int,
                        help='Kepler ID number from the Kepler Input Catalog')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='keptrim.log', type=str)
    args = parser.parse_args()
    keptrim(args.infile, args.column, args.row, args.imsize, args.outfile,
            args.kepid, args.overwrite, args.verbose, args.logfile)
