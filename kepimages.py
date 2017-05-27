import sys
import numpy as np
from astropy.io import fits as pyfits
from . import kepio
from . import kepmsg
from . import kepkey
from . import kepstat

__all__ = ['kepimages']

def kepimages(infile, prefix, imtype='FLUX', ranges='0,0', clobber=True,
              verbose=True, logfile='kepimages.log'):
    """
    kepimages -- create a series of separate FITS image files from a Target
    Pixel File

    infile : str
        Filename for the input Target Pixel File.
    prefix : str
        Prefix for the filenames of the output FITS images. Individual
        filenames will be prefix_BJDddddddd.dddd.fits, where
        ddddddd.dddd is the mid-time of the exposure in units of BJD.
    ranges : str
        The user can choose here specific time ranges of exposures from which
        to export images. Time ranges are supplied as comma-separated pairs of
        Barycentric Julian Dates (BJDs). Multiple ranges are separated by a
        semi-colon. An example containing two time ranges is:
            ``'2455641.658,2455641.740;2455671.658,2455672.740'``
    clobber : bool
        Overwrite the output file? if clobber is False and an existing file has
        the same name as outfile then the task will stop with an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.
    """

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPIMAGES -- '
            'infile={}'.format(infile)
            'prefix={}'.format(outfix)
            'imtype={}'.format(imtype)
            'ranges={}'.format(ranges)
            'clobber={}'.format(clobber)
            'verbose={}'.format(verbose)
            'logfile={}'.format(logfile))
    kepmsg.log(logfile,call+'\n',verbose)


    kepmsg.clock('KEPIMAGES started at',logfile,verbose)

    # open input file
    print (' ')
    instr = pyfits.open(infile, mode='readonly', memmap=True)
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr,file,logfile,verbose)

    # ingest time series data
    time = instr[1].data.field('TIME')[:] + 2454833.0
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

    # choose output image
    if imtype.lower() == 'raw_cnts':
        outim = raw_cnts
    elif imtype.lower() == 'flux_err':
        outim = flux_err
    elif imtype.lower() == 'flux_bkg':
        outim = flux_bkg
    elif imtype.lower() == 'flux_bkg_err':
        outim = flux_bkg_err
    elif imtype.lower() == 'cosmic_rays':
        outim = cosmic_rays
    else:
        outim = flux

    # identify images to be exported
    tim = np.array([])
    dat = np.array([])
    err = np.array([])
    tstart, tstop = kepio.timeranges(ranges, logfile, verbose)
    cadencelis = kepstat.filterOnRange(time, tstart, tstop)

    # provide name for each output file and clobber if file exists
    for cadence in cadencelis:
        outfile = prefix + '_BJD%.4f' % time[cadence] + '.fits'
        if clobber:
            kepio.clobber(outfile, logfile, verbose)
        if kepio.fileexists(outfile):
            errmsg = ('ERROR -- KEPIMAGES: {} exists. Use --clobber'
                      .format(clobber))
            kepmsg.err(logfile, errmsg, True)

    # construct output primary extension
    ncad = 0
    for cadence in cadencelis:
        outfile = prefix + '_BJD%.4f' % time[cadence] + '.fits'
        hdu0 = pyfits.PrimaryHDU()
        for i in range(len(cards0)):
            try:
                if cards0[i].keyword not in hdu0.header.keys():
                    hdu0.header[cards0[i].keyword] = (cards0[i].value,
                                                      cards0[i].comment)
                else:
                    hdu0.header.cards[cards0[i].key].comment = cards0[i].comment
            except:
                pass
        kepkey.history(call, hdu0, outfile, logfile, verbose)
        outstr = pyfits.HDUList(hdu0)

        # construct output image extension
        hdu1 = pyfits.ImageHDU(flux[cadence])
        for i in range(len(cards2)):
            try:
                if cards2[i].keyword not in hdu1.header.keys():
                    hdu1.header[cards2[i].key] = (cards2[i].value,
                                                  cards2[i].comment)
            except:
                pass
        for i in range(len(cards1)):
            if (cards1[i].keyword not in hdu1.header.keys() and
                cards1[i].keyword[:4] not in ['TTYP', 'TFOR', 'TUNI',
                                              'TDIS', 'TDIM', 'WCAX',
                                              '1CTY', '2CTY', '1CRP',
                                              '2CRP', '1CRV', '2CRV',
                                              '1CUN', '2CUN', '1CDE',
                                              '2CDE', '1CTY', '2CTY',
                                              '1CDL', '2CDL', '11PC',
                                              '12PC', '21PC', '22PC',
                                              'WCSN','TFIE']):
                hdu1.header[cards1[i].keyword] = (cards1[i].value,
                                                  cards1[i].comment)
        try:
            int_time = cards1['INT_TIME'].value
        except:
            kepmsg.warn(logfile,
                        'WARNING -- KEPIMAGES: cannot find INT_TIME keyword')
        try:
            frametim = cards1['FRAMETIM'].value
        except:
            kepmsg.warn(logfile,
                        'WARNING -- KEPIMAGES: cannot find FRAMETIM keyword')
        try:
            num_frm = cards1['NUM_FRM'].value
        except:
            kepmsg.warn(logfile,
                        'WARNING -- KEPIMAGES: cannot find NUM_FRM keyword')

        hdu1.header['EXTNAME'] = ('IMAGE','name of extension')

        try:
            hdu1.header['TELAPSE'] = (frametim * num_frm,
                                      '[s] elapsed time for exposure')
        except:
            hdu1.header['TELAPSE'] = (-999, '[s] elapsed time for exposure')
        try:
            hdu1.header['LIVETIME'] = (int_time * num_frm,
                                       '[s] TELASPE multiplied by DEADC')
        except:
            hdu1.header['LIVETIME'] = (-999, '[s] TELASPE multiplied by DEADC')
        try:
            hdu1.header['EXPOSURE'] = (int_time * num_frm,
                                       '[s] time on source')
        except:
            hdu1.header['EXPOSURE'] = (-999, '[s] time on source')
        try:
            hdu1.header['MIDTIME'] = (time[cadence],
                                      '[BJD] mid-time of exposure')
        except:
            hdu1.header['MIDTIME'] = (-999, '[BJD] mid-time of exposure')
        try:
            hdu1.header['TIMECORR'] = (timecorr[cadence],
                                       '[d] barycenter - timeslice correction')
        except:
            hdu1.header['TIMECORR'] = (-999,
                                       '[d] barycenter - timeslice correction')
        try:
            hdu1.header['CADENCEN'] = (cadenceno[cadence],
                                       'unique cadence number')
        except:
            hdu1.header['CADENCEN'] = (-999, 'unique cadence number')
        try:
            hdu1.header['QUALITY'] = (quality[cadence], 'pixel quality flag')
        except:
            hdu1.header['QUALITY'] = (-999, 'pixel quality flag')
        try:
            if True in np.isfinite(cosmic_rays[cadence]):
                hdu1.header['COSM_RAY'] = (True, 'cosmic ray detected?')
            else:
                hdu1.header['COSM_RAY'] = (False, 'cosmic ray detected?')
        except:
            hdu1.header['COSM_RAY'] = (-999, 'cosmic ray detected?')
        try:
            pc1 = str(pos_corr1[cadence])
            pc2 = str(pos_corr2[cadence])
            hdu1.header['POSCORR1'] = (pc1, '[pix] column position correction')
            hdu1.header['POSCORR2'] = (pc2, '[pix] row position correction')
        except:
            hdu1.header['POSCORR1'] = (-999, '[pix] column position correction')
            hdu1.header['POSCORR2'] = (-999, '[pix] row position correction')
        outstr.append(hdu1)

        # write output file
        outstr.writeto(outfile, checksum=True)
        ncad += 1
        txt  = '\r%3d%% ' % (float(ncad) / float(len(cadencelis)) * 100.0)
        txt += '%s ' % outfile
        sys.stdout.write(txt)
        sys.stdout.flush()

    # close input structure
    kepio.closefits(instr, logfile, verbose)
    print ('\n')

    # end time
    kepmsg.clock('KEPIMAGES finished at',logfile,verbose)


def kepimages_main():
    import argparse
    parser = argparse.ArgumentParser(
            description=('Export images within a Target Pixel File to a '
                         'series of FITS image files'))
    parser.add_argument('infile', help='Name of input target pixel file',
                        type=str)
    parser.add_argument('prefix',
                        help='Prefix of name for each output image file',
                        type=str)
    parser.add_argument('--imtype', default='FLUX', dest='imtype',
                        help='Image type', type=str,
                        choices=['RAW_CNTS', 'FLUX', 'FLUX_ERR', 'FLUX_BKG',
                                 'FLUX_BKG_ERR','COSMIC_RAYS'])
    parser.add_argument('--ranges', help='Time ranges for output [BJD]',
                        default='0,0', type=str)
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='kepimages.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepimages(args.infile, args.prefix, args.imtype, args.ranges,
              args.clobber, args.verbose, args.logfile)
