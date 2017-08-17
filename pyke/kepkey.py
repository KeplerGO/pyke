from . import kepmsg
import numpy as np
from astropy.io import fits as pyfits


__all__ = ['get', 'remove', 'new', 'comment', 'history', 'change', 'cadence',
           'getWCSp', 'getWCSs', 'wcs', 'emptykeys']


def get(filename, hdu, keyword, logfile, verbose):
    """get keyword value"""
    try:
        value = hdu.header[keyword]
    except:
        message = ('ERROR -- KEPKEY.GET: Cannot read keyword '
                   + keyword + ' in file ' + filename)
        kepmsg.err(logfile, message, verbose)
    return value


def remove(keyword, hdu, filename, logfile, verbose):
    """delete keyword"""
    try:
        del hdu.header[keyword]
    except:
        message = ('ERROR -- KEPKEY.REMOVE: Cannot delete keyword '
                   + keyword + ' in ' + filename)
        kepmsg.err(logfile, message, verbose)


def new(keyword, value, comment, hdu, filename, logfile, verbose):
    """add new keyword"""
    try:
        hdu.header[keyword] = (value, comment)
    except:
        message = ('ERROR -- KEPKEY.NEW: Cannot create keyword '
                   + keyword + ' in ' + filename)
        kepmsg.err(logfile, message, verbose)

def comment(txt, hdu, filename, logfile, verbose):
    """add comment keyword"""
    try:
        hdu.header.add_comment(txt)
    except:
        message = ('ERROR -- KEPKEY.COMMENT: Cannot create comment keyword'
                   + ' in ' + filename)
        kepmsg.err(logfile, message, verbose)

def history(txt, hdu, filename, logfile, verbose):
    """add history keyword"""
    try:
        hdu.header.add_history(txt)
    except:
        message = ('ERROR -- KEPKEY.HISTORY: Cannot create history keyword'
                    + ' in ' + filename)
        kepmsg.err(logfile, message, verbose)

def change(keyword, value, hdu, filename, logfile, verbose):
    """change existing keyword value"""
    try:
        hdu.header[keyword] = value
    except:
        message = ('ERROR -- KEPKEY.CHANGE: Cannot update keyword '
                   + keyword + ' in ' + filename)
        kepmsg.err(logfile, message, verbose)


def cadence(struct, filename, logfile, verbose):
    """calculate timestamp cadence from keywords"""
    try:
        int_time = get(filename, struct, 'INT_TIME', logfile, verbose)
    except:
        txt = ('ERROR -- KEPKEY.CADENCE: Cannot read keyword INT_TIME in file '
               + filename + '[1]')
        kepmsg.err(logfile, message, verbose)
    try:
        readtime = get(filename, struct, 'READTIME', logfile, verbose)
    except:
        txt = ('ERROR -- KEPKEY.CADENCE: Cannot read keyword READTIME in file '
               + filename + '[1]')
        kepmsg.err(logfile, message, verbose)
    try:
        num_frm = get(filename, struct, 'NUM_FRM', logfile, verbose)
    except:
        txt = ('ERROR -- KEPKEY.CADENCE: Cannot read keyword NUM_FRM in file '
                + filename + '[1]')
        kepmsg.err(logfile, message, verbose)

    cadence = (float(int_time) + float(readtime)) * float(num_frm)
    return cadence


def getWCSp(filename, struct, logfile, verbose):
    """get physical WCS keywords"""
    crpix1p = 0.0
    crpix2p = 0.0
    crval1p = 0.0
    crval2p = 0.0
    cdelt1p = 0.0
    cdelt2p = 0.0

    try:
        crpix1p = get(filename, struct, 'CRPIX1P', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRPIX1P in '
               'file ' + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crpix2p = get(filename, struct, 'CRPIX2P', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRPIX2P in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crval1p = get(filename, struct,'CRVAL1P', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRVAL1P in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crval2p = get(filename, struct,'CRVAL2P',logfile,verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRVAL2P in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        cdelt1p = get(filename, struct, 'CDELT1P', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CDELT1P in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        cdelt2p = get(filename, struct, 'CDELT2P', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSP: Cannot read keyword CDELT2P in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)

    return crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p


def getWCSs(filename, struct, logfile, verbose):
    """# get physical WCS keywords"""

    crpix1 = 0.0; crpix2 = 0.0; crval1 = 0.0; crval2 = 0.0; cdelt1 = 0.0;
    cdelt2 = 0.0
    pc = [0.0, 0.0, 0.0, 0.0]
    try:
        crpix1 = get(filename, struct, 'CRPIX1', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRPIX1 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crpix2 = get(filename, struct, 'CRPIX2', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRPIX2 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crval1 = get(filename, struct, 'CRVAL1', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRVAL1 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        crval2 = get(filename, struct, 'CRVAL2', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRVAL2 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        cdelt1 = get(filename, struct, 'CDELT1', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CDELT1 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        cdelt2 = get(filename, struct, 'CDELT2', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword CDELT2 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        pc1_1 = get(filename, struct, 'PC1_1', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC1_1 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        pc1_2 = get(filename, struct, 'PC1_2', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC1_2 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        pc2_1 = get(filename, struct, 'PC2_1', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC2_1 in file '
                + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        pc2_2 = get(filename, struct, 'PC2_2', logfile, verbose)
    except:
        txt = ('WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC2_2 in file '
               + filename)
        kepmsg.warn(logfile, txt, verbose)
    try:
        pc = np.array([[pc1_1, pc1_2], [pc2_1, pc2_2]])
        pc = np.linalg.inv(pc)
    except:
        pass

    return crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, pc


def wcs(i, crpix, crval, cdelt):
    """calculate coordinates from WCS data"""
    return crval + (float(i + 1) - crpix) * cdelt


def emptykeys(struct, filename, logfile, verbose):
    """remove empty keywords within a FITS file"""
    nhdu = HDUnum(struct)

    for hdu in range(nhdu):
        for keyword in struct[hdu].header.keys():
            head = struct[hdu].header[keyword]
            if 'pyfits' in str(head) and 'Undefined' in str(head):
                remove(keyword, struct[hdu], filename, logfile, verbose)
    return struct

def HDUnum(struct):
    """number of HDU within a FITS structure"""
    ValidHDU = True
    nhdu = 0
    while ValidHDU:
        try:
            struct[nhdu].header[0]
            nhdu += 1
        except:
            ValidHDU = False
    return nhdu

    """why not just len(struct.header[0])?"""
