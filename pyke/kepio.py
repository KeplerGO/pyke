"""
This module contains utility functions for i/o operations.
"""
from . import kepmsg, kepkey
import numpy as np
import os
import glob
import tempfile
import shutil
from astropy.io import fits as pyfits


__all__ = ['delete', 'overwrite', 'openascii', 'closeascii', 'splitfits',
           'readfitstab', 'readfitscol', 'readtimecol', 'readsapcol',
           'readsaperrcol', 'readpdccol', 'readpdcerrcol', 'readcbvcol',
           'readsapqualcol', 'readlctable', 'tabappend', 'readimage',
           'writeimage', 'writefits', 'tmpfile', 'symlink', 'fileexists',
           'move', 'copy', 'parselist', 'createdir', 'createtree',
           'timeranges', 'cadence', 'timekeys', 'filterNaN', 'readTPF',
           'readMaskDefinition', 'readPRFimage']


def delete(filename, logfile, verbose):
    try:
        os.remove(filename)
    except:
        message = 'ERROR -- KEPIO.OBLITERATE: could not delete ' + filename
        kepmsg.err(logfile, message, verbose)

def overwrite(filename, logfile, verbose):
    if (os.path.isfile(filename)):
        try:
            delete(filename, logfile, verbose)
        except:
            message = 'ERROR -- KEPIO.CLOBBER: could not overwrite ' + filename
            kepmsg.err(logfile, message, verbose)

def openascii(filename, mode, logfile, verbose):
    try:
        content = open(filename, mode)
    except:
        message = ('ERROR -- KEPIO.OPENASCII: cannot open ASCII file ' +
                    filename)
        kepmsg.err(logfile, message, verbose)
    return content

def closeascii(file_,logfile,verbose):
    try:
        file_.close()
    except:
        message = ('ERROR - KEPIO.CLOSEASCII: cannot close ASCII file ' +
                   str(file_))
        kepmsg.err(logfile, message, verbose)

def splitfits(fitsfile, logfile, verbose):
    fitsfile = fitsfile.strip()
    if '+' in fitsfile:
        component = fitsfile.split('+')
        filename = str(component[0])
        hdu = int(component[1])
    elif '[' in fitsfile:
        fitsfile = fitsfile.strip(']')
        component = fitsfile.split('[')
        filename = str(component[0])
        hdu = int(component[1])
    else:
        errmsg = ('ERROR -- KEPIO.SPLITFITS: cannot determine HDU number '
                   'from name' + fitsfile)
        kepmsg.err(logfile, errmsg, verbose)
    return filename, hdu

def readfitstab(filename, hdu, logfile, verbose):
    try:
        table = hdu.data
    except:
        message = ('ERROR -- KEPIO.READFITSTAB: could not extract table '
                   'from ' + filename)
        kepmsg.err(logfile, message, verbose)
    return table

def readfitscol(filename, table, column, logfile, verbose):
    try:
        data = table.field(column)
    except:
        message  = ('ERROR -- KEPIO.READFITSCOL: could not extract '
                    + column + 'data from ' + filename)
        kepmsg.err(logfile, message, verbose)
    return data

def readtimecol(filename, table, logfile, verbose):
    try:
        data = table.field('TIME')
    except:
        try:
            data = table.field('barytime')
            if data[0] < 2.4e6 and data[0] > 1.0e4: data += 2.4e6
        except:
            message  = ('ERROR -- KEPIO.READTIMECOL: could not extract '
                        + 'time data from ' + filename)
            kepmsg.err(logfile, message, verbose)
    return data

def readsapcol(filename,table,logfile,verbose):
    try:
        data = table.field('SAP_FLUX')
    except:
        try:
            data = table.field('ap_raw_flux')
        except:
            message  = ('ERROR -- KEPIO.READSAPCOL: could not extract SAP flux'
                        'time series data from ' + filename)
            kepmsg.err(logfile, message, verbose)
    return data

def readsaperrcol(filename, table, logfile, verbose):
    """read FITS SAP error column"""
    try:
        data = table.field('SAP_FLUX_ERR')
    except:
        try:
            data = table.field('ap_raw_err')
        except:
            message  = ('ERROR -- KEPIO.READSAPERRCOL: could not extract SAP '
                        'flux error time series data from ' + filename)
            kepmsg.err(logfile, message, verbose)
    return data

def readpdccol(filename, table, logfile, verbose):
    """read FITS PDC column"""
    try:
        data = table.field('PDCSAP_FLUX')
    except:
        try:
            data = table.field('ap_corr_flux')
        except:
            message  = ('ERROR -- KEPIO.READPDCCOL: could not extract PDCSAP '
                        'flux time series data from ' + filename)
            kepmsg.err(logfile, message, verbose)
    return data

def readpdcerrcol(filename, table, logfile, verbose):
    """read FITS PDC error column"""
    try:
        data = table.field('PDCSAP_FLUX_ERR')
    except:
        try:
            data = table.field('ap_corr_err')
        except:
            message  = ('ERROR -- KEPIO.READPDCERRCOL: could not extract PDC '
                        'flux error time series data from ' + filename)
            kepmsg.err(logfile, message, verbose)
    return data

def readcbvcol(filename, table, logfile, verbose):
    """read FITS CBV column"""
    try:
        data = table.field('CBVSAP_FLUX')
    except:
        message  = ('ERROR -- KEPIO.READCBVCOL: could not extract CBVSAP flux '
                    'time series data from ' + filename)
        kepmsg.err(logfile, message, verbose)
    return data

def readsapqualcol(filename, table, logfile, verbose):
    """read quality column"""
    try:
        data = table.field('SAP_QUALITY')
    except:
        message  = ('ERROR -- KEPIO.READSAPQUALCOL: could not extract SAP '
                    'quality time series data from ' + filename)
        kepmsg.err(logfile, message, verbose)
    return data


def readlctable(infile,instr,logfile,verbose):
    """read all columns within Kepler FITS light curve table"""

    table = instr.data
    barytime        = readfitscol(infile, table, 'barytime', logfile, verbose)
    timcorr         = readfitscol(infile, table, 'timcorr', logfile, verbose)
    cadence_number  = readfitscol(infile, table, 'cadence_number', logfile,
                                  verbose)
    ap_cent_row     = readfitscol(infile, table, 'ap_cent_row', logfile,
                                  verbose)
    ap_cent_r_err   = readfitscol(infile, table, 'ap_cent_r_err', logfile,
                                  verbose)
    ap_cent_col     = readfitscol(infile, table, 'ap_cent_col', logfile,
                                  verbose)
    ap_cent_c_err   = readfitscol(infile, table, 'ap_cent_c_err', logfile,
                                  verbose)
    ap_raw_flux     = readfitscol(infile, table, 'ap_raw_flux', logfile,
                                  verbose)
    ap_raw_err      = readfitscol(infile, table, 'ap_raw_err', logfile,
                                  verbose)
    ap_corr_flux    = readfitscol(infile, table, 'ap_corr_flux', logfile,
                                  verbose)
    ap_corr_err     = readfitscol(infile, table, 'ap_corr_err', logfile,
                                  verbose)
    ap_ins_flux     = readfitscol(infile, table, 'ap_ins_flux', logfile,
                                  verbose)
    ap_ins_err      = readfitscol(infile, table, 'ap_ins_err', logfile,
                                  verbose)
    dia_raw_flux    = readfitscol(infile, table, 'dia_raw_flux', logfile,
                                  verbose)
    dia_raw_err     = readfitscol(infile, table, 'dia_raw_err', logfile,
                                  verbose)
    dia_corr_flux   = readfitscol(infile, table, 'dia_corr_flux', logfile,
                                  verbose)
    dia_corr_err    = readfitscol(infile, table, 'dia_corr_err', logfile,
                                  verbose)
    dia_ins_flux    = readfitscol(infile, table, 'dia_ins_flux', logfile,
                                  verbose)
    dia_ins_err     = readfitscol(infile, table, 'dia_ins_err', logfile,
                                  verbose)

    return [barytime, timcorr, cadence_number, ap_cent_row, ap_cent_r_err,
            ap_cent_col, ap_cent_c_err, ap_raw_flux, ap_raw_err, ap_corr_flux,
            ap_corr_err, ap_ins_flux, ap_ins_err, dia_raw_flux, dia_raw_err,
            dia_corr_flux, dia_corr_err, dia_ins_flux, dia_ins_err]

def tabappend(hdu1, hdu2, logfile, verbose):
    """append two table HDUs"""
    nrows1 = hdu1.data.shape[0]
    nrows2 = hdu2.data.shape[0]
    nrows = nrows1 + nrows2
    out = pyfits.BinTableHDU.from_columns(hdu1.columns, nrows=nrows)
    for name in hdu1.columns.names:
        try:
            out.data.field(name)[nrows1:] = hdu2.data.field(name)
        except:
            errmsg  = ('WARNING -- KEPIO.TABAPPEND: could not append '
                       'column ' + str(name))
            kepmsg.warn(logfile, errmsg, verbose)

    return out

def readimage(image, hdu, logfile,verbose):
    """ read image from HDU structure"""
    try:
        imagedata = image[hdu].data
    except:
        errmsg = ('ERROR -- KEPIO.READIMAGE: cannot read image data from HDU '
                   + str(hdu))
        kepmsg.err(logfile, errmsg, verbose)
    return imagedata

def writeimage(image, hdu, imagedata, logfile, verbose):
    """write image to HDU structure"""
    try:
        image[hdu].data = imagedata
    except:
        errmsg = ('ERROR -- KEPIO.WRITEIMAGE: Cannot write image data to HDU '
                   + str(hdu))
        kepmsg.err(logfile, errmsg, verbose)
    return image

def writefits(hdu, filename, overwrite, logfile, verbose):
    """write new FITS file"""

    if os.path.isfile(filename) and overwrite:
        delete(filename, logfile, verbose)
    try:
        hdu.writeto(filename)
    except:
        errmsg = ('ERROR -- KEPIO.WRITEFITS: Cannot create FITS file '
                  + filename)
        kepmsg.err(logfile, errmsg, verbose)

def tmpfile(path, suffix, logfile, verbose):
    """create a temporary file name"""
    try:
        tempfile.tempdir = path
        tmpfile = tempfile.mktemp() + suffix
    except:
        message = ('ERROR -- KEPIO.TMPFILE: Cannot create temporary file name')
        kepmsg.err(logfile,message,verbose)
    return tmpfile

def symlink(infile,linkfile,overwrite,logfile,verbose):
    """create symbolic link"""

    if os.path.exists(linkfile) and not overwrite:
        errmsg = ('ERROR: KEPIO.SYMLINK -- file ' + linkfile + ' exists, use '
                  'overwrite')
        kepmsg.err(logfile, errmsg, verbose)
    if overwrite:
        try:
            os.remove(linkfile)
        except:
            pass
    try:
        os.symlink(infile,linkfile)
    except:
        errmsg  = ('ERROR: KEPIO.SYMLINK -- could not create symbolic link '
                    'from ' + infile + ' to ' + linkfile)
        kepmsg.err(logfile, message, verbose)

def fileexists(file_):
    """check that a file exists"""
    if not os.path.isfile(file_):
        return False
    return True

def move(file1, file2, logfile, verbose):
    """move file"""
    try:
        shutil.move(file1, file2)
        message = 'KEPIO.MOVE -- moved ' + file1 + ' to ' + file2
        kepmsg.log(logfile, message, verbose)
    except:
        errmsg = ('ERROR -- KEPIO.MOVE: Could not move ' + file1 + ' to '
                  + file2)
        kepmsg.err(logfile, errmsg, verbose)

def copy(file1, file2, logfile, verbose):
    """copy file"""
    try:
        shutil.copy2(file1, file2)
        message = 'KEPIO.COPY -- copied ' + file1 + ' to ' + file2
        kepmsg.log(logfile, message, verbose)
    except:
        errmsg = ('ERROR -- KEPIO.COPY: could not copy ' + file1 +
                  ' to ' + file2)
        kepmsg.err(logfile, errmsg, verbose)

def parselist(inlist, logfile, verbose):
    """reate a list from a file, string or wildcard"""

    inlist.strip()
    if len(inlist) == 0 or inlist.count(' ') > 0:
        errmsg = 'ERROR -- KEPIO.PARSELIST: list not specified'
        kepmsg.err(logfile, errmsg, verbose)

    if inlist[0] == '@':
        infile = inlist.lstrip('@')
        if not os.path.isfile(infile):
            errmsg = ('ERROR -- KEPIO.PARSELIST: input list ' + infile +
                      'doest not exist')
            kepmsg.err(logfile,message,verbose)

    outlist = []
    if inlist[0] == '@':
        line = ' '
        infile = open(inlist.lstrip('@'))
        while line:
            line = infile.readline()
            if len(line.strip()) > 0:
                outlist.append(line.rstrip('\r\n'))
    elif inlist[0] != '@' and inlist.count('*') == 0:
        if inlist.count(',') == 0:
            outlist.append(inlist)
        else:
            list_ = inlist.split(',')
            for listitem in list_:
                outlist.append(listitem)
    elif inlist[0] != '@' and inlist.count('*') > 0:
        outlist = glob.glob(inlist)
    if len(outlist) == 0:
        errmsg = 'ERROR -- KEPIO.PARSELIST: raw input image list is empty'
        kepmsg.err(logfile, errmsg, verbose)

    return outlist


def createdir(path, logfile, verbose):
    """create a directory"""
    path = path.strip()
    if path[-1] != '/':
        path += '/'
    if not os.path.exists(path):
        try:
            os.mkdir(path)
            message = 'KEPIO.CREATEDIR -- Created directory ' + path
            kepmsg.log(logfile, message, verbose)
        except:
            errmsg = ('ERROR -- KEPIO.CREATEDIR: Could not create directory '
                      + path)
            kepmsg.err(logfile, message, verbose)
    else:
        message = 'KEPIO.CREATEDIR -- ' + path + ' directory exists'
        kepmsg.log(logfile, message, verbose)

def createtree(path,logfile,verbose):
    """create a directory tree"""
    path = path.strip()
    if path[-1] != '/':
        path += '/'
    if not os.path.exists(path):
        try:
            os.makedirs(path)
            message = 'KEPIO.CREATETREE -- Created directory tree ' + path
            kepmsg.log(logfile,message,verbose)
        except:
            errmsg = ('ERROR -- KEPIO.CREATETREE: Could not create directory '
                      'tree ' + path)
            kepmsg.err(logfile, errmsg, verbose)
    else:
        message = 'KEPIO.CREATETREE -- ' + path + ' directory exists'
        kepmsg.log(logfile, message, verbose)

def timeranges(ranges, logfile, verbose):
    """read time ranges from ascii file"""

    tstart = []
    tstop = []
    try:
        ranges = ranges.strip().split(';')
        for i in range(len(ranges)):
            tstart.append(float(ranges[i].strip().split(',')[0]))
            tstop.append(float(ranges[i].strip().split(',')[1]))
            if tstart[-1] == 0.0 and tstop[-1] == 0.0: tstop[-1] = 1.0e8
    except:
        tstart = []
        tstop = []
    if len(tstart) == 0 or len(tstop) == 0 or len(tstart) != len(tstop):
        errmsg = ('ERROR -- KEPIO.TIMERANGES: cannot understand time '
                  'ranges provided')
        kepmsg.err(logfile, errmsg, verbose)

    return tstart, tstop

def cadence(instr, infile, logfile, verbose):
    """manual calculation of median cadence within a time series"""
    try:
        intime = instr[1].data.field('barytime')
    except:
        intime = readfitscol(infile, instr[1].data, 'time', logfile, verbose)
    dt = []
    for i in range(1, len(intime)):
        if np.isfinite(intime[i]) and np.isfinite(intime[i-1]):
            dt.append(intime[i] - intime[i-1])
    dt = np.array(dt, dtype='float32')
    cadnce = np.median(dt) * 86400.0

    return intime[0], intime[-1], len(intime), cadnce

def timekeys(instr, filename, logfile, verbose):
    """read time keywords"""
    tstart = 0.0
    tstop = 0.0
    cadence = 0.0

    # BJDREFI
    try:
        bjdrefi = instr[1].header['BJDREFI']
    except:
        bjdrefi = 0.0

    # BJDREFF
    try:
        bjdreff = instr[1].header['BJDREFF']
    except:
        bjdreff = 0.0
    bjdref = bjdrefi + bjdreff

    # TSTART
    try:
        tstart = instr[1].header['TSTART']
    except:
        try:
            tstart = instr[1].header['STARTBJD'] + 2.4e6
        except:
            try:
                tstart = instr[0].header['LC_START'] + 2400000.5
            except:
                try:
                    tstart = instr[1].header['LC_START'] + 2400000.5
                except:
                    errmsg = ('ERROR -- KEPIO.TIMEKEYS: Cannot find TSTART, '
                              'STARTBJD or LC_START in ' + filename)
                    kepmsg.err(logfile, errmsg, verbose)
    tstart += bjdref

    # TSTOP
    try:
        tstop = instr[1].header['TSTOP']
    except:
        try:
            tstop = instr[1].header['ENDBJD'] + 2.4e6
        except:
            try:
                tstop = instr[0].header['LC_END'] + 2400000.5
            except:
                try:
                    tstop = instr[1].header['LC_END'] + 2400000.5
                except:
                    errmsg = ('ERROR -- KEPIO.TIMEKEYS: Cannot find TSTOP, '
                              'STOPBJD or LC_STOP in ' + filename)
                    kepmsg.err(logfile, errmsg, verbose)
    tstop += bjdref

    # OBSMODE
    cadence = 1.0
    try:
        obsmode = instr[0].header['OBSMODE']
    except:
        try:
            obsmode = instr[1].header['DATATYPE']
        except:
            errmsg = ('ERROR -- KEPIO.TIMEKEYS: cannot find keyword OBSMODE '
                      'or DATATYPE in ' + filename)
            kepmsg.err(logfile, errmsg, verbose)
    if 'short' in obsmode:
        cadence = 54.1782
    elif 'long' in obsmode:
        cadence = 1625.35

    return tstart, tstop, bjdref, cadence

def filterNaN(instr, datacol, outfile, logfile, verbose):
    """filter input data table"""
    try:
        nanclean = instr[1].header['NANCLEAN']
    except:
        naxis2 = 0
        for i in range(len(instr[1].columns.names)):
            if 'time' in instr[1].columns.names[i].lower():
                timecol = instr[1].columns.names[i]
        try:
            instr[1].data.field(datacol)
        except:
            msg = ("ERROR -- KEPIO.FILTERNAN: cannot find column {}"
                   "in the infile".format(datacol))
            kepmsg.err(logfile, msg, verbose)
        try:
            for i in range(len(instr[1].data.field(0))):
                if (str(instr[1].data.field(timecol)[i]) != '-inf' and
                    str(instr[1].data.field(datacol)[i]) != '-inf'):
                    instr[1].data[naxis2] = instr[1].data[i]
                    naxis2 += 1
            instr[1].data = instr[1].data[:naxis2]
            comment = 'NaN cadences removed from data'
            kepkey.new('NANCLEAN', True, comment, instr[1], outfile, logfile,
                       verbose)
        except:
            errmsg = ('ERROR -- KEPIO.FILTERNAN: Failed to filter NaNs from '
                      + outfile)
            kepmsg.err(logfile, errmsg, verbose)
    return instr

def readTPF(infile, colname, logfile, verbose):
    """ Read a Target Pixel File (TPF).

    Parameters
    ----------
    infile : str
        target pixel file name.
    colname : str
        name of the column to be read.
    logfile : str
    verbose : bool
    """

    try:
        tpf = pyfits.open(infile, mode='readonly', memmap=True)
    except:
        errmsg = ('ERROR -- KEPIO.OPENFITS: cannot open ' +
                  infile + ' as a FITS file')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        naxis2 = tpf['TARGETTABLES'].header['NAXIS2']
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No NAXIS2 keyword in ' + infile +
                  '[TARGETTABLES]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        kepid = tpf[0].header['KEPLERID']
        kepid = str(kepid)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No KEPLERID keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        channel = tpf[0].header['CHANNEL']
        channel = str(channel)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No CHANNEL keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        skygroup = tpf[0].header['SKYGROUP']
        skygroup = str(skygroup)
    except:
        skygroup = '0'
    try:
        module = tpf[0].header['MODULE']
        module = str(module)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No MODULE keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        output = tpf[0].header['OUTPUT']
        output = str(output)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No OUTPUT keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        quarter = tpf[0].header['QUARTER']
        quarter = str(quarter)
    except:
        try:
            quarter = tpf[0].header['CAMPAIGN']
            quarter = str(quarter)
        except:
            errmsg = ('ERROR -- KEPIO.READTPF: No QUARTER or CAMPAIGN ' +
                      'keyword in ' + infile + '[0]')
            kepmsg.err(logfile, errmsg, verbose)
    try:
        season = tpf[0].header['SEASON']
        season = str(season)
    except:
        season = '0'
    try:
        ra = tpf[0].header['RA_OBJ']
        ra = str(ra)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No RA_OBJ keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        dec = tpf[0].header['DEC_OBJ']
        dec = str(dec)
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: No DEC_OBJ keyword in ' + infile +
                  '[0]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        kepmag = tpf[0].header['KEPMAG']
        kepmag = str(float(kepmag))
    except:
        kepmag = ''
    try:
        tdim5 = tpf['TARGETTABLES'].header['TDIM5']
        xdim = int(tdim5.strip().strip('(').strip(')').split(',')[0])
        ydim = int(tdim5.strip().strip('(').strip(')').split(',')[1])
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: Cannot read TDIM5 keyword in ' +
                  infile + '[TARGETTABLES]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        crv5p1 = tpf['TARGETTABLES'].header['1CRV5P']
        column = crv5p1
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: Cannot read 1CRV5P keyword in ' +
                  infile + '[TARGETTABLES]')
        kepmsg.err(logfile, errmsg, verbose)
    try:
        crv5p2 = tpf['TARGETTABLES'].header['2CRV5P']
        row = crv5p2
    except:
        errmsg = ('ERROR -- KEPIO.READTPF: Cannot read 2CRV5P keyword in ' +
                  infile + '[TARGETTABLES]')
        kepmsg.err(logfile, errmsg, verbose)
    # read and close TPF data pixel image
    try:
        pixels = tpf['TARGETTABLES'].data.field(colname)[:]
    except:
        errmsg = ("\nERROR -- KEPIO.READTPF: Cannot read {0} "
                  "column in {1} '[TARGETTABLES]'".format(colname, infile))
        kepmsg.err(logfile, errmsg, verbose)

    tpf.close()

    # for STSCI_PYTHON v2.12 - convert 3D data array to 2D

    if len(np.shape(pixels)) == 3:
        isize = np.shape(pixels)[0]
        jsize = np.shape(pixels)[1]
        ksize = np.shape(pixels)[2]
        pixels = np.reshape(pixels, (isize, jsize * ksize))

    return (kepid, channel, skygroup, module, output, quarter, season,
            ra, dec, column, row, kepmag, xdim, ydim, pixels)


def readMaskDefinition(infile, logfile, verbose):
    """read target pixel mask data"""

    # open input file
    inf = pyfits.open(infile, 'readonly')

    # read bitmap image
    try:
        img = inf['APERTURE'].data
    except:
        txt = ('WARNING -- KEPIO.READMASKDEFINITION: Cannot read mask '
               'defintion in ' + infile + '[APERTURE]')
        kepwarn.err(txt, logfile)
    try:
        naxis1 = inf['APERTURE'].header['NAXIS1']
    except:
        txt = ('WARNING -- KEPIO.READMASKDEFINITION: Cannot read NAXIS1 '
               'keyword in ' + infile + '[APERTURE]')
        kepwarn.err(txt, logfile)
    try:
        naxis2 = inf['APERTURE'].header['NAXIS2']
    except:
        txt = ('WARNING -- KEPIO.READMASKDEFINITION: Cannot read NAXIS2 '
               'keyword in ' + infile + '[APERTURE]')
        kepwarn.err(txt, logfile)

    # read WCS keywords
    crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p = kepkey.getWCSp(
            infile, inf['APERTURE'], logfile, verbose)

    pixelcoord1 = np.zeros((naxis1, naxis2))
    pixelcoord2 = np.zeros((naxis1, naxis2))
    for j in range(naxis2):
        for i in range(naxis1):
            pixelcoord1[i, j] = kepkey.wcs(i, crpix1p, crval1p, cdelt1p)
            pixelcoord2[i, j] = kepkey.wcs(j, crpix2p, crval2p, cdelt2p)
    # close input file
    inf.close()
    return img, pixelcoord1, pixelcoord2


def readPRFimage(infile, hdu, logfile, verbose):
    """read pixel response file"""

    prf = pyfits.open(infile, 'readonly')

    # read bitmap image
    try:
        img = prf[hdu].data
    except:
        txt = ('ERROR -- KEPIO.READPRFIMAGE: Cannot read PRF image in '
               + infile + '[' + str(hdu) + ']')
        kepmsg.err(logfile, txt, verbose)
    try:
        naxis1 = prf[hdu].header['NAXIS1']
    except:
        txt = ('ERROR -- KEPIO.READPRFIMAGE: Cannot read NAXIS1 keyword in '
               + infile + '[' + str(hdu) + ']')
        kepmsg.err(logfile, txt, verbose)
    try:
        naxis2 = prf[hdu].header['NAXIS2']
    except:
        txt = ('ERROR -- KEPIO.READPRFIMAGE: Cannot read NAXIS2 keyword in '
               + infile + '[' + str(hdu) + ']')
        kepmsg.err(logfile, txt, verbose)

    # read WCS keywords
    crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p = kepkey.getWCSp(
            infile, prf[hdu], logfile, verbose
                                                                         )

    prf.close()
    return img, crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p
