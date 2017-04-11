#!/usr/bin/env python

import kepmsg, kepkey
import sys, tempfile, os, shutil, glob, numpy, warnings
from astropy.io import fits as pyfits

# -----------------------------------------------------------
# delete a file

def delete(file,logfile,verbose):

    status = 0
    try:
	os.remove(file)
    except:
	message = 'ERROR -- KEPIO.DELETE: could not delete ' + file
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# clobber a file

def clobber(file,logfile,verbose):

    status = 0
    if (os.path.isfile(file)):
	try:
	    status = delete(file,logfile,verbose)
	except:
	    message = 'ERROR -- KEPIO.CLOBBER: could not clobber ' + file
	    status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# open ASCII file

def openascii(file,type,logfile,verbose):

    status = 0
    try:
	content = open(file,type)
    except:
	message = 'ERROR -- KEPIO.OPENASCII: cannot open ASCII file ' + file
	status = kepmsg.err(logfile,message,verbose)

    return content, status

# -----------------------------------------------------------
# close ASCII file

def closeascii(file,logfile,verbose):

    status = 0
    try:
	file.close()
    except:
	message = 'ERROR - KEPIO.CLOSEASCII: cannot close ASCII file ' + str(file)
	status = kepmsg.err(logfile,message,verbose)
    
    return status

# -----------------------------------------------------------
# split FITS filename and HDU number

def splitfits(file,logfile,verbose):

    status = 0
    file = file.strip()
    if ('+' in file):
	component = file.split('+')
	filename = str(component[0])
	hdu = int(component[1])
    elif ('[' in file):
	file = file.strip(']')
	component = file.split('[')
	filename = str(component[0])
	hdu = int(component[1])
    else:
	message = 'ERROR -- KEPIO.SPLITFITS: cannot determine HDU number from name' + file
	status = kepmsg.err(logfile,message,verbose)
    return filename, hdu, status

# -----------------------------------------------------------
# open HDU structure

def openfits(file,mode,logfile,verbose):

    status = 0
    try:
        struct = pyfits.open(file,mode=mode)
    except:
	message = 'ERROR -- KEPIO.OPENFITS: cannot open ' + file + ' as a FITS file'
	status = kepmsg.err(logfile,message,verbose)
	struct = None
    return struct, status

# -----------------------------------------------------------
# close HDU structure

def closefits(struct,logfile,verbose):

    status = 0
    try:
	struct.close()
    except:
	message = 'ERROR -- KEPIO.CLOSEFITS: cannot close HDU structure'
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# read FITS table HDU

def readfitstab(file,hdu,logfile,verbose):

    status = 0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            table = hdu.data
    except:
	message = 'ERROR -- KEPIO.READFITSTAB: could not extract table from ' + file
	status = kepmsg.err(logfile,message,verbose)
	table = None
    return table, status

# -----------------------------------------------------------
# read FITS table column

def readfitscol(file,table,column,logfile,verbose):

    status = 0
    try:
	data = table.field(column)
    except:
	message  = 'ERROR -- KEPIO.READFITSCOL: could not extract ' + column 
	message += ' data from ' + file
	status = kepmsg.err(logfile,message,verbose)
	data = None
    return data, status

# -----------------------------------------------------------
# read TIME column

def readtimecol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('TIME')
    except:
        try:
            data = table.field('barytime')
            if data[0] < 2.4e6 and data[0] > 1.0e4: data += 2.4e6
        except:
            message  = 'ERROR -- KEPIO.READTIMECOL: could not extract'
            message += ' time data from ' + file
            status = kepmsg.err(logfile,message,verbose)
            data = None
    return data, status

# -----------------------------------------------------------
# read FITS SAP column

def readsapcol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('SAP_FLUX')
    except:
        try:
            data = table.field('ap_raw_flux')
        except:
            message  = 'ERROR -- KEPIO.READSAPCOL: could not extract SAP flux'
            message += ' time series data from ' + file
            status = kepmsg.err(logfile,message,verbose)
            data = None
    return data, status

# -----------------------------------------------------------
# read FITS SAP error column

def readsaperrcol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('SAP_FLUX_ERR')
    except:
        try:
            data = table.field('ap_raw_err')
        except:
            message  = 'ERROR -- KEPIO.READSAPERRCOL: could not extract SAP flux error'
            message += ' time series data from ' + file
            status = kepmsg.err(logfile,message,verbose)
            data = None
    return data, status

# -----------------------------------------------------------
# read FITS PDC column

def readpdccol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('PDCSAP_FLUX')
    except:
        try:
            data = table.field('ap_corr_flux')
        except:
            message  = 'ERROR -- KEPIO.READPDCCOL: could not extract PDCSAP flux'
            message += ' time series data from ' + file
            status = kepmsg.err(logfile,message,verbose)
            data = None
    return data, status

# -----------------------------------------------------------
# read FITS PDC error column

def readpdcerrcol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('PDCSAP_FLUX_ERR')
    except:
        try:
            data = table.field('ap_corr_err')
        except:
            message  = 'ERROR -- KEPIO.READPDCERRCOL: could not extract PDC flux error'
            message += ' time series data from ' + file
            status = kepmsg.err(logfile,message,verbose)
            data = None
    return data, status

# -----------------------------------------------------------
# read FITS CBV column

def readcbvcol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('CBVSAP_FLUX')
    except:
        message  = 'ERROR -- KEPIO.READCBVCOL: could not extract CBVSAP flux'
        message += ' time series data from ' + file
        status = kepmsg.err(logfile,message,verbose)
        data = None
    return data, status

# -----------------------------------------------------------
# read quality column

def readsapqualcol(file,table,logfile,verbose):

    status = 0
    try:
        data = table.field('SAP_QUALITY')
    except:
        message  = 'ERROR -- KEPIO.READSAPQUALCOL: could not extract SAP quality'
        message += ' time series data from ' + file
        status = kepmsg.err(logfile,message,verbose)
        data, status = None
    return data, status

# -----------------------------------------------------------
# read all columns within Kepler FITS light curve table

def readlctable(infile,instr,logfile,verbose):

    status = 0
    table = instr.data
    barytime, status = readfitscol(infile,table,'barytime',logfile,verbose)
    timcorr, status = readfitscol(infile,table,'timcorr',logfile,verbose)
    cadence_number, status = readfitscol(infile,table,'cadence_number',logfile,verbose)
    ap_cent_row, status = readfitscol(infile,table,'ap_cent_row',logfile,verbose)
    ap_cent_r_err, status = readfitscol(infile,table,'ap_cent_r_err',logfile,verbose)
    ap_cent_col, status = readfitscol(infile,table,'ap_cent_col',logfile,verbose)
    ap_cent_c_err, status = readfitscol(infile,table,'ap_cent_c_err',logfile,verbose)
    ap_raw_flux, status = readfitscol(infile,table,'ap_raw_flux',logfile,verbose)
    ap_raw_err, status = readfitscol(infile,table,'ap_raw_err',logfile,verbose)
    ap_corr_flux, status = readfitscol(infile,table,'ap_corr_flux',logfile,verbose)
    ap_corr_err, status = readfitscol(infile,table,'ap_corr_err',logfile,verbose)
    ap_ins_flux, status = readfitscol(infile,table,'ap_ins_flux',logfile,verbose)
    ap_ins_err, status = readfitscol(infile,table,'ap_ins_err',logfile,verbose)
    dia_raw_flux, status = readfitscol(infile,table,'dia_raw_flux',logfile,verbose)
    dia_raw_err, status = readfitscol(infile,table,'dia_raw_err',logfile,verbose)
    dia_corr_flux, status = readfitscol(infile,table,'dia_corr_flux',logfile,verbose)
    dia_corr_err, status = readfitscol(infile,table,'dia_corr_err',logfile,verbose)
    dia_ins_flux, status = readfitscol(infile,table,'dia_ins_flux',logfile,verbose)
    dia_ins_err, status = readfitscol(infile,table,'dia_ins_err',logfile,verbose)

    return [barytime, timcorr, cadence_number, ap_cent_row, ap_cent_r_err, \
                ap_cent_col, ap_cent_c_err, ap_raw_flux, ap_raw_err, ap_corr_flux, \
                ap_corr_err, ap_ins_flux, ap_ins_err, dia_raw_flux, dia_raw_err, \
                dia_corr_flux, dia_corr_err, dia_ins_flux, dia_ins_err], status

# -----------------------------------------------------------
# append two table HDUs

def tabappend(hdu1,hdu2,logfile,verbose):

    status = 0
    nrows1 = hdu1.data.shape[0]
    nrows2 = hdu2.data.shape[0]
    nrows = nrows1 + nrows2
    out = pyfits.new_table(hdu1.columns,nrows=nrows)
    for name in hdu1.columns.names:
	try:
	    out.data.field(name)[nrows1:] = hdu2.data.field(name)
	except:
	    message  = 'WARNING -- KEPIO.TABAPPEND: could not append column '
	    message += str(name)
	    status = kepmsg.warn(logfile,message,verbose)

    return out, status

# -----------------------------------------------------------
# read image from HDU structure

def readimage(struct,hdu,logfile,verbose):

    status = 0
    try:
	imagedata = struct[hdu].data
    except:
	message = 'ERROR -- KEPIO.READIMAGE: cannot read image data from HDU ' + str(hdu)
	status = kepmsg.err(logfile,message,verbose)
    return imagedata, status

# -----------------------------------------------------------
# write image to HDU structure

def writeimage(struct,hdu,imagedata,logfile,verbose):

    status = 0
    try:
	struct[hdu].data = imagedata
    except:
	message = 'ERROR -- KEPIO.WRITEIMAGE: Cannot write image data to HDU ' + str(hdu)
	status = kepmsg.err(logfile,message,verbose)
    return struct, status

# -----------------------------------------------------------
# write new FITS file

def writefits(hdu,filename,clobber,logfile,verbose):

    status = 0
    if (os.path.isfile(filename) and clobber):
	delete(filename,logfile,verbose)
    try:
	hdu.writeto(filename)
    except:
	message = 'ERROR -- KEPIO.WRITEFITS: Cannot create FITS file ' + filename
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# create a temporary file name

def tmpfile(path,suffix,logfile,verbose):

    status = 0
    try:
	tempfile.tempdir = path
        file = tempfile.mktemp() + suffix
    except:
	message = 'ERROR -- KEPIO.TMPFILE: Cannot create temporary file name'
	status = kepmsg.err(logfile,message,verbose)
    return file, status

# -----------------------------------------------------------
# create symbolic link

def symlink(infile,linkfile,clobber,logfile,verbose):

# delete file if one of the same name already exists

    status = 0
    if (os.path.exists(linkfile) and not clobber):
	message = 'ERROR: KEPIO.SYMLINK -- file ' + linkfile + ' exists, use clobber'
	status = kepmsg.err(logfile,message,verbose)
    if (status == 0 and clobber):
	try:
	    os.remove(linkfile)
	except:
	    status = 0

# create symbolic link

    if (status == 0):
	try:
	    os.symlink(infile,linkfile)
	except:
	    message  = 'ERROR: KEPIO.SYMLINK -- could not create symbolic link from '
	    message += infile + ' to ' + linkfile
	    status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# check that a file exists

def fileexists(file):

    status = True
    if not os.path.isfile(file):
	status = False
    return status

# -----------------------------------------------------------
# move file
   
def move(file1,file2,logfile,verbose):

    status = 0
    message = 'KEPIO.MOVE -- moved ' + file1 + ' to ' + file2
    try:
        shutil.move(file1,file2)
        kepmsg.log(logfile,message,verbose)
    except:
	message = 'ERROR -- KEPIO.MOVE: Could not move ' + file1 + ' to ' + file2
	status = kepmsg.err(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# copy file
   
def copy(file1,file2,logfile,verbose):

    status = 0
    message = 'KEPIO.COPY -- copied ' + file1 + ' to ' + file2
    try:
	shutil.copy2(file1,file2)
	kepmsg.log(logfile,message,verbose)
    except:
	message = 'ERROR -- KEPIO.COPY: could not copy ' + file1 + ' to ' + file2
	status = kepmsg.err(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# create a list from a file, string or wildcard

def parselist(inlist,logfile,verbose):

# test input name list
   
    status = 0
    inlist.strip()
    if (len(inlist) == 0 or inlist.count(' ') > 0):
	message = 'ERROR -- KEPIO.PARSELIST: list not specified'
	status = kepmsg.err(logfile,message,verbose)

# test @filelist exists

    if (inlist[0] == '@'):
	infile = inlist.lstrip('@')
	if not os.path.isfile(infile):
	    message = 'ERROR -- KEPIO.PARSELIST: input list '+infile+' does not exist'
	    status = kepmsg.err(logfile,message,verbose)

# parse wildcard and comma-separated lists

    outlist = []
    if (status == 0 and inlist[0] == '@'):
	line = ' '
        infile = open(inlist.lstrip('@'))
	while line:
	    line = infile.readline()
	    if (len(line.strip()) > 0):
		outlist.append(line.rstrip('\r\n'))
    elif (status == 0 and inlist[0] != '@' and inlist.count('*') == 0):
	if (inlist.count(',') == 0):
	    outlist.append(inlist)
        else:
	    list = inlist.split(',')
	    for listitem in list:
		outlist.append(listitem)
    elif (status == 0 and inlist[0] != '@' and inlist.count('*') > 0):
	outlist = glob.glob(inlist)
    if (status == 0 and len(outlist) == 0):
	message = 'ERROR -- KEPIO.PARSELIST: raw input image list is empty'
	status = kepmsg.err(logfile,message,verbose)

    return outlist, status

# -----------------------------------------------------------
# create a directory

def createdir(path,logfile,verbose):

    status = 0
    path = path.strip()
    message = 'KEPIO.CREATEDIR -- Created directory ' + path
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
	try:
            os.mkdir(path)
	    kepmsg.log(logfile,message,verbose)
	except:
	    message  = 'ERROR -- KEPIO.CREATEDIR: Could not create '
	    message += 'directory ' + path
	    status = kepmsg.err(logfile,message,verbose)
    else:
	message = 'KEPIO.CREATEDIR -- ' + path + ' directory exists'
	kepmsg.log(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# create a directory tree

def createtree(path,logfile,verbose):

    status = 0
    path = path.strip()
    message = 'KEPIO.CREATETREE -- Created directory tree ' + path
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
	try:
            os.makedirs(path)
	    kepmsg.log(logfile,message,verbose)
	except:
	    message  = 'ERROR -- KEPIO.CREATETREE: Could not create '
	    message += 'directory tree ' + path
	    status = kepmsg.err(logfile,message,verbose)
    else:
	message = 'KEPIO.CREATETREE -- ' + path + ' directory exists'
	kepmsg.log(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# number of HDU within a FITS structure

def HDUnum(struct):

    ValidHDU = True
    nhdu = 0
    while ValidHDU:
        try:
            struct[nhdu].header[0]
            nhdu += 1
        except:
            ValidHDU = False

    return nhdu

# -----------------------------------------------------------
# read time ranges from ascii file

def timeranges(ranges,logfile,verbose):

    status = 0; tstart = []; tstop = []
    if '@' in ranges:
        try:
            lines, status = openascii(ranges[1:],'r',logfile,verbose)
        except:
            txt = 'ERROR -- KEPIO.TIMERANGES: cannot open file ' + ranges[1:]
            status = kepmsg.err(logfile,txt,verbose)
            return tstart, tstop, status
        for line in lines:
            line = line.strip().split(',')
            if len(line) == 2:
                try:
                    float(line[0])
                    float(line[1])
                    tstart.append(float(line[0]))
                    tstop.append(float(line[1]))
                    if tstart[-1] == 0.0 and tstop[-1] == 0.0: tstop[-1] = 1.0e8
                except:
                    continue
        status = closeascii(lines,logfile,verbose)
        if len(tstart) == 0 or len(tstop) == 0 or len(tstart) != len(tstop) or status > 0:
            txt = 'ERROR -- KEPIO.TIMERANGES: cannot understand content of ' + ranges[1:]
            status = kepmsg.err(logfile,txt,verbose)
            return tstart, tstop, status
    else:
        try:
            ranges = ranges.strip().split(';')
            for i in range(len(ranges)):
                tstart.append(float(ranges[i].strip().split(',')[0]))
                tstop.append(float(ranges[i].strip().split(',')[1]))
                if tstart[-1] == 0.0 and tstop[-1] == 0.0: tstop[-1] = 1.0e8
        except:
            tstart = []; tstop = []
        if len(tstart) == 0 or len(tstop) == 0 or len(tstart) != len(tstop) or status > 0:
            txt = 'ERROR -- KEPIO.TIMERANGES: cannot understand time ranges provided'
            status = kepmsg.err(logfile,txt,verbose)
            return tstart, tstop, status
        
    return tstart, tstop, status

## -----------------------------------------------------------
## manual calculation of median cadence within a time series

def cadence(instr,infile,logfile,verbose,status):
    
    try:
        intime = instr[1].data.field('barytime')
    except:
        intime, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
    dt = []
    for i in range(1,len(intime)):
        if numpy.isfinite(intime[i]) and numpy.isfinite(intime[i-1]):
            dt.append(intime[i] - intime[i-1])
    dt = numpy.array(dt,dtype='float32')
    cadnce = numpy.median(dt) * 86400.0

    return intime[0], intime[-1], len(intime), cadnce, status

# -----------------------------------------------------------
# read time keywords

def timekeys(instr,file,logfile,verbose,status):

    tstart = 0.0; tstop = 0.0; cadence = 0.0

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
            tstart = instr[1].header['STARTBJD']
            tstart += 2.4e6
	except:
            try:
                tstart = instr[0].header['LC_START']
                tstart += 2400000.5
            except:
                try:
                    tstart = instr[1].header['LC_START']
                    tstart += 2400000.5
                except:
                    message  = 'ERROR -- KEPIO.TIMEKEYS: Cannot find TSTART, STARTBJD or '
                    message += 'LC_START in ' + file
                    status = kepmsg.err(logfile,message,verbose)
    tstart += bjdref

# TSTOP

    try:
        tstop = instr[1].header['TSTOP']
    except:
	try:
            tstop = instr[1].header['ENDBJD']
            tstop += 2.4e6
	except:
            try:
                tstop = instr[0].header['LC_END']
                tstop += 2400000.5
            except:
                try:
                    tstop = instr[1].header['LC_END']
                    tstop += 2400000.5
                except:
                    message  = 'ERROR -- KEPIO.TIMEKEYS: Cannot find TSTOP, STOPBJD or '
                    message += 'LC_STOP in ' + file 
                    status = kepmsg.err(logfile,message,verbose)
    tstop += bjdref

# OBSMODE

    cadence = 1.0
    try:
        obsmode = instr[0].header['OBSMODE']
    except:
        try:
            obsmode = instr[1].header['DATATYPE']
        except:
            message  = 'ERROR -- KEPIO.TIMEKEYS: cannot find keyword OBSMODE '
            message += 'or DATATYPE in ' + file
            status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        if 'short' in obsmode: # and bjdref == 0.0:
            cadence = 54.1782
        elif 'long' in obsmode: # and bjdref == 0.0:
            cadence = 1625.35

    return tstart, tstop, bjdref, cadence, status


# -----------------------------------------------------------
# filter input data table

def filterNaN(instr,datacol,outfile,logfile,verbose):

    status = 0
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
            msg = 'ERROR -- KEPIO.FILTERNAN: cannot find column ' + datacol + ' in the infile'
            status = kepmsg.err(logfile,msg,verbose)

        if status == 0:
            try:
                for i in range(len(instr[1].data.field(0))):
                    if str(instr[1].data.field(timecol)[i]) != '-inf' and \
                            str(instr[1].data.field(datacol)[i]) != '-inf':
                        instr[1].data[naxis2] = instr[1].data[i]
                        naxis2 += 1
                instr[1].data = instr[1].data[:naxis2]
                comment = 'NaN cadences removed from data'
                status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)
            except:
                msg = 'ERROR -- KEPIO.FILTERNAN: Failed to filter NaNs from '+ outfile
                status = kepmsg.err(logfile,msg,verbose)

    return instr, status


# -----------------------------------------------------------
# read target pixel data file

def readTPF(infile,colname,logfile,verbose):

    status = 0
    tpf = pyfits.open(infile,mode='readonly',memmap=True)
    if status == 0:
        try:
            naxis2 = tpf['TARGETTABLES'].header['NAXIS2']
        except:
            txt = 'ERROR -- KEPIO.READTPF: No NAXIS2 keyword in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            kepid = tpf[0].header['KEPLERID']
            kepid = str(kepid)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No KEPLERID keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            channel = tpf[0].header['CHANNEL']
            channel = str(channel)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No CHANNEL keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            skygroup = tpf[0].header['SKYGROUP']
            skygroup = str(skygroup)
        except:
            skygroup = '0'
    if status == 0:
        try:
            module = tpf[0].header['MODULE']
            module = str(module)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No MODULE keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            output = tpf[0].header['OUTPUT']
            output = str(output)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No OUTPUT keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            quarter = tpf[0].header['QUARTER']
            quarter = str(quarter)
        except:
            try:
                quarter = tpf[0].header['CAMPAIGN']
                quarter = str(quarter)
            except:
                txt = 'ERROR -- KEPIO.READTPF: No QUARTER or CAMPAIGN keyword in ' + infile + '[0]'
                status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            season = tpf[0].header['SEASON']
            season = str(season)
        except:
            season = '0'
    if status == 0:
        try:
            ra = tpf[0].header['RA_OBJ']
            ra = str(ra)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No RA_OBJ keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            dec = tpf[0].header['DEC_OBJ']
            dec = str(dec)
        except:
            txt = 'ERROR -- KEPIO.READTPF: No DEC_OBJ keyword in ' + infile + '[0]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            kepmag = tpf[0].header['KEPMAG']
            kepmag = str(float(kepmag))
        except:
            kepmag = ''
    if status == 0:
        try:
            tdim5 = tpf['TARGETTABLES'].header['TDIM5']
            xdim = int(tdim5.strip().strip('(').strip(')').split(',')[0])
            ydim = int(tdim5.strip().strip('(').strip(')').split(',')[1])
        except:
            txt = 'ERROR -- KEPIO.READTPF: Cannot read TDIM5 keyword in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            crv5p1 = tpf['TARGETTABLES'].header['1CRV5P']
            column = crv5p1
        except:
            txt = 'ERROR -- KEPIO.READTPF: Cannot read 1CRV5P keyword in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            crv5p2 = tpf['TARGETTABLES'].header['2CRV5P']
            row = crv5p2
        except:
            txt = 'ERROR -- KEPIO.READTPF: Cannot read 2CRV5P keyword in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,verbose)

# read and close TPF data pixel image

    if status == 0:
        try:
            pixels = tpf['TARGETTABLES'].data.field(colname)[:]
        except:
            pixels = None
            txt = '\nWARNING -- KEPIO.READTPF: Cannot read ' + colname + ' column in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        status = closefits(tpf,logfile,verbose)

# for STSCI_PYTHON v2.12 - convert 3D data array to 2D

    if status == 0 and len(numpy.shape(pixels)) == 3:
        isize = numpy.shape(pixels)[0]
        jsize = numpy.shape(pixels)[1]
        ksize = numpy.shape(pixels)[2]
        pixels = numpy.reshape(pixels,(isize,jsize*ksize))
    
    return kepid, channel, skygroup, module, output, quarter, season, \
        ra, dec, column, row, kepmag, xdim, ydim, pixels, status

# -----------------------------------------------------------
# read target pixel mask data

def readMaskDefinition(infile,logfile,verbose):

    status = 0

# open input file

    inf, status = openfits(infile,'readonly',logfile,verbose)

# read bitmap image

    if status == 0:
        try:
            img = inf['APERTURE'].data
        except:
            txt = 'WARNING -- KEPIO.READMASKDEFINITION: Cannot read mask defintion in ' + infile + '[APERTURE]'
            kepwarn.err(txt,logfile)
            status = 1
    if status == 0:
        try:
            naxis1 = inf['APERTURE'].header['NAXIS1']
        except:
            txt = 'WARNING -- KEPIO.READMASKDEFINITION: Cannot read NAXIS1 keyword in ' + infile + '[APERTURE]'
            kepwarn.err(txt,logfile)
            status = 1
    if status == 0:
        try:
            naxis2 = inf['APERTURE'].header['NAXIS2']
        except:
            txt = 'WARNING -- KEPIO.READMASKDEFINITION: Cannot read NAXIS2 keyword in ' + infile + '[APERTURE]'
            kepwarn.err(txt,logfile)
            status = 1

# read WCS keywords

    if status == 0:
        crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status = \
            kepkey.getWCSp(infile,inf['APERTURE'],logfile,verbose)     
    if status == 0:
        pixelcoord1 = numpy.zeros((naxis1,naxis2))
        pixelcoord2 = numpy.zeros((naxis1,naxis2))
        for j in range(naxis2):
            for i in range(naxis1):
                pixelcoord1[i,j] = kepkey.wcs(i,crpix1p,crval1p,cdelt1p)
                pixelcoord2[i,j] = kepkey.wcs(j,crpix2p,crval2p,cdelt2p)

# close input file

    if status == 0:
        status = closefits(inf,logfile,verbose)

    return img, pixelcoord1, pixelcoord2, status

# -----------------------------------------------------------
# read pixel response file

def readPRFimage(infile,hdu,logfile,verbose):

    status = 0

# open input file

    prf, status = openfits(infile,'readonly',logfile,verbose)

# read bitmap image

    if status == 0:
        try:
            img = prf[hdu].data
        except:
            txt = 'ERROR -- KEPIO.READPRFIMAGE: Cannot read PRF image in ' + infile + '[' + str(hdu) + ']'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            naxis1 = prf[hdu].header['NAXIS1']
        except:
            txt = 'ERROR -- KEPIO.READPRFIMAGE: Cannot read NAXIS1 keyword in ' + infile + '[' + str(hdu) + ']'
            status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            naxis2 = prf[hdu].header['NAXIS2']
        except:
            txt = 'ERROR -- KEPIO.READPRFIMAGE: Cannot read NAXIS2 keyword in ' + infile + '[' + str(hdu) + ']'
            status = kepmsg.err(logfile,txt,verbose)

# read WCS keywords

    if status == 0:
        crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status = \
            kepkey.getWCSp(infile,prf[hdu],logfile,verbose)     

# close input file

    if status == 0:
        status = closefits(prf,logfile,verbose)

    return img, crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status
