#!/usr/bin/env python

import kepmsg, kepio
import numpy
from astropy.io import fits as pyfits

# -----------------------------------------------------------
# get keyword value

def get(file,hdu,keyword,logfile,verbose):

    status = 0
    try:
	value = hdu.header[keyword]
    except:
	message = 'ERROR -- KEPKEY.GET: Cannot read keyword ' + keyword
	message += ' in file ' + file
	status = kepmsg.err(logfile,message,verbose)
	value = None
    return value, status

# -----------------------------------------------------------
# delete keyword

def delete(keyword,hdu,file,logfile,verbose):

    status = 0
    try:
	del hdu.header[keyword]
    except:
	message = 'ERROR -- KEPKEY.DELETE: Cannot delete keyword ' + keyword
	message += ' in ' + file
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# add new keyword

def new(keyword,value,comment,hdu,file,logfile,verbose):

    status = 0
    try:
        hdu.header.update(keyword,value,comment)
    except:
	message = 'ERROR -- KEPKEY.NEW: Cannot create keyword ' + keyword
	message += ' in ' + file
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# add comment keyword

def comment(txt,hdu,file,logfile,verbose):

    status = 0
    try:
	hdu.header.add_comment(txt)
    except:
	message = 'ERROR -- KEPKEY.COMMENT: Cannot create comment keyword'
	message += ' in ' + file
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# add history keyword

def history(txt,hdu,file,logfile,verbose):

    status = 0
    try:
	hdu.header.add_history(txt)
    except:
	message = 'ERROR -- KEPKEY.HISTORY: Cannot create history keyword'
	message += ' in ' + file
	status = kepmsg.err(logfile,message,verbose)
    return status

# -----------------------------------------------------------
# change existing keyword value

def change(keyword,value,hdu,file,logfile,verbose):

    status = 0
    try:
	hdu.header.update(keyword,value)
    except:
	message = 'ERROR -- KEPKEY.CHANGE: Cannot update keyword ' + keyword
	message += ' in ' + file
	status = kepmsg.err(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# calculate timestamp cadence from keywords

def cadence(struct,file,logfile,verbose):

# get keyword data

    status = 0
    try:
        int_time, status = get(file,struct,'INT_TIME',logfile,verbose)
    except:
        txt = 'ERROR -- KEPKEY.CADENCE: Cannot read keyword INT_TIME in file ' + file + '[1]'
	status = kepmsg.err(logfile,message,verbose)
    try:
        readtime, status = get(file,struct,'READTIME',logfile,verbose)
    except:
        txt = 'ERROR -- KEPKEY.CADENCE: Cannot read keyword READTIME in file ' + file + '[1]'
	status = kepmsg.err(logfile,message,verbose)
    try:
        num_frm, status = get(file,struct,'NUM_FRM',logfile,verbose)
    except:
        txt = 'ERROR -- KEPKEY.CADENCE: Cannot read keyword NUM_FRM in file ' + file + '[1]'
	status = kepmsg.err(logfile,message,verbose)

# calculate cadence

    cadence = (float(int_time) + float(readtime)) * float(num_frm)

    return cadence, status

# -----------------------------------------------------------
# get physical WCS keywords

def getWCSp(file,struct,logfile,verbose):

    status = 0
    crpix1p = 0.0
    crpix2p = 0.0
    crval1p = 0.0
    crval2p = 0.0
    cdelt1p = 0.0
    cdelt2p = 0.0
    try:
        crpix1p, status = get(file,struct,'CRPIX1P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRPIX1P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crpix2p, status = get(file,struct,'CRPIX2P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRPIX2P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crval1p, status = get(file,struct,'CRVAL1P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRVAL1P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crval2p, status = get(file,struct,'CRVAL2P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CRVAL2P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        cdelt1p, status = get(file,struct,'CDELT1P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CDELT1P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        cdelt2p, status = get(file,struct,'CDELT2P',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSP: Cannot read keyword CDELT2P in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1

    return crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status

# -----------------------------------------------------------
# get physical WCS keywords

def getWCSs(file,struct,logfile,verbose):

    status = 0
    crpix1 = 0.0
    crpix2 = 0.0
    crval1 = 0.0
    crval2 = 0.0
    cdelt1 = 0.0
    cdelt2 = 0.0
    pc = [0.0,0.0,0.0,0.0]
    try:
        crpix1, status = get(file,struct,'CRPIX1',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRPIX1 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crpix2, status = get(file,struct,'CRPIX2',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRPIX2 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crval1, status = get(file,struct,'CRVAL1',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRVAL1 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        crval2, status = get(file,struct,'CRVAL2',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CRVAL2 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        cdelt1, status = get(file,struct,'CDELT1',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CDELT1 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        cdelt2, status = get(file,struct,'CDELT2',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword CDELT2 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        pc1_1, status = get(file,struct,'PC1_1',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC1_1 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        pc1_2, status = get(file,struct,'PC1_2',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC1_2 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        pc2_1, status = get(file,struct,'PC2_1',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC2_1 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        pc2_2, status = get(file,struct,'PC2_2',logfile,verbose)
    except:
        txt = 'WARNING -- KEPKEY.GETWCSS: Cannot read keyword PC2_2 in file ' + file
        kepmsg.warn(logfile,txt)
        status = 1
    try:
        pc = numpy.array([[pc1_1,pc1_2],[pc2_1,pc2_2]])
        pc = numpy.linalg.inv(pc)
    except:
        pass
    

    return crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, pc, status

# -----------------------------------------------------------
# calculate coordinates from WCS data

def wcs(i,crpix,crval,cdelt):
    return crval + (float(i + 1) - crpix) * cdelt

# -----------------------------------------------------------
# remove empty keywords within a FITS file

def emptykeys(struct,file,logfile,verbose):

# determine number of HDU from keyword search

    nhdu = kepio.HDUnum(struct)

# delete empty keywords

    for hdu in range(nhdu):
        for keyword in struct[hdu].header.keys():
            head = struct[hdu].header[keyword]
            if ('pyfits' in str(head) and 'Undefined' in str(head)):
                status = delete(keyword,struct[hdu],file,logfile,verbose)
        
    return struct

