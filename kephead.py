import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey
import sys, time, re

# -----------------------------------------------------------
# core code

def kephead(infile, outfile, keyname, clobber, verbose, logfile, status):

# input arguments

    status = 0
    np.seterr(all="ignore")

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPHEAD -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'keyname='+keyname+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time
    kepmsg.clock('KEPHEAD started at: ',logfile,verbose)

# test log file
    logfile = kepmsg.test(logfile)

# clobber output file
    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPHEAD: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# Is there an output file?
    if status == 0:
        if outfile.lower() == 'none':
            ofile = False
        else:
            ofile = outfile

# open input FITS file
    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)

# number of HDU in input file
    if status == 0:
        nhdu = kepio.HDUnum(instr)

# loop through each HDU in infile
    if status == 0:
        kepmsg.log(ofile,'',True)
        for hdu in range(nhdu):
            # how many keywords in the HDU?
            keylist = instr[hdu].header.cards
            nkeys = len(keylist)
            # print header number
            prhead = False
            for i in range(nkeys):
                if (keyname.upper() == 'ALL' or
                    keyname.upper() in instr[hdu].header.keys()[i]):
                    prhead = True
            if prhead:
                dashes = ''
                title = infile + '[' + str(hdu) + ']'
                for j in range(len(title)):
                    dashes += '-'
                kepmsg.log(ofile,dashes,True)
                kepmsg.log(ofile,title,True)
                kepmsg.log(ofile,dashes,True)
                kepmsg.log(ofile,'',True)

            # print keywords
            for i in range(nkeys):
                if ((keyname.upper() == 'ALL' or \
                         keyname.upper() in instr[hdu].header.keys()[i]) and \
                        'COMMENT' not in instr[hdu].header.keys()[i]):
                    kepmsg.log(ofile,str(keylist[i]),True)
            kepmsg.log(ofile,'',True)

    # stop time
    kepmsg.clock('KEPHEAD ended at: ',logfile,verbose)

    return
# -----------------------------------------------------------
# main
if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Search for and list FITS keywords in Kepler data files')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)
    parser.add_argument('keyname', help='Snippet of keyword name', type=str)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kephead.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    kephead(args.infile, args.outfile, args.keyname, args.clobber, args.verbose, args.logfile,args.status)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kephead.par")
    t = iraf.IrafTaskFactory(taskname="kephead", value=parfile, function=kephead)
