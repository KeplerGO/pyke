"""
Name: keparith.py
Written by: Tom Barclay
Date released: Nov 10 2011

Changelog:
1.0 released
2.0 made it possible to run from the command line if there are arguements given to the call to keparith
2.1 fixed a bug when dividing a number and using PDCSAP_DATA - thanks to Joseph Long for spotting this
"""

import sys
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import kepio, kepmsg, kepkey, kepfit, kepstat
from numpy import median, subtract, maximum, ones, multiply, float32, shape, absolute, mean, std, isfinite, where, nan

def MAD(xx,minSd=1E-16):
    """Median Absolute Deviation
    """
    med = median(xx, 0)
    absdev = np.absolute(subtract(xx, med))
    mad = median(absdev, 0)
    mad = maximum(mad, np.multiply(ones(mad.shape, np.float32),
                       (minSd / 1.48)))
    return mad

def kepaddconstant(infile, outfile, datacol, constant, constantval, sign,
                   clobber, verbose, logfile, status):
    # log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPARITH -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'constant='+str(constant)+' '
    call += 'constantval='+str(constantval)+' '
    call += 'sign'+str(sign)+' '
    overwrite = 'n'
    if clobber: overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if verbose: chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

    # start time
    kepmsg.clock('KEPARITH started at',logfile,verbose)

    # test log file
    logfile = kepmsg.test(logfile)

    # clobber output file
    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPARITH: ' + outfile + ' exists. Use clobber=yes'
        status = kepmsg.err(logfile, message, verbose)

    ## open input file
    instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        try:
            test = str(instr[0].header['FILEVER'])
            version = 2
        except KeyError:
            version = 1
        lc_flux = instr[1].data.field(datacol)
        if datacol == 'SAP_FLUX':
            errcol = 'SAP_FLUX_ERR'
            try:
                lc_err = instr[1].data.field(errcol)
                haveerr = True
            except:
                haveerr = False
        elif datacol == 'PDCSAP_FLUX':
            errcol = 'PDCSAP_FLUX_ERR'
            try:
                lc_err = instr[1].data.field(errcol)
                haveerr = True
            except:
                haveerr = False
        else:
            errcol = datacol + '_ERR'
            try:
                lc_err = instr[1].data.field(errcol)
                haveerr = True
            except:
                try:
                    errcol = datacol + '_err'
                    lc_err = instr[1].data.field(errcol)
                    haveerr = True
                except:
                    haveerr = False
        #subtractor he just refers to the number that will be added/subtracted
        #divided or multiplied
        if isinstance(constantval,(long,int,float)) and constant == 'None':
            subtractor = float(constantval)
        elif constant.lower() == 'median':
            subtractor = float(median(lc_flux[isfinite(lc_flux)]))
        elif constant.lower() == 'mean':
            subtractor = float(mean(lc_flux[isfinite(lc_flux)]))
        elif constant.lower() == 'std':
            subtractor = float(std(lc_flux[isfinite(lc_flux)]))
        elif constant.lower() == 'mad':
            subtractor = float(MAD(lc_flux[isfinite(lc_flux)]))
        elif constant.lower() == 'max':
            subtractor = float(max(lc_flux[isfinite(lc_flux)]))
        elif constant.lower() == 'range':
            subtractor = float(max(lc_flux[isfinite(lc_flux)]) - min(lc_flux[isfinite(lc_flux)]))
        elif str(constant).lower() == 'none':
            subtractor = 0.
            message = 'No operation will be performed if you select None for the function and do not give a constant value'
            status = kepmsg.err(logfile,message,verbose)
        else:
            message = 'Your constant term is not in the list of possible functions'
            status = kepmsg.err(logfile,message,verbose)
    if subtractor == 0. and sign == 'divide' and status == 0:
        message = 'You are trying to divide by zero: not a good idea.'
        status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        if sign.lower() == 'add':
            instr[1].data.field(datacol)[:] = where(isfinite(instr[1].data.field(datacol)[:]),
                                                    (lc_flux + subtractor),nan)
        elif sign.lower() == 'subtract':
            instr[1].data.field(datacol)[:] = where(isfinite(instr[1].data.field(datacol)[:]),
                                                    (lc_flux - subtractor),nan)
        elif sign.lower() == 'divide':
            instr[1].data.field(datacol)[:] = where(isfinite(instr[1].data.field(datacol)[:]),
                                                    (lc_flux / subtractor),nan)
            if haveerr:
                instr[1].data.field(errcol)[:] = where(isfinite(instr[1].data.field(errcol)[:]),
                                                       (lc_err / subtractor),nan)
        elif sign.lower() == 'multiply':
                instr[1].data.field(datacol)[:] = where(isfinite(instr[1].data.field(datacol)[:]),
                                                        (lc_flux * subtractor),nan)
            if haveerr:
                instr[1].data.field(errcol)[:] = where(isfinite(instr[1].data.field(errcol)[:]),
                                                       (lc_err * subtractor),nan)
        else:
            message = 'Your operation need to be one of: add, subtract, divide or multiply'
            status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        instr.writeto(outfile)
    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)
    ## end time
    if status == 0:
        message = 'KEPARITH completed at'
    else:
        message = '\nKEPARITH aborted at'
    kepmsg.clock(message,logfile,verbose)

# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Time invariant algebra on light curve data')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)
    parser.add_argument('--datacol', '-d', default='SAP_FLUX', dest='datacol',
                        help='Name of the column containing the flux time series',
                        type=str)
    parser.add_argument('--constantfunc', '-f', default='None', dest='constantfunc',
                        help='A value calculated from a function to be used in operation', type=str,
                        choices=['None','median','mean','MAD','std','max','range'])
    parser.add_argument('--constantval', '-v', default='None', dest='constantval',
                        help='Use a given number in operation', type=str)
    parser.add_argument('--operation', '-o', default='add', dest='operation',
                        help='Add, subtract, multiply, divide flux by a constant', type=str,
                        choices=['add','','subtract','multiply','divide'])
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    kepaddconstant(args.infile, args.outfile, args.datacol, args.constantfunc,
                   args.constantval, args.operation, args.clobber,
                   args.verbose, args.logfile, args.status)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$keparith.par")
    t = iraf.IrafTaskFactory(taskname="keparith", value=parfile, function=kepaddconstant)
