import numpy as np
import matplotlib.pyplot as plt
import pyfits
import kepio
import kepmsg
import sys
import datetime
import sys

def keptimefix(infile,outfile,clobber,verbose,logfile,status,cmdLine=False):

    """
    All Kepler light curve and target pixel files with version numbers 5.0 contain an error in
    the time stamps. This was fixed in the light curve with version 5.0 (at MAST after May 2013).
    The timescale for fixing the target pixel files is unclear but in the mean time this script will
    fix the target pixel file time stamps and make the times consistent with the light curve files.
    The error in Q0-13 can be corrected by adding 66.184s. During Q14 there was a leap second added
    Q15+ can be corrected by adding 67.184s. This tool fixes the time stamp accordingly.

    inputs
    infile - the name of the input target pixel file
    output - the name of the output target pixel file

    optional
    clobber (default=False) - overwrite a file with the stame name as outfile.
    """

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPTIMEFIX -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPTIMEFIX started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

    if clobber:
        status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPTIMEFIX: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

    instr, status = kepio.openfits(infile,'readonly',logfile,verbose)

    creator =  instr[0].header['CREATOR']
    if creator.find('TargetPixelExporterPipelineModule') < 0:
        message = 'ERROR -- KEPTIMEFIX: this file is not a target pixel file'
        status = kepmsg.err(logfile,message,verbose)

    if status == 0:
        header_ext1 = instr[1].header.cards
        data_ext1 = instr[1].data

        fileversion = instr[0].header['FILEVER']
        if float(fileversion) > 4.0:
            message = 'ERROR -- KEPTIMEFIX: no time fix needed for this file. FILEVER > 4.0'
            status = kepmsg.err(logfile,message,verbose)
            sys.exit(0)

        quarter = instr[0].header['QUARTER']
        if instr[0].header['OBSMODE'] == 'long cadence':
            cadencetype = 'L'
        elif instr[0].header['OBSMODE'] == 'short cadence':
            cadencetype = 'S'

        TIME_wrong = data_ext1.field('TIME')
        CADNUM = data_ext1.field('CADENCENO')
        TIMECORR_old = data_ext1.field('TIMECORR')

        ## update headers
        ##TSTART, TSTART, EXPOSURE, TELAPSE, LIVETIME
        ##DATE-OBS, DATE-END
        if cadencetype == 'L':
            offset = np.where(CADNUM <= 57139,66.184,67.184) / 86400.
        elif cadencetype == 'S':
            offset = np.where(CADNUM <= 1702663,66.184,67.184) / 86400.

        TIME_right = TIME_wrong + offset
        TIMECORR_new = TIMECORR_old + offset
        #tcol = pyfits.Column(name='TIME',format='D14.7',
        #    array=TIME_right, unit = 'BJD - 2454833', disp='D14.7')

        #instr[1].columns.change_attrib('TIME',array,TIME_right)

        #cols = instr[1].data.columns + tcol

        instr[1].data['TIME'][:] = TIME_right

        #we decided not to use the updated timecorr because
        #it is different from the LC FITS files by ~1 ms.
        instr[1].data['TIMECORR'][:] = np.nan * np.empty(len(TIMECORR_old))


        #instr[1] = pyfits.new_table(cols,header=instr[1].header)


        #now to fix the header
        tstart_right = instr[1].header['TSTART'] + offset[0]
        tstop_right = instr[1].header['TSTOP'] + offset[-1]

        telapse_right = tstop_right - tstart_right

        instr[1].header['TSTART'] = tstart_right
        instr[1].header['TSTOP'] = tstop_right
        instr[1].header['TELAPSE'] = telapse_right

        deadc = instr[1].header['DEADC']
        instr[1].header['LIVETIME'] = telapse_right * deadc

        #get the date-obs
        dstart = instr[1].header['DATE-OBS']
        dend = instr[1].header['DATE-END']

        ## This datetime stuff is not nessessary!!!!

        #dts = datetime.datetime.strptime(dstart, "%Y-%m-%dT%H:%M:%S.%fZ")
        #dte = datetime.datetime.strptime(dend, "%Y-%m-%dT%H:%M:%S.%fZ")

        #offset_s1 = datetime.timedelta(seconds=66.184)
        #offset_s2 = datetime.timedelta(seconds=67.184)

        #if quarter <14:
        #    date_obs_new = dts + offset_s1
        #    date_end_new = dte + offset_s1
        #if quarter > 14:
        #    date_obs_new = dts + offset_s2
        #    date_end_new = dte + offset_s2
        #if quarter == 14:
        #    date_obs_new = dts + offset_s1
        #    date_end_new = dte + offset_s2

        #instr[1].header['DATE-OBS'] = str(date_obs_new)[:-3] + 'Z'
        #instr[1].header['DATE-END'] = str(date_end_new)[:-3] + 'Z'

        instr.writeto(outfile)


        #if quarter == 14:
        #    livtime = instr[1].header['LIVTIME']
        #    livtime += 1.
        #    exposure += 1.

# end time

    if (status == 0):
        message = 'KEPTIMEFIX completed at'
    else:
        message = '\nKEPTIMEFIX aborted at'
    kepmsg.clock(message,logfile,verbose)

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Fix the time error in the target pixel files')

    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of FITS input target pixel file', type=str)
    parser.add_argument('outfile', help='Name of FITS target pixel file to output', type=str)

    parser.add_argument('--clobber', action='store_true', help='overwrite a file with the same name as outfile?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='keptimefix.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    cmdLine=True

    args = parser.parse_args()
    cmdLine=True

    keptimefix(args.infile,args.outfile,args.clobber,args.verbose,
        args.logfile,args.status,cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$keptimefix.par")
    t = iraf.IrafTaskFactory(taskname="keptimefix", value=parfile, function=keptimefix)
