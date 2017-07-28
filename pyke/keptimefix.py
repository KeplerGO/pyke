from .utils import PyKEArgumentHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from . import kepio, kepmsg


__all__ = ['keptimefix']


def keptimefix(infile, outfile=None, overwrite=False, verbose=False,
               logfile='keptimefix.log'):
    """
    keptimefix -- Correct a time stamp error in the target pixel files

    All Kepler light curve and target pixel files with version numbers 5.0
    contain an error in the time stamps. This was fixed in the light curve with
    version 5.0 (at MAST after May 2013). The timescale for fixing the target
    pixel files is unclear but in the mean time this script will fix the target
    pixel file time stamps and make the times consistent with the light curve
    files. The error in Q0-13 can be corrected by adding 66.184s. During Q14
    there was a leap second added Q15+ can be corrected by adding 67.184s. This
    tool fixes the time stamp accordingly.

    Parameters
    ----------
    infile : str
        The name of a Kepler target pixel file obtained from the MAST.
    outfile : str
        The name of the output FITS target pixel filefile. outfile will be a
        direct copy of infile but with the TIME column updates to be correct
        and consistent with the Kepler light curve files
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.
    """
    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPTIMEFIX -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPTIMEFIX started at', logfile, verbose)

    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPTIMEFIX: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    instr = pyfits.open(infile, 'readonly')

    creator =  instr[0].header['CREATOR']
    if creator.find('TargetPixelExporterPipelineModule') < 0:
        errmsg = 'ERROR -- KEPTIMEFIX: this file is not a target pixel file'
        kepmsg.err(logfile, errmsg, verbose)

    header_ext1 = instr[1].header.cards
    data_ext1 = instr[1].data

    fileversion = instr[0].header['FILEVER']
    if float(fileversion) > 4.0:
        errmsg = ('ERROR -- KEPTIMEFIX: no time fix needed for this file.'
                  ' FILEVER > 4.0')
        kepmsg.err(logfile, errmsg, verbose)

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
        offset = np.where(CADNUM <= 57139, 66.184, 67.184) / 86400.
    elif cadencetype == 'S':
        offset = np.where(CADNUM <= 1702663, 66.184, 67.184) / 86400.

    TIME_right = TIME_wrong + offset
    TIMECORR_new = TIMECORR_old + offset
    instr[1].data['TIME'][:] = TIME_right
    #we decided not to use the updated timecorr because
    #it is different from the LC FITS files by ~1 ms.
    instr[1].data['TIMECORR'][:] = np.nan * np.empty(len(TIMECORR_old))
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
    print("Writing output file {}...".format(outfile))
    instr.writeto(outfile)
    # end time
    kepmsg.clock('KEPTIMEFIX completed at', logfile, verbose)


def keptimefix_main():
    import argparse
    parser = argparse.ArgumentParser(
             description='Fix the time error in the target pixel files',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of FITS input target pixel file',
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS target pixel file to output.'
                              ' If None, outfile is infile-keptimefix.'),
                        default=None)
    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite a file with the same name as outfile?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='keptimefix.log', type=str)
    args = parser.parse_args()
    keptimefix(args.infile, args.outfile, args.overwrite, args.verbose,
               args.logfile)
