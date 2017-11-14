from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey
import numpy as np
from astropy.io import fits as pyfits

__all__ = ['kepclean']

def kepclean(infile, outfile=None, zero=True, overwrite=False, verbose=False,
             logfile='kepclean.log'):
    """
    Remove NaN values from a kepler light curve or TPF fits file. If passed a TPF
    only cadences where the postage stamp is ALL NaN values will be removed.

    Parameters
    ----------
    infile: str
         The name of a MAST standard format FITS file containing a Kepler light
         curve within the first data extension.
    outfile: str
         The name of the output FITS file with unwanted data removed.
    zero: bool
        Zero any remaining NaNs? (When passing TPFs, lone NaNs can exist in the
        corners of files.)
    overwrite: bool
        Overwrite the output file?
    verbose: bool
        Print informative messages and warnings to the shell and logfile?
    logfile: str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepclean kplr002436324-2009259160929_llc.fits
          --verbose --overwrite
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPCLEAN -- '
            + 'infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' zero={}'.format(zero)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPCLEAN started at',logfile,verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)

    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPCLIP: ' + outfile + ' exists. Use --overwrite'
        kepmsg.err(logfile, errmsg, verbose)

    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    #loop over all extentions in the file
    for ext in range(len(instr)):
        try:
            kepmsg.log(logfile, 'KEPCLEAN - Cleaning extention {} ({})'.format(ext,instr[ext].header['EXTNAME']), verbose)
        except:
            kepmsg.log(logfile, 'KEPCLEAN - Cleaning extention {}'.format(ext), verbose)
        #Dig out the header names that are data and find those that are 'FLUX'
        hnames = []
        flux_keys = []
        for k in (instr[ext].header):
            if k.startswith('TTYPE'):
                hname = instr[ext].header[k]
                try:
                    dat = instr[ext].data[instr[ext].header[k]]
                    hnames.append(hname)
                    if 'FLUX' in hname:
                        flux_keys.append(True)
                    else:
                        flux_keys.append(False)
                except:
                    continue

        flux_keys = np.asarray(flux_keys, dtype=bool)
        if np.any([(len(flux_keys) == 0), (len(np.where(flux_keys == True)[0]) == 0)]) :
            continue

        hnames = np.asarray(hnames)
        #if there are no FLUX keys to correct by, just move on
        kepmsg.log(logfile, 'KEPCLEAN - Found flux keys {}'.format(hnames[flux_keys]), verbose)
        #Loop through and remove indicies where any 'FLUX' keys have NaNs
        length = len(instr[ext].data[hnames[flux_keys][0]])
        fin = set(list(np.arange(length)))

        for y in hnames[flux_keys]:
            ydat = instr[ext].data[y]
            if len(np.shape(ydat)) > 1:
                for i in range(len(np.shape(ydat))-1):
                    bad = np.all(~np.isfinite(ydat),axis=-1)
                    ydat = np.nansum(ydat, axis=-1)
                    ydat[bad] = np.nan
            bad = set(list(np.where(~np.isfinite(ydat))[0]))
            if len(np.asarray(list(fin - bad))) == 0:
                kepmsg.log(logfile, 'WARNING - All {} values are nan. Cannot clean {}.'.format(y, y),
                    verbose)
                continue
            fin -= bad
        fin = np.asarray(list(fin))

        if len(fin) == 0:
            kepmsg.log(logfile, 'WARNING - All values are nan. Cannot clean.', verbose)
            continue

        kepmsg.log(logfile, 'KEPCLEAN - Removed {} nan entries out of {}'.format(length - len(fin), length), verbose)

        #Replace the extention
        if zero:
            for d in instr[ext].data.names:
                y = instr[ext].data[d]
                y[np.isfinite(y) == False] = 0
                instr[ext].data[d] = y

        instr[ext].data = instr[ext].data[fin]
        comment = 'NaN cadences removed from data'
        #Add a cleaned header keyword
        kepkey.new('NANCLEAN', True, comment, instr[ext], outfile, logfile, verbose)
    # write output file
    kepmsg.log(logfile, "KEPCLEAN - Writing output file {}...".format(outfile), verbose)
    # history keyword in output file
    kepmsg.log(logfile, "Writing output file {}...".format(outfile), verbose)
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    instr.writeto(outfile)
    instr.close()
    kepmsg.clock('KEPCLEAN completed at', logfile, verbose)

def kepclean_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Remove unwanted time'
                          ' ranges from Kepler time series data'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepclean.'),
                        default=None)
    parser.add_argument('--zero', action='store_true',
                        help='Zero any remaining NaN values?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepclean.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepclean(args.infile, args.outfile, args.zero, args.overwrite, args.verbose,
    args.logfile)
