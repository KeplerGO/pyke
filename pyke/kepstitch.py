from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from . import kepio, kepmsg, kepkey, kepstat


__all__ = ['kepstitch']


def kepstitch(infiles, outfile=None, overwrite=False, verbose=False,
              logfile='kepstich.log'):
    """
    kepstitch -- Append short cadence months and/or long cadence quarters

    ``kepstitch`` takes a list of FITS light curves obtained from the MAST
    archive and concatenates them into one long table. Be aware that this task
    currently does not sort the table by time. If users must ensure that the
    input file list is ordered by increasing time. Also no attempt is made to
    correct for aperture and PSF differences across quarterly rolls. Light
    curves spanning multiple quarters are likely to be discontinuous at the
    quarterly boundaries. Both issues will be addressed in the future.

    Parameters
    ----------
    infiles : list of str
        A list of individual FITS files. Most commonly these will contain data
        from the same target obtained over different months or quarters. Short
        cadence and long cadence data can be mixed within this list.
        Alternatively an ASCII file containing file names can parsed to this
        argument. The ASCII file must be formatted to have one FITS filename
        per line.
    outfile : str
        The name of the output FITS file with concatenated time series data.
    overwrite : bool
        Overwrite the output file? if ``overwrite = False`` and an existing
        file has the same name as outfile then the task will stop with an
        error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepstitch kplr012557548-2011271113734_llc.fits kplr012557548-2012004120508_llc.fits

    .. image:: ../_static/images/api/kepstitch/kplr012557548-2011271113734_llc.png
        :align: center

    .. image:: ../_static/images/api/kepstitch/kplr012557548-2012004120508_llc.png
        :align: center

    .. image:: ../_static/images/api/kepstitch/kplr012557548-stitch.png
        :align: center

    """
    if outfile is None:
        outfile = infiles[0].split('.')[0] + "-{}.fits".format(__all__[0])

    # startup parameters
    lct, bjd = [], []

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPSTITCH -- '
            + ' infiles={}'.format(infiles)
            + ' outfile={}'.format(outfile)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPSTITCH started at', logfile, verbose)
    # parse input file list
    try:
        infiles = kepio.parselist(infiles, logfile, verbose)
    except AttributeError:
        pass

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPSTITCH: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    # open output file
    outstr = pyfits.open(infiles[0], 'readonly')
    nrows1 = outstr[1].data.shape[0]
    head0 = outstr[0].header
    head1 = outstr[1].header

    # open input files
    nfiles = 0
    for infile in infiles:
        instr = pyfits.open(infile, 'readonly')
        # append table data
        if nfiles > 0:
            nrows2 = instr[1].data.shape[0]
            nrows = nrows1 + nrows2
            outtab = pyfits.BinTableHDU.from_columns(outstr[1].columns, nrows=nrows)
            for name in outstr[1].columns.names:
                try:
                    outtab.data.field(name)[nrows1:]=instr[1].data.field(name)
                except:
                    warnmsg = ('ERROR -- KEPSTITCH: column {} missing from'
                               ' some files.'.format(name))
                    kepmsg.warn(logfile, warnmsg, verbose)
                    pass
            outstr[1] = outtab
            outstr[0].header = head0
            outstr[1].header = head1
            nrows1 = nrows

        # start and stop times of data
        fitsvers = 1.0
        lc_start = kepkey.get(infile, instr[1], 'LC_START', logfile, verbose)
        lc_end = kepkey.get(infile, instr[1], 'LC_END', logfile, verbose)
        try:
            startbjd = instr[1].header['STARTBJD']
        except:
            startbjd = kepkey.get(infile, instr[1], 'TSTART', logfile, verbose)
            fitsvers = 2.0
        try:
            endbjd = instr[1].header['ENDBJD']
        except:
            endbjd = kepkey.get(infile, instr[1], 'TSTOP', logfile, verbose)
            fitsvers = 2.0
        lct.append(lc_start)
        lct.append(lc_end)
        bjd.append(startbjd)
        bjd.append(endbjd)

        # close input files
        instr.close()
        nfiles += 1

    # maxmimum and minimum times in file sample
    lc_start = np.min(lct)
    lc_end = np.max(lct)
    startbjd = np.min(bjd)
    endbjd = np.max(bjd)
    kepkey.change('LC_START', lc_start, outstr[1], outfile, logfile, verbose)
    kepkey.change('LC_END', lc_end, outstr[1], outfile, logfile, verbose)

    if fitsvers == 1.0:
        kepkey.change('STARTBJD', startbjd, outstr[1], outfile, logfile,
                      verbose)
        kepkey.change('ENDBJD', endbjd, outstr[1], outfile, logfile, verbose)
    else:
        kepkey.change('TSTART', startbjd, outstr[1], outfile, logfile, verbose)
        kepkey.change('TSTOP', endbjd, outstr[1], outfile, logfile, verbose)

    # comment keyword in output file
    kepkey.comment(call, outstr[0], outfile, logfile, verbose)
    # close output file
    print("Writing output file {}...".format(outfile))
    outstr.writeto(outfile)
    outstr.close()
    ## end time
    kepmsg.clock('KEPSTITCH completed at', logfile, verbose)

def kepstitch_main():
    import argparse

    parser = argparse.ArgumentParser(
             description=('Append multiple month short cadence and/or'
                          ' multiple quarter long cadence data'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infiles', help='List of input files', nargs='+',
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepstitch.'),
                        default=None)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepstitch.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepstitch(args.infiles, args.outfile, args.overwrite, args.verbose,
              args.logfile)
