from .utils import PyKEArgumentHelpFormatter
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from . import kepio, kepmsg, kepkey


__all__ = ['kephead']


def kephead(infile, keyname, outfile=None, overwrite=False, verbose=False,
            logfile='kephead.log'):
    """
    kephead -- search for and list FITS keywords in Kepler data files

    The Kepler archive contains science data in Flexible Image Transport System
    (FITS) format. Images are stored in binary arrays, time series data in
    binary tables and scalar data as binary keywords. Tools such as ``kepmask``
    and ``kepconvert`` read Kepler image and table data. ``kephead`` displays
    scalar data or converts it from binary to ASCII format. ``kephead`` can
    display all keywords in a file, display individual keywords, or search for
    character combinations in keywords names.

    Parameters
    ----------
    infile : str
        The name of a standard format FITS file.
    keyname : str
        The name of a keyword, value and description to print or store within
        the output file. Partial completions for the keyword are allowed. For
        example ``keyname='mag'`` will return all instances of magnitude
        keywords (BMAG, VMAG, KEPMAG etc) and their values within the light
        curves stored at MAST. ``keyname='all'`` will return every keyword and
        their values stored in all extensions within the input FITS file.
    outfile : str
        The name of the output ASCII file for storing search results.
    overwrite : bool
        Overwrite the output file? if overwrite is False and an existing file
        has the same name as outfile then the task will stop with an error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kephead kplr002436324-2009259160929_llc.fits 'mag' --verbose --overwrite
        --------------------------------------------------------------
        KEPHEAD --  infile=kplr002436324-2009259160929_llc.fits outfile=.txt keyname=mag
        overwrite=True verbose=True logfile=kephead.log

        KEPHEAD started at: : Tue Jun 13 15:26:58 2017

        ---------------------------------------
        kplr002436324-2009259160929_llc.fits[0]
        ---------------------------------------

        GMAG    =               16.334 / [mag] SDSS g band magnitude
        RMAG    =               16.789 / [mag] SDSS r band magnitude
        IMAG    =               17.112 / [mag] SDSS i band magnitude
        ZMAG    =                      / [mag] SDSS z band magnitude
        D51MAG  =                      / [mag] D51 magnitude,
        JMAG    =                8.176 / [mag] J band magnitude from 2MASS
        HMAG    =                7.262 / [mag] H band magnitude from 2MASS
        KMAG    =                6.874 / [mag] K band magnitude from 2MASS
        KEPMAG  =               16.879 / [mag] Kepler magnitude (Kp)


        KEPHEAD ended at: : Tue Jun 13 15:26:58 2017
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.txt".format(__all__[0])
    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPHEAD -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' keyname={}'.format(keyname)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPHEAD started at: ',logfile,verbose)
    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPHEAD: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    # open input FITS file
    instr = pyfits.open(infile)
    # number of HDU in input file
    nhdu = kepkey.HDUnum(instr)
    # loop through each HDU in infile
    kepmsg.log(outfile, '', True)
    print("Writing output file {}...".format(outfile))
    for hdu in tqdm(range(nhdu)):
        # how many keywords in the HDU?
        keylist = instr[hdu].header.cards
        nkeys = len(keylist)
        # print header number
        prhead = False
        for i in range(nkeys):
            if (keyname.upper() == 'ALL'
                or keyname.upper() in list(instr[hdu].header.keys())[i]):
                prhead = True
        if prhead:
            dashes = ''
            title = infile + '[' + str(hdu) + ']'
            for j in range(len(title)):
                dashes += '-'
            kepmsg.log(outfile, dashes, True)
            kepmsg.log(outfile, title, True)
            kepmsg.log(outfile, dashes, True)
            kepmsg.log(outfile, '', True)
        # print keywords
        for i in range(nkeys):
            if ((keyname.upper() == 'ALL'
                  or keyname.upper() in list(instr[hdu].header.keys())[i])
                and 'COMMENT' not in list(instr[hdu].header.keys())[i]):
                kepmsg.log(outfile, str(keylist[i]), True)
        kepmsg.log(outfile, '', True)
    # stop time
    kepmsg.clock('KEPHEAD ended at: ', logfile, verbose)

def kephead_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Search for and list FITS keywords in Kepler data'
                         ' files'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('keyname', help='Snippet of keyword name', type=str)
    parser.add_argument('--outfile',
                        help=('The name of the output ASCII file for storing search results.'
                              ' If None, outfile is infile-kephead.'),
                        default=None)
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='kephead.log', type=str)
    args = parser.parse_args()
    kephead(args.infile, args.keyname, args.outfile, args.overwrite,
            args.verbose, args.logfile)
