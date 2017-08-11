from .utils import PyKEArgumentHelpFormatter
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm
from . import kepio, kepmsg, kepkey, kepfunc


__all__ = ['kepsmooth']


def kepsmooth(infile, outfile=None, datacol='SAP_FLUX', function='flat',
              fscale=1.0, plot=False, overwrite=False, verbose=False,
              logfile='kepsmooth.log'):
    """
    kepsmooth -- Smooth Kepler light curve data by convolution with a choice
    of analytical functions.

    The smoothed data is copied to a new FITS file with the same structure as
    the input file.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler light
        curve data within the first data extension.
    outfile : str
        The name of the output FITS file. The output file is identical in
        format to the input file. The data to be smoothed will be overwritten
        in the output file by its convolved version.
    datacol : str
        The name of the data column in the input FITS file to be smoothed, e.g.
        ``SAP_FLUX``, ``PDCSAP_FLUX``, ``MOM_CENTR1`` etc. A full list of
        archived data columns is provided in the Kepler Archive Manual.
    function : str
        The form of the smoothing function. A flat function is a moving
        average. The options are:

            * flat
            * hanning
            * hamming
            * bartlett
            * blackman

    fscale : float
        The width of the convolution function in units of days.
    plot : bool
        Plot the original light curve and the result of the smoothing?
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
         Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepsmooth kplr005110407-2009259160929_llc.fits --plot

    .. image:: ../_static/images/api/kepsmooth.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = ('KEPSMOOTH -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' function={}'.format(function)
            + ' fscale={}'.format(fscale)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    ## start time
    kepmsg.clock('KEPSMOOTH started at', logfile, verbose)

    ## overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPSMOOTH: {} exists. Use overwrite=True'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    ## open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    if cadence == 0.0:
        tstart, tstop, ncad, cadence = kepio.cadence(instr, infile, logfile,
                                                     verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence

    ## fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    ## read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # read time and flux columns
    barytime = kepio.readtimecol(infile, table, logfile, verbose)
    flux = kepio.readfitscol(infile, instr[1].data, datacol, logfile, verbose)
    # filter input data table
    try:
        nanclean = instr[1].header['NANCLEAN']
    except:
        naxis2 = 0
        for i in tqdm(range(len(table.field(0)))):
            if (np.isfinite(barytime[i]) and np.isfinite(flux[i])
                and flux[i] != 0.0):
                table[naxis2] = table[i]
                naxis2 += 1
        instr[1].data = table[:naxis2]
        comment = 'NaN cadences removed from data'
        kepkey.new('NANCLEAN', True, comment, instr[1], outfile, logfile,
                   verbose)

    ## read table columns
    try:
        intime = instr[1].data.field('barytime')
    except:
        intime = kepio.readfitscol(infile, instr[1].data, 'time', logfile,
                                   verbose)
    indata = kepio.readfitscol(infile, instr[1].data, datacol, logfile,
                               verbose)
    intime = intime + bjdref
    indata = indata / cadenom

    ## smooth data
    outdata = kepfunc.smooth(indata, fscale/(cadence/86400), function)
    ## comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    ## clean up x-axis unit
    intime0 = float(int(tstart / 100) * 100.0)
    if intime0 < 2.4e6:
        intime0 += 2.4e6
    ptime = intime - intime0
    xlab = 'BJD $-$ {0}'.format(intime0)
    ## clean up y-axis units
    pout = indata * 1.0
    pout2 = outdata * 1.0
    nrm = len(str(int(np.nanmax(pout)))) - 1
    pout = pout / 10 ** nrm
    pout2 = pout2 / 10 ** nrm
    ylab = '10$^{0}$ {1}'.format(nrm, 'e$^-$ s$^{-1}$')
    ## data limits
    xmin = np.nanmin(ptime)
    xmax = np.nanmax(ptime)
    ymin = np.min(pout)
    ymax = np.nanmax(pout)
    xr = xmax - xmin
    yr = ymax - ymin
    ptime = np.insert(ptime, [0], [ptime[0]])
    ptime = np.append(ptime, [ptime[-1]])
    pout =  np.insert(pout, [0], [0.0])
    pout =  np.append(pout, 0.0)
    pout2 = np.insert(pout2, [0], [0.0])
    pout2 = np.append(pout2, 0.0)
    ## plot light curve
    if plot:
        plt.figure()

        # delete any fossil plots in the matplotlib window
        plt.clf()

        # position axes inside the plotting window
        ax = plt.subplot(111)
        plt.subplots_adjust(0.06, 0.1, 0.93, 0.88)

        # force tick labels to be absolute rather than relative
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # rotate y labels by 90 deg
        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90)
        plt.plot(ptime[1:-1], pout[1:-1], color='#ff9900', linestyle='-',
                 linewidth=1.0)
        plt.fill(ptime, pout, color='#ffff00', linewidth=0.0, alpha=0.2)
        plt.plot(ptime, pout2, color='#0000ff', linestyle='-', linewidth=1.0*4.0)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
        if ymin >= 0.0:
            plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)
        else:
            plt.ylim(1.0e-10, ymax+yr*0.01)
        plt.grid()

        # render plot
        plt.show()

    ## write output file
    print("Writing output file {}...".format(outfile))
    for i in tqdm(range(len(outdata))):
        instr[1].data.field(datacol)[i] = outdata[i]
    instr.writeto(outfile)
    ## close input file
    instr.close()
    ## end time
    kepmsg.clock('KEPSMOOTH completed at', logfile, verbose)

def kepsmooth_main():
    import argparse
    parser = argparse.ArgumentParser(
             description='Smooth Kepler light curve data by convolution',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepsmooth.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column to plot', type=str)
    parser.add_argument('--function', default='hanning',
                        help='Type of convolution filter', type=str,
                        choices=['flat', 'hanning', 'hamming', 'bartlett',
                                 'blackman'])
    parser.add_argument('--fscale', default=1.0,
                        help=('Characteristic width of smoothing function'
                              ' [days]'), type=float)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='kepsmooth.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepsmooth(args.infile, args.outfile, args.datacol, args.function,
              args.fscale, args.plot, args.overwrite, args.verbose,
              args.logfile)
