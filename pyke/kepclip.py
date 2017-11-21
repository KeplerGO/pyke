from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from tqdm import tqdm


__all__ = ['kepclip']


def kepclip(infile, ranges, outfile=None, datacol='SAP_FLUX', plot=False,
            overwrite=False, verbose=False, logfile='kepclip.log'):
    """
    Remove unwanted time ranges from Kepler time series data.

    Parameters
    ----------
    infile: str
         The name of a MAST standard format FITS file containing a Kepler light
         curve within the first data extension.
    outfile: str
         The name of the output FITS file with unwanted data removed.
    ranges: str
       Time ranges are supplied as comma-separated pairs of Barycentric Julian
       Dates (BJDs). Multiple ranges are separated by a semi-colon. An example
       containing two time ranges is:
       ``'2455012.48517,2455014.50072;2455022.63487,2455025.08231'``
    plot: bool
        Plot the output data?
    datacol: str
        Types of photometry stored in the input file.
    overwrite: bool
        Overwrite the output file?
    verbose: bool
        Print informative messages and warnings to the shell and logfile?
    logfile: str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ kepclip kplr002436324-2009259160929_llc.fits
          '2455012.48517,2455018.50072;2455022.63487,2455060.08231'
          --verbose --plot --overwrite

    .. image:: ../_static/images/api/kepclip.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPCLIP -- '
            + 'infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' ranges={}'.format(ranges)
            + ' plot={}'.format(plot)
            + ' datacol={}'.format(datacol)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPCLIP started at', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPCLIP: ' + outfile + ' exists. Use --overwrite'
        kepmsg.err(logfile, errmsg, verbose)

    # time ranges for region
    t1 = []; t2 = []
    t1, t2 = kepio.timeranges(ranges, logfile, verbose)

    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile, logfile,
                                                    verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # input data
    table = instr[1].data
    # read time and flux columns
    barytime = kepio.readtimecol(infile, table, logfile, verbose)
    barytime = barytime + bjdref

    #Test file type is LC or TPF:
    if len(set(instr[1].data.columns.names) & set(['SAP_FLUX'])) == 1:
        filetype='LC'
    if len(set(instr[1].data.columns.names) & set(['SAP_FLUX'])) == 0:
        filetype='TPF'

    if filetype == 'LC':
        message = 'KEPCLIP clipping a Light Curve'
        kepmsg.clock(message, logfile, verbose)
        flux = kepio.readfitscol(infile, table, datacol, logfile, verbose)
        if 'flux' in datacol.lower():
            flux = flux / cadenom
        # filter input data table
        naxis2 = 0
        work1 = np.array([], 'float64')
        work2 = np.array([], 'float32')
        for i in tqdm(range(len(barytime))):
            if (np.isfinite(barytime[i]) and np.isfinite(flux[i])
                and flux[i] != 0.0):
                reject = True
                for j in range(len(t1)):
                    if (barytime[i] >= t1[j] and barytime[i] <= t2[j]):
                        reject = False
                if not reject:
                    table[naxis2] = table[i]
                    work1 = np.append(work1, barytime[i])
                    work2 = np.append(work2, flux[i])
                    naxis2 += 1

        # comment keyword in output file
        kepkey.history(call, instr[0], outfile, logfile, verbose)
        # comment keyword in output file
        kepmsg.log(logfile, "Writing output file {}...".format(outfile), verbose)
        # write output file
        instr[1].data = table[:naxis2]
        comment = 'NaN cadences removed from data'
        kepkey.new('NANCLEAN', True, comment, instr[1], outfile, logfile, verbose)
        instr.writeto(outfile)
        # clean up x-axis unit
        barytime0 = (tstart // 100) * 100.0
        barytime = work1 - barytime0
        xlab = 'BJD $-$ {}'.format(barytime0)
        # clean up y-axis units
        try:
            nrm = len(str(int(work2.max()))) - 1
        except:
            nrm = 0
        flux = work2 / 10 ** nrm
        ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm
        # data limits
        xmin = barytime.min()
        xmax = barytime.max()
        ymin = flux.min()
        ymax = flux.max()
        xr = xmax - xmin
        yr = ymax - ymin
        # clear window, plot box
        if plot:
            plt.figure()
            plt.clf()
            ax = plt.axes()
            # force tick labels to be absolute rather than relative
            plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
            # rotate y labels by 90 deg
            # plot line data
            ltime = [barytime[0]]; ldata = [flux[0]]
            for i in range(1, len(flux)):
                if barytime[i - 1] > barytime[i] - 0.025:
                    ltime.append(barytime[i])
                    ldata.append(flux[i])
                else:
                    ltime = np.array(ltime, dtype=np.float)
                    ldata = np.array(ldata, dtype=np.float)
                    plt.plot(ltime, ldata, color='#0000ff', linestyle='-',
                             linewidth=1.0)
                    ltime = []; ldata = []
            ltime = np.array(ltime, dtype=np.float)
            ldata = np.array(ldata, dtype=np.float)
            plt.plot(ltime, ldata, color='#0000ff', linestyle='-', linewidth=1.0)
            # plot fill data
            barytime = np.insert(barytime, [0], [barytime[0]])
            barytime = np.append(barytime, [barytime[-1]])
            flux = np.insert(flux, [0], [0.0])
            flux = np.append(flux, [0.0])
            plt.fill(barytime, flux, fc='#ffff00', linewidth=0.0, alpha=0.2)
            plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
            if ymin - yr * 0.01 <= 0.0:
                plt.ylim(1.0e-10, ymax + yr * 0.01)
    if filetype == 'TPF':
        message = 'KEPCLIP clipping a Target Pixel File'
        kepmsg.clock(message, logfile, verbose)
        # filter input data table
        naxis2 = 0

        for i in tqdm(range(len(barytime))):
            if (np.isfinite(barytime[i])):
                reject = True
                for j in range(len(t1)):
                    if (barytime[i] >= t1[j]) and (barytime[i] < t2[j]):
                        reject = False
                if not reject:
                    table[naxis2] = table[i]
                    naxis2 += 1

        # comment keyword in output file
        kepkey.history(call, instr[0], outfile, logfile, verbose)

        # write output file
        instr[1].data = table[:naxis2]
        comment = 'Trimmed TPF'
        kepkey.new('CLIP_TPF', True, comment, instr[1], outfile, logfile, verbose)
        instr.writeto(outfile)
        if plot:
            message = 'KEPCLIP no plotting available for Target Pixel Files'
            kepmsg.clock(message, logfile, verbose)
    # close input file
    instr.close()

    # end time
    message = 'KEPCLIP completed at'
    kepmsg.clock(message, logfile, verbose)

def kepclip_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Remove unwanted time'
                          ' ranges from Kepler time series data'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('ranges',
                        help='List of time domain ranges to be excluded',
                        type=str)
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepclip.'),
                        default=None)
    parser.add_argument('--datacol', help='Data column to plot',
                        default='SAP_FLUX', type=str)
    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepclip.log', dest='logfile', type=str)
    args = parser.parse_args()
    kepclip(args.infile, args.ranges, args.outfile, args.datacol, args.plot,
            args.overwrite, args.verbose, args.logfile)
