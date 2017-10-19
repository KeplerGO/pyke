from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey, kepfunc, kepstat
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm


__all__ = ['kepfilter']


def kepfilter(infile, passband, outfile=None, datacol='SAP_FLUX', function='boxcar',
              cutoff=1.0, plot=False, overwrite=False, verbose=False,
              logfile='kepfilter.log'):
    """
    kepfilter -- bandpass filtering of Kepler light curve data

    ``kepfilter`` applies a bandpass filter to Kepler light curve data. In the
    low bandpass option, the data is convolved with a function of
    user-specified width. Choices of convolution function are **boxcar**,
    **Gaussian** or **sinc**. In the high bandpass option the convolution minus
    the median of the convolution is subtracted from the original data. The
    filtered data is copied to a new FITS file with the same structure as the
    input file.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing Kepler light
        curve data within the first data extension.
    passband : str
        The type of filter to be applied. A low bandpass filter will suppress
        high-frequency signal shorter than the cutoff. A high bandpass filter
        will suppress low-frequency signal longer than the cutoff.
        The options are:

        * low
        * high
    outfile : str
        The name of the output FITS file. The output file is identical in
        format to the input file. The data to be filtered will be overwritten
        in the output file by its filtered version.
    datacol : str
        The name of the data column in the input FITS file to be filtered, e.g.
        SAP_FLUX, PDCSAP_FLUX, MOM_CENTR1 etc. A full list of
        archived data columns is provided in the Kepler Archive Manual.
    function : string
        The functional form of the bandpass convolution function.
        The options are:

        * boxcar
        * gauss
        * sinc
    cutoff : float
        The frequency of the bandpass cutoff in units of days-1.
    plot : bool
        Plot the original light curve and the result of the filter?
    overwrite : bool
        Overwrite the output file? if overwrite is **False** and an existing
        file has the same name as outfile then the task will stop with an
        error.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------

    .. code-block :: bash

        $ kepfilter kplr002436324-2009259160929_llc.fits --datacol 'SAP_FLUX' --function 'boxcar'
        --plot --verbose --overwrite

    .. image :: ../_static/images/api/kepfilter.png
        :align: center
    """
    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.fits".format(__all__[0])
    ## log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPFILTER -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' function={}'.format(function)
            + ' cutoff={}'.format(cutoff)
            + ' passband={}'.format(passband)
            + ' plot={}'.format(plot)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    ## start time
    kepmsg.clock('KEPFILTER started at',logfile,verbose)
    ## overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPFILTER: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, message, verbose)

    ## open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)
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
    flux= kepio.readsapcol(infile, table, logfile, verbose)
    # filter input data table
    try:
        nanclean = instr[1].header['NANCLEAN']
    except:
        naxis2 = 0
        for i in range(len(table.field(0))):
            if (np.isfinite(barytime[i]) and np.isfinite(flux[i])
                and flux[i] != 0.0):
                table[naxis2] = table[i]
                naxis2 += 1
        instr[1].data = table[:naxis2]
        kepkey.new('NANCLEAN', True, 'NaN cadences removed from data',
                   instr[1], outfile, logfile, verbose)

    ## read table columns
    intime = (kepio.readtimecol(infile, instr[1].data, logfile, verbose)
              + bjdref)
    indata = kepio.readfitscol(infile, instr[1].data, datacol, logfile,
                               verbose) / cadenom
    ## define data sampling
    tr = 1.0 / (cadence / 86400)
    timescale = 1.0 / (cutoff / tr)
    ## define convolution function
    if function == 'boxcar':
        filtfunc = np.ones(int(np.ceil(timescale)))
    elif function == 'gauss':
        timescale /= 2
        dx = np.ceil(timescale * 10 + 1)
        filtfunc = kepfunc.gauss([1.0, dx / 2 - 1.0, timescale],
                                 np.linspace(0, dx - 1, dx))
    elif function == 'sinc':
        dx = np.ceil(timescale * 12 + 1)
        fx = (np.linspace(0, dx - 1, dx) - dx / 2 + 0.5) / timescale
        filtfunc = np.sinc(fx)

    filtfunc /= np.sum(filtfunc)
    ## pad time series at both ends with noise model
    ave, sigma = (np.mean(indata[:len(filtfunc)]),
                  np.std(indata[:len(filtfunc)]))
    padded = np.append(kepstat.randarray(np.ones(len(filtfunc)) * ave,
                       np.ones(len(filtfunc)) * sigma), indata)
    ave, sigma = (np.mean(indata[-len(filtfunc):]),
                  np.std(indata[-len(filtfunc):]))
    padded = np.append(padded, kepstat.randarray(np.ones(len(filtfunc)) * ave,
                       np.ones(len(filtfunc)) * sigma))
    ## convolve data
    convolved = np.convolve(padded,filtfunc,'same')
    ## remove padding from the output array
    if function == 'boxcar':
        outdata = convolved[len(filtfunc):-len(filtfunc)]
    else:
        outdata = convolved[len(filtfunc):-len(filtfunc)]
    ## subtract low frequencies
    if passband == 'high':
        outmedian = np.median(outdata)
        outdata = indata - outdata + outmedian
    ## comment keyword in output file
    kepkey.history(call, instr[0], outfile, logfile, verbose)
    ## clean up x-axis unit
    intime0 = float(int(tstart / 100) * 100.0)
    if intime0 < 2.4e6: intime0 += 2.4e6
    ptime = intime - intime0
    xlab = 'BJD $-$ {}'.format(intime0)
    ## clean up y-axis units
    pout = indata * 1.0
    pout2 = outdata * 1.0
    nrm = len(str(int(np.nanmax(pout)))) - 1
    pout = pout / 10 ** nrm
    pout2 = pout2 / 10 ** nrm
    ylab = '10$^{}$ {}'.format(nrm, 'e$^-$ s$^{-1}$')
    ## data limits
    xmin = ptime.min()
    xmax = ptime.max()
    ymin = np.nanmin(pout)
    ymax = np.nanmax(pout)
    xr = xmax - xmin
    yr = ymax - ymin
    ptime = np.insert(ptime, [0], [ptime[0]])
    ptime = np.append(ptime, [ptime[-1]])
    pout = np.insert(pout, [0], [0.0])
    pout = np.append(pout, 0.0)
    pout2 = np.insert(pout2, [0], [0.0])
    pout2 = np.append(pout2, 0.0)
    ## plot light curve
    if plot:
        plt.figure()
        plt.clf()

        ## plot filtered data
        ax = plt.axes([0.06, 0.1, 0.93, 0.87])
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.plot(ptime, pout, color='#ff9900', linestyle='-', linewidth=1.0)
        plt.fill(ptime, pout, color='#ffff00', linewidth=0.0, alpha=0.2)
        if passband == 'low':
            plt.plot(ptime[1:-1], pout2[1:-1], color='#0000ff', linestyle='-',
                     linewidth=1.0)
        else:
            plt.plot(ptime, pout2, color='#0000ff', linestyle='-',
                     linewidth=1.0)
            plt.fill(ptime, pout2, color='#0000ff', linewidth=0.0, alpha=0.2)
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(ylab, {'color' : 'k'})
        plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0:
            plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            plt.ylim(1.0e-10,ymax+yr*0.01)
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
    kepmsg.clock('KEPFILTER completed at', logfile, verbose)

def kepfilter_main():
    import argparse

    parser = argparse.ArgumentParser(
             description='Low bandpass or high bandpass signal filtering',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--passband', help='low- or high-bandpass filter',
                        type=str, choices=['low','high'])
    parser.add_argument('--outfile',
                        help=('Name of FITS file to output.'
                              ' If None, outfile is infile-kepfilter.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of data column', type=str)
    parser.add_argument('--function', default='boxcar',
                        help='The bandpass convolution function', type=str,
                        choices=['boxcar','gauss','sinc'])
    parser.add_argument('--cutoff', default=1.0,
                        help='Characteristic frequency cutoff of filter [1/days]',
                        type=float)
    parser.add_argument('--plot', action='store_true',
                        help='Plot result?')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file',
                        default='kepfilter.log', type=str)
    args = parser.parse_args()
    kepfilter(args.infile, args.passband, args.outfile, args.datacol,
              args.function, args.cutoff, args.plot, args.overwrite,
              args.verbose, args.logfile)
