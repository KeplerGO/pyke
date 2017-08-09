from .utils import PyKEArgumentHelpFormatter
import re
import numpy as np
import sys
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from . import kepio, kepmsg, kepkey


__all__ = ['kepdraw']


def kepdraw(infile, outfile=None, datacol=None, ploterr=False,
            errcol='SAP_FLUX_ERR', quality=False, lcolor='#363636', lwidth=1.0,
            fcolor='#a8a7a7', falpha=0.2, labelsize=20, ticksize=14, xsize=18.,
            ysize=6., fullrange=False, chooserange=False, y1=0, y2=1e4,
            plotgrid=False, ylabel='e$^-$ s$^{-1}$', plottype='fast',
            noninteractive=False, verbose=False, logfile='kepdraw.log'):
    """
    kepdraw -- Interactive plotting of time series data

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : string
        The name of the output PNG plot file.
    datacol : string
        The FITS column within extension 1 of infile to be plotted against
        time.
    ploterr : bool
        Plot error bars with the data?
    errcol : str
        The uncertainty data coupled to datacol. If ploterr is ``True`` then
        these data will be plotted as error bars
    quality : bool
        Ignore bad data and cadence gaps in datacol?
    lcolor : str
        The color of the line plot. These can be chosen in html notation, e.g.
        black is #000000 and white is #ffffff.
    lwidth : float
        The width of the plot line in arbitrary units. The plotting default
        width is 1.0.
    fcolor : str
        The color of the fill under the plot line. These can be chosen in html
        notation, e.g. black is #000000 and white is #ffffff. The html color
        scheme can be browsed in multiple online charts.
    falpha : float
        The transparency of the fill color. 0.0 is transparent, 1.0 is opaque.
    labelsize : float
        The font size of the plot labels.
    ticksize : float
        The axes tick sizes within the plot.
    xsize : float
        The length of the plotting window in arbitrary units.
        A typical value is 16.0
    ysize : float
        The height of the plotting window in arbitrary units.
        A typical value is 8.0
    fullrange : bool
        Plot the flux range between zero and the maximum and 101% of the flux
        maximum? Otherwise plot the range 99% of the flux minimum to 101% of
        the flux maximum.
    plotgrid : bool
        Plot a grid over the data on the major tick marks?
    plottype : str
        The look and feel of the plot. 'fast' provides a rapid quick look plot
        using discrete points to represent the data. 'pretty' provides a
        continuous curve between data points, except at natural data gaps or
        where quality issues occur. Pretty curves are significantly slower to
        build. Options:

            * fast
            * pretty
    non-interactive : bool
        If True, prevents the matplotlib window to pop up.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------

    .. code-block:: bash

        $ kepdraw kplr010666592-2009131110544_slc.fits

    .. image:: ../_static/images/api/kepdraw.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.png".format(__all__[0])

    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPDRAW -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' datacol={}'.format(datacol)
            + ' ploterr={}'.format(ploterr)
            + ' errcol={}'.format(errcol)
            + ' quality={}'.format(quality)
            + ' lcolor={}'.format(lcolor)
            + ' lwidth={}'.format(lwidth)
            + ' fcolor={}'.format(fcolor)
            + ' falpha={}'.format(falpha)
            + ' labelsize={}'.format(labelsize)
            + ' ticksize={}'.format(ticksize)
            + ' xsize={}'.format(xsize)
            + ' ysize={}'.format(ysize)
            + ' fullrange={}'.format(fullrange)
            + ' chooserange={}'.format(chooserange)
            + ' ymin={}'.format(y1)
            + ' ymax={}'.format(y2)
            + ' plotgrid={}'.format(plotgrid)
            + ' ylabel={}'.format(ylabel)
            + ' plottype={}'.format(plottype)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)

    # start time
    kepmsg.clock('KEPDRAW started at', logfile, verbose)

    # open input file
    struct = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(struct, infile,
                                                    logfile, verbose)

    # read table structure
    table = kepio.readfitstab(infile, struct[1], logfile, verbose)

    # read table columns
    intime = kepio.readtimecol(infile, table, logfile, verbose) + bjdref

    if datacol is None:
        default_priority_columns = ['DETSAP_FLUX', 'PDCSAP_FLUX', 'SAP_FLUX']
        flux_columns_found = [c for c in struct[1].data.columns.names
                              if 'FLUX' in c]

        if len(flux_columns_found) == 0:
            raise ValueError("No flux columns found.")
        print("Found the following flux columns: {}".format(", ".join(flux_columns_found)))

        for colname in default_priority_columns:
            if colname in flux_columns_found:
                datacol = colname
                print("Using data column {} on the plot...".format(datacol))
                break
        else:
            datacol = flux_columns_found[0]
            print("Using data column {} on the plot...".format(datacol))

    indata = kepio.readfitscol(infile, table, datacol, logfile, verbose)
    indataerr = kepio.readfitscol(infile, table, errcol, logfile, verbose)
    qualty = kepio.readfitscol(infile, table, 'SAP_QUALITY', logfile, verbose)

    # close infile
    struct.close()

    # remove infinities and bad data
    if np.isnan(np.nansum(indataerr)):
        indataerr[:] = 1.0e-5
    work1 = np.array([intime, indata, indataerr, qualty],dtype='float64')
    work1 = np.rot90(work1, 3)
    work1 = work1[~np.isnan(work1).any(1)]
    work1 = work1[~np.isinf(work1).any(1)]
    if quality:
        work1 = work1[work1[:, 0] == 0.0]
    barytime = np.array(work1[:, 3], dtype='float64')
    data = np.array(work1[:, 2], dtype='float32')
    dataerr = np.array(work1[:, 1], dtype='float32')
    if len(barytime) == 0:
        message = 'ERROR -- KEPDRAW: Plotting arrays are full of NaN'
        kepmsg.err(logfile, message, verbose)

    # clean up x-axis unit
    barytime0 = float(int(tstart / 100) * 100.0)
    barytime -= barytime0
    xlab = 'BJD $-$ %d' % barytime0

    # clean up y-axis units
    nrm = 0
    try:
        nrm = len(str(int(np.nanmax(data))))-1
    except:
        nrm = 0
    data = data / 10**nrm
    if 'e$^-$ s$^{-1}$' in ylabel or 'default' in ylabel:
        if nrm == 0:
            ylab1 = 'e$^-$ s$^{-1}$'
        else:
            ylab1 = '10$^{%d}$ e$^-$ s$^{-1}$' % nrm
    else:
        ylab1 = re.sub('_','-',ylabel)

    # data limits
    xmin = np.nanmin(barytime)
    xmax = np.nanmax(barytime)
    ymin = np.nanmin(data)
    ymax = np.nanmax(data)
    xr = xmax - xmin
    yr = ymax - ymin
    barytime = np.insert(barytime,[0],[barytime[0]])
    barytime = np.append(barytime,[barytime[-1]])
    data = np.insert(data,[0],[-10000.0])
    data = np.append(data,-10000.0)
    # define size of plot on monitor screen
    plt.figure(figsize=[xsize,ysize])
    # delete any fossil plots in the matplotlib window
    plt.clf()
    # position axes inside the plotting window
    ax = plt.subplot(111)
    plt.subplots_adjust(0.06,0.15,0.92,0.83)
    # force tick labels to be absolute rather than relative
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    # rotate y labels by 90 deg
    labels = ax.get_yticklabels()
    # ifplot type is 'fast' plot data time series as points
    if plottype == 'fast':
        plt.plot(barytime, data, 'o', color=lcolor, markersize=lwidth)
    # ifplot type is 'pretty' plot data time series as an unbroken line,
    # retaining data gaps
    else:
        ltime = np.array([], dtype='float64')
        ldata = np.array([], dtype='float32')
        dt = 0
        work1 = 2.0 * cadence / 86400
        for i in range(1, len(data) - 1):
            dt = barytime[i] - barytime[i - 1]
            if dt < work1:
                ltime = np.append(ltime, barytime[i])
                ldata = np.append(ldata, data[i])
            else:
                plt.plot(ltime, ldata, color=lcolor, linestyle='-', linewidth=lwidth)
                ltime = np.array([], dtype='float64')
                ldata = np.array([], dtype='float32')
        plt.plot(ltime, ldata, color=lcolor, linestyle='-', linewidth=lwidth)

    # plot the fill color below data time series, with no data gaps
    plt.fill(barytime, data, fc=fcolor, linewidth=0.0, alpha=falpha)

    # define plot x and y limits
    plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
    if ymin-yr*0.01 <= 0.0 or fullrange:
        plt.ylim(1.0e-10,ymax+yr*0.01)
    else:
        plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
    if chooserange:
        plt.ylim(y1,y2)

    # plot labels
    plt.xlabel(xlab, {'color' : 'k'})
    try:
        plt.ylabel(ylab1, {'color' : 'k'})
    except:
        ylab1 = '10**%d e-/s' % nrm
        plt.ylabel(ylab1, {'color' : 'k'})

    if plotgrid:
        plt.grid()

    ax.minorticks_on()
    # save plot to file
    if outfile is not None:
        print("Writing output file {}...".format(outfile))
        plt.savefig(outfile)

    if not noninteractive:
        plt.show()
    # end time
    kepmsg.clock('KEPDRAW completed at' , logfile, verbose)


def kepdraw_main():
    import argparse

    parser = argparse.ArgumentParser(
             description='Interactive plotting of Kepler time series data',
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile', default=None,
                        help='name of output PNG file')
    parser.add_argument('--datacol', default=None,
                        help='Name of data column to plot')
    parser.add_argument('--ploterr', action='store_true',
                        help='Plot data error bars?')
    parser.add_argument('--errcol', default='SAP_FLUX_ERR',
                        help='Name of data error column', type=str)
    parser.add_argument('--quality', action='store_true',
                        help=('Ignore cadences where the data quality is'
                              ' questionable?'))
    parser.add_argument('--lcolor', default='#363636',
                        help='HTML color of data line within plot', type=str)
    parser.add_argument('--lwidth', default=1.0,
                        help='type of image intensity scale', type=float)
    parser.add_argument('--fcolor', default='#a8a7a7',
                        help='HTML color to fill the area below the data', type=str)
    parser.add_argument('--falpha', default=0.2,
                        help='type of image intensity scale', type=float)
    parser.add_argument('--labelsize', default=20.,
                        help='Fontsize of axis labels', type=float)
    parser.add_argument('--ticksize', default=14.,
                        help='Fontsize of numeric tick labels', type=float)
    parser.add_argument('--xsize', default=18.,
                        help='X-dimension size of plot', type=float)
    parser.add_argument('--ysize', default=6.,
                        help='Y-dimension size of plot', type=float)
    parser.add_argument('--fullrange', action='store_true',
                        help='Plot flux range from 0.0 e-/sec?')
    parser.add_argument('--chooserange', action='store_true',
                        help='Choose Y-axis range?')
    parser.add_argument('--ymin', default=0.,
                        help='Low limit of the Y-axis range to plot [e-/s]',
                        type=float)
    parser.add_argument('--ymax', default=1e4,
                        help='High limit of the Y-axis range to plot [e-/s]',
                        type=float)
    parser.add_argument('--plotgrid', action='store_true',
                        help='Plot axis grid?')
    parser.add_argument('--ylabel', default='e$^-$ s$^{-1}$',
                        help='Plot axis label', type=str)
    parser.add_argument('--plottype', default='fast', help='plot type',
                        type=str, choices=['fast','pretty'])
    parser.add_argument('--non-interactive', action='store_true',
                        help='Pop up matplotlib plot window?',
                        dest='noninteractive')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepdraw.log', type=str)
    args = parser.parse_args()
    kepdraw(args.infile, args.outfile, args.datacol, args.ploterr, args.errcol,
            args.quality, args.lcolor, args.lwidth, args.fcolor, args.falpha,
            args.labelsize, args.ticksize, args.xsize, args.ysize,
            args.fullrange, args.chooserange, args.ymin, args.ymax,
            args.plotgrid, args.ylabel, args.plottype, args.noninteractive,
            args.verbose, args.logfile)
