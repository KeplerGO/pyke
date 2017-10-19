from .utils import PyKEArgumentHelpFormatter
from . import kepio, kepmsg, kepkey
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt

# global variables

instr = ''; cadence = 1800.0; barytime0 = 0
nrm = 1; barytime = []; flux = []; xmin = 0; xmax = 1
ymin = 0; ymax = 1; xr = 1; yr = 1; xlab = ''; ylab = ''
mask = []; aid = None; bid = None; cid = None; did = None; eid = None; fid = None
clobb = True; outf = ''; verb = True; logf = ''; rinf = ''


__all__ = ['keprange']


def keprange(infile, outfile=None, datacol='SAP_FLUX', rinfile='',
             overwrite=False, verbose=False, logfile='keprange.log'):
    """
    keprange -- Define time ranges interactively for use with other PyKE tasks.

    A number of PyKE tasks, e.g. kepdetrend, kepoutlier, require the user to
    specify ranges in time over which to operate. keprange provides a visual
    and interactive tool with which to define time ranges and store them within
    an ASCII file. Choices are made using a GUI. Use the left-button of your
    mouse to select ranges. An existing ASCII file can be loaded, a new ASCII
    file can be written, the list of times can be cleared or printed using the
    buttons on the GUI.

    Parameters
    ----------
    infile : str
        The name of a MAST standard format FITS file containing a Kepler light
        curve within the first data extension.
    outfile : str
        The name of the output ASCII file storing time ranges for future use in
        other PyKE tools.
    rinfile : str
        An existing ASCII file containing time ranges in Barycentric Julian
        Date (BJD) can be uploaded into the task. This can be used as a basis
        for a new set of time ranges. This argument is optional and is not
        prompted for automatically. If no ascii file will be input then
        rinfile=None will clear the argument buffer after a previous use.
    datacol : str
        The datacol name containing data stored within extension 1 of infile.
        This data will be plotted against time so that the user can choose
        appropriate time ranges.
    overwrite : bool
        Overwrite the output file?
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.

    Examples
    --------
    .. code-block:: bash

        $ keprange kplr002436324-2009259160929_llc.fits --verbose

    .. image:: ../_static/images/api/keprange.png
        :align: center
    """

    if outfile is None:
        outfile = infile.split('.')[0] + "-{}.txt".format(__all__[0])

    # startup parameters
    global instr, cadence, barytime0, nrm, barytime, flux
    global xmin, xmax, ymin, ymax, xr, yr, xlab, ylab
    global clobb, outf, verb, logf, rinf, col, bjdref, cade

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPRANGE -- '
            + ' infile={}'.format(infile)
            + ' outfile={}'.format(outfile)
            + ' rinfile={}'.format(rinfile)
            + ' datacol={}'.format(datacol)
            + ' overwrite={}'.format(overwrite)
            + ' verbose={}'.format(verbose)
            + ' logfile={}'.format(logfile))
    kepmsg.log(logfile, call+'\n', verbose)
    clobb = overwrite
    outf = outfile
    verb = verbose
    logf = logfile
    rinf = rinfile

    # start time
    kepmsg.clock('KEPRANGE started at: ', logfile, verbose)

    # overwrite output file
    if overwrite:
        kepio.overwrite(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        errmsg = 'ERROR -- KEPRANGE: {} exists. Use --overwrite'.format(outfile)
        kepmsg.err(logfile, errmsg, verbose)

    ## open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence
    cade = cadenom

    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # input data
    table = instr[1].data

    # filter out NaNs
    work1 = []; work2 = []
    col = datacol
    barytime = kepio.readtimecol(infile, table, logfile, verbose)
    try:
        flux = instr[1].data.field(col)
    except:
        errmsg = ('ERROR -- KEPRANGE: no datacol named {} in table {} [1]'
                  .format(col, infile))
        kepmsg.err(infile, errmsg, verbose)
    for i in range(len(barytime)):
        if (np.isfinite(barytime[i]) and np.isfinite(flux[i])
            and flux[i] != 0.0):
            work1.append(barytime[i] + bjdref)
            work2.append(flux[i])
    barytime = np.array(work1, dtype=np.float64)
    flux = np.array(work2, dtype=np.float32) / cadenom

    # clean up x-axis unit
    barytime0 = float(int(tstart / 100) * 100.0)
    barytime = barytime - barytime0
    xlab = 'BJD $-$ {}'.format(barytime0)

    # clean up y-axis units
    nrm = len(str(int(flux.max()))) - 1
    flux = flux / 10 ** nrm
    ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

    # data limits
    xmin = barytime.min()
    xmax = barytime.max()
    ymin = flux.min()
    ymax = flux.max()
    xr = xmax - xmin
    yr = ymax - ymin
    flux[0] = 0.0
    flux[-1] = 0.0

    # plot new light curve
    plt.rcParams['figure.dpi'] = 80
    plt.figure(figsize=[17, 7])
    plotlc()


def plotlc():
    global aid, bid, cid, did, eid, fid, mask

    # load button
    plt.clf()
    plt.axes([0.06, 0.02, 0.22, 0.1])
    plt.text(0.5, 0.5, 'LOAD', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    aid = plt.connect('button_press_event', clicker1)

    # save button
    plt.axes([0.2933, 0.02, 0.22, 0.1])
    plt.text(0.5, 0.5, 'SAVE', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    bid = plt.connect('button_press_event', clicker2)

    # clear button
    plt.axes([0.5266, 0.02, 0.22, 0.1])
    plt.text(0.5, 0.5, 'CLEAR', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    cid = plt.connect('button_press_event', clicker3)

    # print button
    plt.axes([0.76, 0.02, 0.22, 0.1])
    plt.text(0.5, 0.5, 'PRINT', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    did = plt.connect('button_press_event', clicker4)

    # light curve
    plt.axes([0.06, 0.213, 0.92, 0.77])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ltime = []; ldata = []
    work1 = instr[1].data.field(0) + bjdref
    work2 = instr[1].data.field(col) / cade
    for i in range(len(work1)):
        if np.isfinite(work1[i]) or np.isfinite(work2[i]):
            ltime.append(work1[i])
            ldata.append(work2[i])
        else:
            ltime = np.array(ltime, dtype=np.float64) - barytime0
            ldata = np.array(ldata, dtype=np.float64) / 10 ** nrm
            plt.plot(ltime,ldata,color='#0000ff',linestyle='-',linewidth=1.0)
            ltime = []; ldata = []
    ltime = np.array(ltime, dtype=np.float64) - barytime0
    ldata = np.array(ldata, dtype=np.float64) / 10 ** nrm
    plt.plot(ltime, ldata, color='#0000ff', linestyle='-', linewidth=1.0)
    plt.fill(barytime, flux, fc='#ffff00', linewidth=0.0, alpha=0.2)
    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab, {'color' : 'k'})
    plt.grid()

    # plt.plot masks
    for i in range(len(mask)):
        t = float(mask[i])
        plt.plot([t, t], [ymin, ymax], color='g', linestyle='-', linewidth=0.5)
    nt = 0
    for i in range(int(len(mask) / 2)):
        t1 = float(mask[nt])
        t2 = float(mask[nt + 1])
        nt += 2
        plt.fill([t1, t1, t2, t2, t1], [ymin, ymax, ymax, ymin, ymin],
                 fc='g', linewidth=0.0, alpha=0.5)
    # plot ranges
    plt.xlim(xmin - xr * 0.01, xmax + xr * 0.01)
    if ymin - yr * 0.01 <= 0.0:
        plt.ylim(1.0e-10, ymax + yr * 0.01)
    else:
        plt.ylim(ymin - yr * 0.01, ymax + yr * 0.01)

    # ranges
    eid = plt.connect('button_press_event', clicker6)

    # render plot
    plt.show()

def clicker1(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 83 and event.x < 383 and
                event.y > 12 and event.y < 68):
                if kepio.fileexists(rinf):
                    mask = []
                    lines = kepio.openascii(rinf, 'r', logf, verb)
                    for line in lines:
                        line = line.strip().split(',')
                        try:
                            float(line[0])
                            float(line[1])
                            if barytime0 > 2.4e6:
                                mask.append(float(line[0]) - barytime0)
                                mask.append(float(line[1]) - barytime0)
                            else:
                                mask.append(float(line[0]) - barytime0 - 2.4e6)
                                mask.append(float(line[1]) - barytime0 - 2.4e6)
                        except:
                            message = 'ERROR -- KEPRANGE: ascii format of ranges '
                            message += 'file not recognized.'
                            kepmsg.err(logf, message, False)
                    plt.disconnect(aid)
                    plt.disconnect(bid)
                    plt.disconnect(cid)
                    plt.disconnect(did)
                    plt.disconnect(eid)
                    plt.disconnect(fid)
                    plotlc()
                else:
                    print('WARNING -- KEPRANGE: input ranges file does not'
                          ' exist or was not provided')
    return

# -----------------------------------------------------------
# save mask to ascii file

def clicker2(event):

    global mask, aid, bid, cid, did, eid, fid, clobb

    if clobb:
        kepio.overwrite(outf, logf, verb)
    if kepio.fileexists(outf):
        message = 'ERROR -- KEPRANGE: {} exists. Use --overwrite'.format(outf)
        kepmsg.err(logf, message, verb)
    else:
        if event.inaxes:
            if event.button == 1:
                if (event.x > 402 and event.x < 702 and
                    event.y > 12 and event.y < 68):
                    nt = 0; txt = ''
                    for i in range(int(len(mask)/2)):
                        t1 = float(mask[nt]) + barytime0
                        t2 = float(mask[nt+1]) + barytime0
                        if t1 < 2.4e6: t1 += 2.4e6
                        if t2 < 2.4e6: t2 += 2.4e6
                        txt += str(t1) + ',' + str(t2) + '\n'
                        nt += 2
                    txt = txt.strip()
                    kepmsg.log(outf, txt, True)
                    print('\nWrote ASCII file ' + outf)
                    plotlc()

# -----------------------------------------------------------
# clear time domain mask

def clicker3(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 723 and event.x < 1022 and
                event.y > 12 and event.y < 68):
                mask = []
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                plotlc()

# -----------------------------------------------------------
# print time domain mask

def clicker4(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 1042 and event.x < 1342 and
                event.y > 12 and event.y < 68):
                nt = 0; txt = ''
                for i in range(int(len(mask)/2)):
                    t1 = float(mask[nt]) + barytime0
                    t2 = float(mask[nt+1]) + barytime0
                    if t1 < 2.4e6: t1 += 2.4e6
                    if t2 < 2.4e6: t2 += 2.4e6
                    txt += str(t1) + ',' + str(t2) + '\n'
                    nt += 2
                txt = txt.strip()
                print('\n' + txt)
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                plotlc()

# -----------------------------------------------------------
# left-click create time ranges

def clicker6(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 83 and event.x < 1337 and
                event.y > 122 and event.y < 558):
                if len(mask) % 2 == 0:
                    mask.append(event.xdata)
                else:
                    if event.xdata > mask[-1]:
                        mask.append(event.xdata)
                    else:
                        mask.append(mask[-1])
                        mask[-2] = event.xdata
                plotlc()

def keprange_main():
    import argparse
    parser = argparse.ArgumentParser(
             description=('Interactively define and store time ranges via a'
                          ' GUI'),
             formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--outfile',
                        help=('Name of ASCII file to output time ranges.'
                              ' If None, outfile is infile-keprange.'),
                        default=None)
    parser.add_argument('--datacol', default='SAP_FLUX',
                        help='Name of diagnostic FITS column', type=str)
    parser.add_argument('--rinfile', default='',
                        help='Name of input ASCII time ranges file')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='keprange.log', dest='logfile', type=str)
    args = parser.parse_args()
    keprange(args.infile, args.outfile, args.datacol, args.rinfile,
             args.overwrite, args.verbose, args.logfile)
