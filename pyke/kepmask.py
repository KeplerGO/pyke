from . import kepio
from . import kepmsg
from . import kepplot
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from copy import copy


__all__ = ['kepmask']


""" global variables """

infile = False; aperfile = False; maskfile = 'mask.txt'
plotfile = 'kepmask.png'; pxdim = 0; pydim = 0; pimg = None; mask = []
zscale = False; xmin = 0.0; xmax = 1000.0; ymin = 0.0; ymax = 1000.0
zmin = False; zmax = False; kepid = ''; ra = ''; dec = ''; kepmag = ''
season = ''; quarter = -1; skygroup = ''; channel = ''; module = ''
output = ''; column = ''; row = ''; colmap='jet'; aid = None; bid = None
cid = None; did = None; eid = None; fid = None; pkepmag = None; pkepid = None
pra = None; pdec = None

def kepmask(infile, nframe, maskfile='mask.txt', plotfile='kepmask.png',
            imin=False, imax=False, iscale=False, cmap='jet', verbose=True,
            logfile='kepmask.log'):
    """
    kepmask plots, creates or edits custom target masks for target pixel
    files. The product from this task is a target mask definition file which
    can be used by kepextract to extract a light curve from target pixel data.
    This tool is a GUI interface for defining a pixel mask by moving a mouse
    over image pixels and selecting them by pressing the 'x' key on the
    keyboard.

    Parameters
    ----------
    infile : str
        The name of a target pixel file from the MAST Kepler archive,
        containing a standard mask definition image in the second data
        extension.
    nframe : int
    maskfile : str
        The name of an ASCII mask definition file. This is either the name of
        a file to be plotted, a file to be created, or a file to be edited.
    plotfile : str
        The name of a PNG plot file containing a record of the mask defined or
        uploaded by this task.
    imin : float (optional)
        Minimum intensity value (in electrons per cadence) for the image
        display. The default minimum intensity level is the median of the
        faintest 10% of pixels in the image.
    imax : float (optional)
        Maximum intensity value (in electrons per cadence) for the image
        display. The default maximum intensity level is the median of the
        brightest 10% of pixels in the image.
    iscale : str (optional)
        Type of intensity scaling for the image display.
            * linear
            * logarithmic
            * squareroot
    cmap : str (optional)
        Color intensity scheme for the image display.
    verbose : bool
        Print informative messages and warnings to the shell and logfile?
    logfile : str
        Name of the logfile containing error and warning messages.
    """

    global pimg, zscale, zmin, zmax, xmin, xmax, ymin, ymax, quarter
    global pxdim, pydim, kepmag, skygroup, season, channel
    global module, output, row, column, mfile, pfile
    global pkepid, pkepmag, pra, pdec, colmap

# input arguments

    zmin = imin; zmax = imax; zscale = iscale; colmap = cmap
    mfile = maskfile; pfile = plotfile

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile, hashline, verbose)
    call = ('KEPMASK -- '
            'infile=' + infile + ' ' + 'maskfile =' + mfile + ' '
            'plotfile =' + pfile + ' ' + 'nframe =' + str(nframe) + ' '
            'imin =' + str(imin) +' ' + 'imax =' + str(imax) + ' '
            'iscale =' + str(iscale)+' ' + 'cmap =' + str(cmap) + ' '
            'verbose=' + str(verbose) + ' ' + 'logfile=' + logfile)

    kepmsg.log(logfile, call + '\n', verbose)
    kepmsg.clock('KEPMASK started at', logfile, verbose)

# open TPF FITS file and check whether or not nframe exists
    try:
        tpf = pyfits.open(infile, mode='readonly')
    except:
        errmsg = ('ERROR -- KEPIO.OPENFITS: cannot open ' +
                  infile + ' as a FITS file')
        kepmsg.err(logfile, errmsg, verbose)
        raise

    try:
        naxis2 = tpf['TARGETTABLES'].header['NAXIS2']
    except:
        errmsg = ('ERROR -- KEPMASK: No NAXIS2 keyword in ' + infile +
                  '[TARGETTABLES]')
        kepmsg.err(logfile, errmsg, verbose)
        raise

    if nframe > naxis2:
        errmsg = ('ERROR -- KEPMASK: nframe is too large. There are ' +
                  str(naxis2) + ' rows in the table.')
        kepmsg.err(logfile, errmsg, verbose)

    kepio.closefits(tpf, logfile, verbose)

# read TPF data pixel image

    kepid, channel, skygroup, module, output, quarter, season, \
    ra, dec, column, row, kepmag, xdim, ydim, pixels = \
        kepio.readTPF(infile,'FLUX',logfile,verbose)
    img = pixels[nframe]
    pkepid = copy(kepid)
    pra = copy(ra)
    pdec = copy(dec)
    pkepmag = copy(kepmag)
    pxdim = copy(xdim)
    pydim = copy(ydim)
    pimg = copy(img)

# print target data

    print('')
    print('      KepID:  {}'.format(kepid))
    print(' RA (J2000):  {}'.format(ra))
    print('Dec (J2000):  {}'.format(dec))
    print('     KepMag:  {}'.format(kepmag))
    print('   SkyGroup:  {}'.format(skygroup))
    print('     Season:  {}'.format(str(season)))
    print('    Channel:  {}'.format(channel))
    print('     Module:  {}'.format(module))
    print('     Output:  {}'.format(output))
    print('')

# subimage of channel for plot

    ymin = copy(row)
    ymax = ymin + ydim
    xmin = copy(column)
    xmax = xmin + xdim

# intensity scale

    pimg, imin, imax = kepplot.intScale1D(pimg, zscale)
    if zmin and zmax and 'log' in zscale:
        zmin = math.log10(zmin)
        zmax = math.log10(zmax)
    elif zmin and zmax and 'sq' in zscale:
        zmin = math.sqrt(zmin)
        zmax = math.sqrt(zmax)
    elif zmin and zmax and 'li' in zscale:
        zmin *= 1.0
        zmax *= 1.0
    else:
        zmin = copy(imin)
        zmax = copy(imax)

# plot limits

    ymin = float(ymin) - 0.5
    ymax = float(ymax) - 0.5
    xmin = float(xmin) - 0.5
    xmax = float(xmax) - 0.5

# plot style

    plt.rcParams['figure.dpi'] = 80
    plt.figure(figsize=[10, 7])
    plotimage()

# plot channel image

def plotimage():

    global aid, bid, cid, did, eid, fid

# print image and source location data on plot

    plt.ion()
    plt.clf()
    plt.axes([0.73, 0.09, 0.25, 0.4])
    plt.text(0.1, 1.0,'      KepID: {}'.format(pkepid, fontsize=12))
    plt.text(0.1, 0.9,' RA (J2000): {}'.format(pra, fontsize=12))
    plt.text(0.1, 0.8,'Dec (J2000): {}'.format(pdec, fontsize=12))
    plt.text(0.1, 0.7,'     KepMag: {}'.format(pkepmag, fontsize=12))
    plt.text(0.1, 0.6,'   SkyGroup: {}'.format(skygroup, fontsize=12))
    plt.text(0.1, 0.5,'     Season: {}'.format(str(season), fontsize=12))
    plt.text(0.1, 0.4,'    Channel: {}'.format(channel, fontsize=12))
    plt.text(0.1, 0.3,'     Module: {}'.format(module, fontsize=12))
    plt.text(0.1, 0.2,'     Output: {}'.format(output, fontsize=12))
    plt.text(0.1, 0.1,'     Column: {}'.format(column, fontsize=12))
    plt.text(0.1, 0.0,'        Row: {}'.format(row, fontsize=12))
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.xlim(0.0, 1.0)
    plt.ylim(-0.05, 1.12)

# clear button

    plt.axes([0.73, 0.86, 0.25, 0.11])
    plt.text(0.5, 0.5, 'CLEAR', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[],
             yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    aid = plt.connect('button_press_event', clicker1)

# load mask button

    plt.axes([0.73, 0.74, 0.25, 0.11])
    plt.text(0.5, 0.5, 'LOAD', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    bid = plt.connect('button_press_event', clicker2)

# dump custom aperture to file button

    plt.axes([0.73, 0.62, 0.25, 0.11])
    plt.text(0.5, 0.5, 'DUMP', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    cid = plt.connect('button_press_event', clicker3)

# print window to png file button

    plt.axes([0.73, 0.50, 0.25, 0.11])
    plt.text(0.5, 0.5, 'PRINT', fontsize=24, weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[], yticks=[])
    plt.fill([0.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0], '#ffffee')
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    did = plt.connect('button_press_event', clicker4)

# set the image window location and size

    ax = plt.axes([0.07, 0.09, 0.63, 0.88])

# force tick labels to be absolute rather than relative

    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.subplots_adjust(0.06, 0.1, 0.93, 0.88)
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)

# plot the image window
    imgsum = pimg.reshape((pydim, pxdim))
    plt.imshow(imgsum, aspect='auto', interpolation='nearest', origin='lower',
               extent=(xmin, xmax, ymin, ymax), cmap=colmap, vmin=zmin,
               vmax=zmax)
    plt.gca().set_autoscale_on(False)
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    plt.ylabel('Pixel Row Number', {'color' : 'k'})

# plot the mask

    if colmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                'gray','pink','RdGy']:
        sqcol = 'g'
        alpha = 0.5
    else:
        sqcol = '#ffffee'
        alpha = 0.8
    for pixel in mask:
        m = int(pixel.split(',')[0])
        n = int(pixel.split(',')[1])
        x = [m - 0.5, m + 0.5, m + 0.5, m - 0.5, m - 0.5]
        y = [n - 0.5, n - 0.5, n + 0.5, n + 0.5, n - 0.5]
        plt.fill(x, y, sqcol, alpha=alpha, ec=sqcol)
    fid = plt.connect('key_press_event', clicker6)

# render plot
    plt.ion()
    plt.show()
# -----------------------------------------------------------
# clear all pixels from pixel mask

def clicker1(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 601 and event.x < 801 and
                event.y > 492 and event.y < 522):
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                mask = []
                plt.clf()
                plotimage()
    return

# -----------------------------------------------------------
# load mask from file

def clicker2(event):

    global mask, aid, bid, cid, did, eid, fid, done

    if event.inaxes:
        if event.button == 1:
            if (event.x > 601 and event.x < 801 and
                event.y > 422 and event.y < 482):
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                try:
                    lines = kepio.openascii(mfile, 'r', None, False)
                    for line in lines:
                        mask = []
                        work = line.strip().split('|')
                        y0 = int(work[3])
                        x0 = int(work[4])
                        work = work[5].split(';')
                        for i in range(len(work)):
                            y = int(work[i].split(',')[0]) + y0
                            x = int(work[i].split(',')[1]) + x0
                            mask.append(str(x) + ',' + str(y))
                        plt.clf()
                        plotimage()
                except:
                    errmsg = ('ERROR -- KEPMASK: Cannot open or read mask '
                              'file ' + mfile)
                    kepmsg.err('kepmask.log', errmsg, True)
    return

# -----------------------------------------------------------
# dump custom aperture definition file

def clicker3(event):

    global aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 601 and event.x < 801 and
                event.y > 354 and event.y < 415):
                masktxt  = 'NEW|'
                masktxt += skygroup + '|'
                masktxt += str(pkepid)
                masktxt += ',TAD_NO_HALO,TAD_NO_UNDERSHOOT_COLUMN|'
                masktxt += str(int(row)) + '|'
                masktxt += str(int(column)) + '|'
                for coord in sorted(set(mask)):
                    masktxt += str(int(coord.split(',')[1]) - int(row)) + ','
                    masktxt += str(int(coord.split(',')[0]) - int(column)) + ';'
                if os.path.isfile(mfile):
                    os.remove(mfile)
                out = open(mfile, 'a')
                out.write(masktxt[:-1] + '\n')
                out.close()
                print('Wrote custom aperture definition file ' + mfile)
    return

# -----------------------------------------------------------
# print plot to png with left-mouse click

def clicker4(event):

    if event.inaxes:
        if event.button == 1:
            if (event.x > 601 and event.x < 801 and
                event.y > 285 and event.y < 347):
                plt.savefig(pfile)
                print('Wrote plot hardcopy file ' + pfile)
    return

# -----------------------------------------------------------
# this function will be called with every click of the mouse

def clicker6(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.key == 'x':
            if (event.x > 75 and event.x < 580 and
                event.y > 53 and event.y < 550):
                if colmap in ['Greys', 'binary', 'bone', 'gist_gray',
                              'gist_yarg', 'gray', 'pink', 'RdGy']:
                    sqcol = 'g'
                    alpha = 0.5
                else:
                    sqcol = '#ffffee'
                    alpha = 0.8
                m = event.xdata + 0.5
                n = event.ydata + 0.5
                txt = str(int(m)) + ',' + str(int(n))
                if txt in mask:
                    tmpmask = []
                    for pixel in mask:
                        if pixel != txt:
                            tmpmask.append(pixel)
                    mask = tmpmask
                else:
                    mask.append(txt)
                plotimage()

def kepmask_main():
    import argparse
    parser = argparse.ArgumentParser(
                   description=("Plot, create or edit custom light curve "
                                "extraction masks for target pixel files ")
                   )
    parser.add_argument('infile', help='name of input target pixel FITS file',
                        type=str)
    parser.add_argument('maskfile', default='maskfile.txt',
                         help='name of ASCII custom aperture definition file',
                         type=str)
    parser.add_argument('--plotfile', default='',
                        help='name of output PNG plot file', type=str)
    parser.add_argument('--nframe', default=2177,
                        help='The number of the frame to plot', type=int)
    parser.add_argument('--imin', default=1.5e5,
                        help='minimum of image intensity scale [e-]',
                        type=float)
    parser.add_argument('--imax', default=5.0e5,
                        help='maximum of image intensity scale [e-]',
                        type=float)
    parser.add_argument('--iscale', default='logarithmic',
                        help='type of image intensity scale',
                        type=str,
                        choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap',
                        type=str)
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepmask.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    kepmask(args.infile, args.maskfile, args.plotfile, args.nframe, args.imin,
            args.imax, args.iscale, args.cmap, args.verbose, args.logfile)
