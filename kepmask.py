import numpy as np
from matplotlib import pyplot as plt
from copy import copy
import kepio, kepmsg, kepkey, kepplot
import sys, time, re, math, os

# global variables

infile = False; aperfile = False; maskfile = 'mask.txt'
plotfile = 'kepmask.png'; pxdim = 0; pydim = 0; pimg = None; mask = []; zscale = False
xmin = 0.0; xmax = 1000.0; ymin = 0.0; ymax = 1000.0; zmin = False; zmax = False
kepid = ''; ra = ''; dec = ''; kepmag = ''; season = ''; quarter = -1
skygroup = ''; channel = ''; module = ''; output = ''; column = ''; row = ''
colmap='jet'; aid = None; bid = None; cid = None; did = None; eid = None; fid = None
pkepmag = None; pkepid = None; pra = None; pdec = None
cmdLine = False

# -----------------------------------------------------------
# core code

def kepmask(infile,mfile,pfile,tabrow,imin,imax,iscale,cmap,verbose,logfile,status,cLine=False):

    global pimg, zscale, zmin, zmax, xmin, xmax, ymin, ymax, quarter
    global pxdim, pydim, kepmag, skygroup, season, channel
    global module, output, row, column, maskfile, plotfile
    global pkepid, pkepmag, pra, pdec, colmap, cmdLine

# input arguments

    status = 0
    np.seterr(all="ignore")
    zmin = imin; zmax = imax; zscale = iscale; colmap = cmap
    maskfile = mfile; plotfile = pfile
    cmdLine = cLine

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPMASK -- '
    call += 'infile='+infile+' '
    call += 'maskfile='+mfile+' '
    call += 'plotfile='+pfile+' '
    call += 'tabrow='+str(tabrow)+' '
    call += 'imin='+str(imin)+' '
    call += 'imax='+str(imax)+' '
    call += 'iscale='+str(iscale)+' '
    call += 'cmap='+str(cmap)+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPMASK started at',logfile,verbose)

# reference color map

    if cmap == 'browse':
        status = cmap_plot()

# open TPF FITS file and check tabrow exists

    if status == 0:
        tpf, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        try:
            naxis2 = tpf['TARGETTABLES'].header['NAXIS2']
        except:
            txt = 'ERROR -- KEPMASK: No NAXIS2 keyword in ' + infile + '[TARGETTABLES]'
            status = kepmsg.err(logfile,txt,True)
    if status == 0 and tabrow > naxis2:
        txt = 'ERROR -- KEPMASK: tabrow is too large. There are ' + str(naxis2) + ' rows in the table.'
        status = kepmsg.err(logfile,txt,True)
    if status == 0:
        status = kepio.closefits(tpf,logfile,verbose)

# read TPF data pixel image

    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, pixels, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
        img = pixels[tabrow]
        pkepid = copy(kepid)
        pra = copy(ra)
        pdec = copy(dec)
        pkepmag = copy(kepmag)
        pxdim = copy(xdim)
        pydim = copy(ydim)
        pimg = copy(img)

# print target data

    if status == 0:
        print ''
        print '      KepID:  %s' % kepid
        print ' RA (J2000):  %s' % ra
        print 'Dec (J2000):   %s' % dec
        print '     KepMag:   %s' % kepmag
        print '   SkyGroup:   %2s' % skygroup
        print '     Season:   %2s' % str(season)
        print '    Channel:   %2s' % channel
        print '     Module:   %2s' % module
        print '     Output:    %1s' % output
        print ''

# subimage of channel for plot

    if status == 0:
        ymin = copy(row)
        ymax = ymin + ydim
        xmin = copy(column)
        xmax = xmin + xdim

# intensity scale

    if status == 0:
        pimg, imin, imax = kepplot.intScale1D(pimg,zscale)
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

    if status == 0:
        plt.figure(figsize=[10,7])
        plotimage(cmdLine)

# -----------------------------------------------------------
# plot channel image

def plotimage(cmdLine=False):

    global aid, bid, cid, did, eid, fid

# print image and source location data on plot

    plt.ion()
    plt.clf()
    plt.axes([0.73,0.09,0.25,0.4])
    plt.text(0.1,1.0,'      KepID: %s' % pkepid,fontsize=12)
    plt.text(0.1,0.9,' RA (J2000): %s' % pra,fontsize=12)
    plt.text(0.1,0.8,'Dec (J2000): %s' % pdec,fontsize=12)
    plt.text(0.1,0.7,'     KepMag: %s' % pkepmag,fontsize=12)
    plt.text(0.1,0.6,'   SkyGroup: %2s' % skygroup,fontsize=12)
    plt.text(0.1,0.5,'     Season: %2s' % str(season),fontsize=12)
    plt.text(0.1,0.4,'    Channel: %2s' % channel,fontsize=12)
    plt.text(0.1,0.3,'     Module: %2s' % module,fontsize=12)
    plt.text(0.1,0.2,'     Output: %1s' % output,fontsize=12)
    plt.text(0.1,0.1,'     Column: %4s' % column,fontsize=12)
    plt.text(0.1,0.0,'        Row: %4s' % row,fontsize=12)
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.xlim(0.0,1.0)
    plt.ylim(-0.05,1.12)

# clear button

    plt.axes([0.73,0.86,0.25,0.11])
    plt.text(0.5,0.5,'CLEAR',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    aid = plt.connect('button_press_event',clicker1)

# load mask button

    plt.axes([0.73,0.74,0.25,0.11])
    plt.text(0.5,0.5,'LOAD',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    bid = plt.connect('button_press_event',clicker2)

# dump custom aperture to file button

    plt.axes([0.73,0.62,0.25,0.11])
    plt.text(0.5,0.5,'DUMP',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    cid = plt.connect('button_press_event',clicker3)

# print window to png file button

    plt.axes([0.73,0.50,0.25,0.11])
    plt.text(0.5,0.5,'PRINT',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    did = plt.connect('button_press_event',clicker4)

# set the image window location and size

    ax = plt.axes([0.07,0.09,0.63,0.88])

# force tick labels to be absolute rather than relative

    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.subplots_adjust(0.06,0.1,0.93,0.88)
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)

# plot the image window

    imgsum = np.empty((pydim,pxdim))
    n = 0
    for i in range(pydim):
        for j in range(pxdim):
            imgsum[i,j] = pimg[n]
            n += 1
    plt.imshow(imgsum,aspect='auto',interpolation='nearest',origin='lower',
           extent=(xmin,xmax,ymin,ymax),cmap=colmap,vmin=zmin,vmax=zmax)
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
        x = [m-0.5,m+0.5,m+0.5,m-0.5,m-0.5]
        y = [n-0.5,n-0.5,n+0.5,n+0.5,n-0.5]
        plt.fill(x,y,sqcol,alpha=alpha,ec=sqcol)
    fid = plt.connect('key_press_event',clicker6)

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
                plotimage(cmdLine)

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
                    lines, status = kepio.openascii(maskfile,'r',None,False)
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
                        plotimage(cmdLine)
                except:
                    txt = 'ERROR -- KEPMASK: Cannot open or read mask file ' + maskfile
                    kepmsg.err(logfile,txt,True)
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
                if (os.path.isfile(maskfile)):
                    os.remove(maskfile)
                out = open(maskfile,'a')
                out.write(masktxt[:-1]+'\n')
                out.close()
                print 'Wrote custom aperture definition file ' + maskfile
    return

# -----------------------------------------------------------
# print plot to png with left-mouse click

def clicker4(event):

    if event.inaxes:
        if event.button == 1:
            if (event.x > 601 and event.x < 801 and
                event.y > 285 and event.y < 347):
                plt.savefig(plotfile)
                print 'Wrote plot hardcopy file ' + plotfile
    return

# -----------------------------------------------------------
# this function will be called with every click of the mouse

def clicker6(event):

    global mask, aid, bid, cid, did, eid, fid

    if event.inaxes:
        if event.key == 'x':
            if (event.x > 75 and event.x < 580 and
                event.y > 53 and event.y < 550):
                if colmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                              'gray','pink','RdGy']:
                    sqcol = 'g'
                    alpha = 0.5
                else:
                    sqcol = '#ffffee'
                    alpha = 0.8
                m = event.xdata + 0.5
                n = event.ydata + 0.5
                txt = str(int(m))+','+str(int(n))
                if txt in mask:
                    tmpmask = []
                    for pixel in mask:
                        if pixel != txt:
                            tmpmask.append(pixel)
                    mask = tmpmask
                else:
                    mask.append(txt)
                plotimage(cmdLine)

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot():

    plt.figure(figsize=[5,10])
    plt.ion()
    a=outer(ones(10),arange(0,1,0.01))
    plt.subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        plt.subplot(l,1,i+1)
        plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        plt.imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        plt.text(100.85,0.5,m,fontsize=10)

# render plot

    plt.ion()
    plt.show()

    status = 1
    return status

# -----------------------------------------------------------
# main
if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Plot, create or edit custom light curve extraction masks for target pixel files')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='name of input target pixel FITS file', type=str)
    parser.add_argument('maskfile', help='name of ASCII custom aperture definition file', type=str)
    parser.add_argument('--plotfile', default='', help='name of output PNG plot file', type=str)
    parser.add_argument('--tabrow', default=2177, help='The table row containing the image to plot', type=int)
    parser.add_argument('--imin', default=1.5e5, help='minimum of image intensity scale [e-]', type=float)
    parser.add_argument('--imax', default=5.0e5, help='maximum of image intensity scale [e-]', type=float)
    parser.add_argument('--iscale', default='logarithmic', help='type of image intensity scale',
        type=str, choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap', type=str)
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True

    kepmask(args.infile, args.maskfile, args.plotfile, args.tabrow, args.imin, args.imax, args.iscale,
        args.cmap, args.verbose, args.logfile, args.status, cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepmask.par")
    t = iraf.IrafTaskFactory(taskname="kepmask", value=parfile, function=kepmask)
