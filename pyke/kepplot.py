from . import kepmsg, kepstat
import math
import numpy as np
from matplotlib import pyplot as plt

def location(shape):
    """shape the window, enforce absolute scaling, rotate the labels"""

    # position first axes inside the plotting window
    ax = plt.axes(shape)
    # force tick labels to be absolute rather than relative
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    # rotate y labels by 90 deg
    labels = ax.get_yticklabels()

    return ax

def plot1d(x, y, cadence, lcolor, lwidth, fcolor, falpha, underfill):
    """plot a 1d distribution"""

    # pad first and last points in case a fill is required
    x = np.insert(x, [0], [x[0]])
    x = np.append(x, [x[-1]])
    y = np.insert(y, [0], [-1.0e10])
    y = np.append(y, -1.0e10)

    # plot data so that data gaps are not spanned by a line
    ltime = np.array([], dtype='float64')
    ldata = np.array([], dtype='float32')
    for i in range(1, len(x)-1):
        if x[i] - x[i - 1] < 2.0 * cadence / 86400:
            ltime = np.append(ltime, x[i])
            ldata = np.append(ldata, y[i])
        else:
            plt.plot(ltime, ldata, color=lcolor, linestyle='-',
                     linewidth=lwidth)
            ltime = np.array([], dtype='float64')
            ldata = np.array([], dtype='float32')
    plt.plot(ltime, ldata, color=lcolor, linestyle='-', linewidth=lwidth)

    # plot the fill color below data time series, with no data gaps
    if underfill:
        plt.fill(x, y, fc=fcolor, linewidth=0.0, alpha=falpha)

def RangeOfPlot(x, y, pad, origin):
    """determine data limits"""
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    xr = xmax - xmin
    yr = ymax - ymin
    plt.xlim(xmin - xr * pad, xmax + xr * pad)
    plt.ylim(ymin - yr * pad, ymax + yr * pad)
    if origin:
        if ymin - yr * pad <= 0.0:
            plt.ylim(1.0e-10, ymax + yr * pad)
        else:
            plt.ylim(ymin - yr * pad, ymax + yr * pad)

def cleanx(time, logfile, verbose):
    """clean up x-axis of plot"""

    try:
        time0 = float(int(time[0] / 100) * 100.0)
        if time0 < 2.4e6:
            time0 += 2.4e6
        timeout = time - time0
        label = "BJD $-$ {}".format(time0)
    except:
        txt = ("ERROR -- KEPPLOT.CLEANX: cannot calculate plot scaling in "
               "x dimension")
        kepmsg.err(logfile, txt, verbose)

    return timeout, label

def cleany(signal, cadence, logfile, verbose):
    """clean up y-axis of plot"""
    try:
        signal /= cadence
        nrm = math.ceil(math.log10(np.nanmax(signal))) - 1.0
        signal = signal / 10 ** nrm
        if nrm == 0:
            label = 'Flux (e$^-$ s$^{-1}$)'
        else:
            label = "Flux ($10^%d$" % nrm + "e$^-$ s$^{-1}$)"
    except:
        txt = ("ERROR -- KEPPLOT.CLEANY: cannot calculate plot scaling in "
               "y dimension")
        kepmsg.err(logfile, txt, verbose)

    return signal, label


def limits(x, y, logfile, verbose):
    """plot limits"""

    try:
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xr = xmax - xmin
        yr = ymax - ymin
        x = np.insert(x, [0], [x[0]])
        x = np.append(x, [x[-1]])
        y = np.insert(y, [0], [0.0])
        y = np.append(y, 0.0)
    except:
        txt = 'ERROR -- KEPPLOT.LIMITS: cannot calculate plot limits'
        kepmsg.err(logfile, txt, verbose)

    return x, y, xmin,  xmax, xr, ymin, ymax, yr

def labels(xlab, ylab, labcol, fs):
    """plot labels"""

    plt.xlabel(xlab, fontsize=fs, color=labcol)
    plt.ylabel(ylab, fontsize=fs, color=labcol)


def intScale1D(image, imscale):
    """intensity scale limits of 1d array"""

    nstat = 2; work2 = []
    image = np.ma.array(image, mask=np.isnan(image))
    work1 = np.array(np.sort(image), dtype=np.float32)
    for i in range(len(work1)):
        if 'nan' not in str(work1[i]).lower():
            work2.append(work1[i])
    work2 = np.array(work2, dtype=np.float32)
    if int(float(len(work2)) / 10 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 10 + 0.5)
    zmin = np.median(work2[:nstat])
    zmax = np.median(work2[-nstat:])
    if imscale == 'logarithmic':
        if zmin < 0.0:
            zmin = 100.0
        if np.any(image <= 0):
            image = np.log10(image + abs(image.min()) + 1)
        else:
            image = np.log10(image)
        zmin = math.log10(zmin)
        zmax = math.log10(zmax)
    if imscale == 'squareroot':
        if zmin < 0.0:
            zmin = 100.0
        if np.any(image < 0):
            image = np.sqrt(image + abs(image.min()))
        else:
            image = np.sqrt(image)
        zmin = math.sqrt(zmin)
        zmax = math.sqrt(zmax)

    return image, zmin, zmax

def intScale2D(image, imscale):
    """intensity scale limits of 2d array"""
    nstat = 2
    work1 = np.array([], dtype=np.float32)
    (ysiz, xsiz) = np.shape(image)
    for i in range(ysiz):
        for j in range(xsiz):
            if np.isfinite(image[i, j]) and image[i, j] > 0.0:
                work1 = np.append(work1, image[i, j])
    work2 = np.array(np.sort(work1))
    if int(float(len(work2)) / 1000 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 1000 + 0.5)
    zmin = np.median(work2[:nstat])
    zmax = np.median(work2[-nstat:])
    if imscale == 'logarithmic':
        image = np.log10(image)
        zmin = math.log10(zmin)
        zmax = math.log10(zmax)
    if imscale == 'squareroot':
        image = np.sqrt(image)
        zmin = math.sqrt(zmin)
        zmax = math.sqrt(zmax)

    return image, zmin, zmax

def borders(maskimg, xdim, ydim, pixcoord1, pixcoord2, bit, lcolor, lstyle, lwidth):
    """plot mask borders in CCD coordinates"""
    for i in range(1, ydim):
        for j in range(1, xdim):
            if (kepstat.bitInBitmap(maskimg[i, j], bit) and not
                kepstat.bitInBitmap(maskimg[i - 1, j], bit)):
                x = np.array([pixcoord1[j - 1, i], pixcoord1[j, i]]) + 0.5
                y = np.array([pixcoord2[j, i], pixcoord2[j , i]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
            if (not kepstat.bitInBitmap(maskimg[i, j], bit) and
                kepstat.bitInBitmap(maskimg[i - 1, j], bit)):
                x = np.array([pixcoord1[j - 1, i], pixcoord1[j, i]]) + 0.5
                y = np.array([pixcoord2[j, i], pixcoord2[j, i]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
            if (kepstat.bitInBitmap(maskimg[i, j], bit) and not
                kepstat.bitInBitmap(maskimg[i, j - 1], bit)):
                x = np.array([pixcoord1[j, i], pixcoord1[j, i]]) - 0.5
                y = np.array([pixcoord2[j, i - 1], pixcoord2[j, i]]) + 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
            if (not kepstat.bitInBitmap(maskimg[i, j], bit) and
                    kepstat.bitInBitmap(maskimg[i, j - 1], bit)):
                x = np.array([pixcoord1[j, i], pixcoord1[j, i]]) - 0.5
                y = np.array([pixcoord2[j, i - 1],pixcoord2[j, i]]) + 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    # corner cases
    for j in range(ydim):
        try:
            if (kepstat.bitInBitmap(maskimg[j, 0], bit) and not
                kepstat.bitInBitmap(maskimg[j - 1,0], bit)):
                x = np.array([pixcoord1[0, j], pixcoord1[1, j]]) - 0.5
                y = np.array([pixcoord2[0, j], pixcoord2[0, j]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        try:
            if (not kepstat.bitInBitmap(maskimg[j + 1, 0], bit) and
                    kepstat.bitInBitmap(maskimg[j,0],bit)):
                x = np.array([pixcoord1[0, j], pixcoord1[1, j]]) - 0.5
                y = np.array([pixcoord2[0, j], pixcoord2[0, j]]) + 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j, 0], bit):
            x = np.array([pixcoord1[0, j], pixcoord1[0, j]]) - 0.5
            try:
                y = np.array([pixcoord2[0, j], pixcoord2[0, j + 1]]) - 0.5
            except:
                y = np.array([pixcoord2[0, j - 1], pixcoord2[0, j]]) + 0.5
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j, xdim - 1], bit):
            x = np.array([pixcoord1[xdim - 1, j], pixcoord1[xdim - 1, j]]) + 0.5
            try:
                y = (np.array([pixcoord2[xdim - 1, j],
                               pixcoord2[xdim - 1, j + 1]]) - 0.5)
            except:
                y = (np.array([pixcoord2[xdim - 1, j - 1],
                              pixcoord2[xdim - 1, j]]) + 0.5)
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
    for i in range(xdim):
        try:
            if (kepstat.bitInBitmap(maskimg[0, i], bit) and not
                kepstat.bitInBitmap(maskimg[0, i - 1], bit)):
                x = np.array([pixcoord1[i, 0], pixcoord1[i, 0]]) - 0.5
                y = np.array([pixcoord2[i, 0], pixcoord2[i, 1]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)
        except:
            pass
        try:
            if (not kepstat.bitInBitmap(maskimg[0, i + 1], bit) and
                    kepstat.bitInBitmap(maskimg[0, i], bit)):
                x = np.array([pixcoord1[i, 0], pixcoord1[i, 0]]) + 0.5
                y = np.array([pixcoord2[i, 0], pixcoord2[i, 1]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0, i], bit):
            try:
                x = np.array([pixcoord1[i, 0], pixcoord1[i + 1, 0]]) - 0.5
            except:
                x = np.array([pixcoord1[i - 1, 0], pixcoord1[i, 0]]) + 0.5
            y = np.array([pixcoord2[i, 0], pixcoord2[i, 0]]) - 0.5
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim - 1, i], bit):
            try:
                x = (np.array([pixcoord1[i, ydim - 1],
                               pixcoord1[i + 1, ydim - 1]]) - 0.5)
            except:
                x = (np.array([pixcoord1[i - 1, ydim - 1],
                               pixcoord1[i, ydim - 1]]) - 0.5)
            y = np.array([pixcoord2[i, ydim - 1], pixcoord2[i, ydim - 1]]) + 0.5
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[ydim - 1, xdim - 1], bit):
        x = (np.array([pixcoord1[xdim - 2, ydim - 1],
                       pixcoord1[xdim - 1, ydim - 1]]) + 0.5)
        y = (np.array([pixcoord2[xdim - 1, ydim - 1],
                       pixcoord2[xdim - 1, ydim - 1]]) + 0.5)
        plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0, xdim - 1], bit):
        x = np.array([pixcoord1[xdim - 1, 0], pixcoord1[xdim - 1, 0]]) + 0.5
        y = np.array([pixcoord2[xdim - 1, 0], pixcoord2[xdim - 1, 1]]) - 0.5
        plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    return

def PrfBorders(maskimg,xdim,ydim,pixcoord1,pixcoord2,bit,lcolor,lstyle,lwidth):
    """plot mask borders in CCD coordinates"""
    for i in range(1, ydim):
        for j in range(1, xdim):
            if (kepstat.bitInBitmap(maskimg[i, j], bit) and not
                kepstat.bitInBitmap(maskimg[i - 1, j], bit)):
                x = np.array([pixcoord1[j - 1, i], pixcoord1[j, i]]) + 0.5
                y = np.array([pixcoord2[j, i], pixcoord2[j, i]]) - 0.5
                plt.plot(x*50, y*50, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)
            if (not kepstat.bitInBitmap(maskimg[i, j], bit) and
                    kepstat.bitInBitmap(maskimg[i - 1, j], bit)):
                x = np.array([pixcoord1[j - 1, i], pixcoord1[j, i]]) + 0.5
                y = np.array([pixcoord2[j , i], pixcoord2[j, i]]) - 0.5
                plt.plot(x*50, y*50, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)
            if (kepstat.bitInBitmap(maskimg[i, j], bit) and not
                kepstat.bitInBitmap(maskimg[i, j - 1], bit)):
                x = np.array([pixcoord1[j, i], pixcoord1[j, i]]) - 0.5
                y = np.array([pixcoord2[j, i - 1], pixcoord2[j, i]]) + 0.5
                plt.plot(x*50, y*50, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)
            if (not kepstat.bitInBitmap(maskimg[i, j], bit) and
                    kepstat.bitInBitmap(maskimg[i, j - 1], bit)):
                x = np.array([pixcoord1[j, i], pixcoord1[j, i]]) - 0.5
                y = np.array([pixcoord2[j, i - 1], pixcoord2[j, i]]) + 0.5
                plt.plot(x*50, y*50, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)

    # corner cases
    for j in range(ydim):
        try:
            if (kepstat.bitInBitmap(maskimg[j, 0], bit) and not
                kepstat.bitInBitmap(maskimg[j - 1, 0], bit)):
                x = np.array([pixcoord1[0, j], pixcoord1[1, j]]) - 0.5
                y = np.array([pixcoord2[0, j], pixcoord2[0, j]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        try:
            if (not kepstat.bitInBitmap(maskimg[j + 1, 0], bit) and
                    kepstat.bitInBitmap(maskimg[j, 0], bit)):
                x = np.array([pixcoord1[0, j], pixcoord1[1, j]]) - 0.5
                y = np.array([pixcoord2[0, j], pixcoord2[0, j]]) + 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j, 0], bit):
            x = np.array([pixcoord1[0,j],pixcoord1[0,j]]) - 0.5
            try:
                y = np.array([pixcoord2[0, j], pixcoord2[0, j + 1]]) - 0.5
            except:
                y = np.array([pixcoord2[0, j - 1], pixcoord2[0, j]]) + 0.5
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j, xdim - 1], bit):
            x = np.array([pixcoord1[xdim - 1, j], pixcoord1[xdim - 1, j]]) + 0.5
            try:
                y = (np.array([pixcoord2[xdim - 1, j],
                               pixcoord2[xdim - 1, j + 1]]) - 0.5)
            except:
                y = (np.array([pixcoord2[xdim - 1, j - 1],
                               pixcoord2[xdim - 1, j]]) + 0.5)
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
    for i in range(xdim):
        try:
            if (kepstat.bitInBitmap(maskimg[0, i], bit) and not
                kepstat.bitInBitmap(maskimg[0, i - 1], bit)):
                x = np.array([pixcoord1[i, 0], pixcoord1[i, 0]]) - 0.5
                y = np.array([pixcoord2[i, 0], pixcoord2[i, 1]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        except:
            pass
        try:
            if (not kepstat.bitInBitmap(maskimg[0, i + 1], bit) and
                    kepstat.bitInBitmap(maskimg[0, i], bit)):
                x = np.array([pixcoord1[i, 0], pixcoord1[i, 0]]) + 0.5
                y = np.array([pixcoord2[i, 0], pixcoord2[i, 1]]) - 0.5
                plt.plot(x, y, color=lcolor, linestyle=lstyle,
                         linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0, i], bit):
            try:
                x = np.array([pixcoord1[i, 0], pixcoord1[i + 1, 0]]) - 0.5
            except:
                x = np.array([pixcoord1[i - 1, 0], pixcoord1[i, 0]]) + 0.5
            y = np.array([pixcoord2[i,0],pixcoord2[i,0]]) - 0.5
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim - 1, i], bit):
            try:
                x = (np.array([pixcoord1[i, ydim - 1],
                               pixcoord1[i + 1, ydim-1]]) - 0.5)
            except:
                x = (np.array([pixcoord1[i - 1, ydim - 1],
                               pixcoord1[i, ydim - 1]]) - 0.5)
            y = (np.array([pixcoord2[i, ydim - 1],
                           pixcoord2[i, ydim - 1]]) + 0.5)
            plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[ydim - 1, xdim -1], bit):
        x = (np.array([pixcoord1[xdim - 2, ydim - 1],
                       pixcoord1[xdim - 1, ydim - 1]]) + 0.5)
        y = (np.array([pixcoord2[xdim - 1, ydim - 1],
                       pixcoord2[xdim - 1, ydim - 1]]) + 0.5)
        plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0, xdim - 1], bit):
        x = np.array([pixcoord1[xdim - 1, 0], pixcoord1[xdim - 1, 0]]) + 0.5
        y = np.array([pixcoord2[xdim - 1, 0], pixcoord2[xdim - 1, 1]]) - 0.5
        plt.plot(x, y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
