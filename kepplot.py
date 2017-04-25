#!/usr/bin/env python
import kepmsg, kepstat
import numpy, scipy
import scipy.stats
from numpy import *
from matplotlib import pyplot as plt

# -----------------------------------------------------------
# shape the window, enforce absolute scaling, rotate the labels

def location(shape):

# position first axes inside the plotting window

    ax = plt.axes(shape)
# force tick labels to be absolute rather than relative
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
# rotate y labels by 90 deg
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)

    return ax

# -----------------------------------------------------------
# plot a 1d distribution

def plot1d(x,y,cadence,lcolor,lwidth,fcolor,falpha,underfill):

# pad first and last points in case a fill is required

    x = insert(x,[0],[x[0]])
    x = append(x,[x[-1]])
    y = insert(y,[0],[-1.0e10])
    y = append(y,-1.0e10)

# plot data so that data gaps are not spanned by a line

    ltime = array([],dtype='float64')
    ldata = array([],dtype='float32')
    for i in range(1,len(x)-1):
        if (x[i] - x[i-1]) < 2.0 * cadence / 86400:
            ltime = append(ltime,x[i])
            ldata = append(ldata,y[i])
        else:
            plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)
            ltime = array([],dtype='float64')
            ldata = array([],dtype='float32')
    plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)

# plot the fill color below data time series, with no data gaps

    if underfill:
        plt.fill(x,y,fc=fcolor,linewidth=0.0,alpha=falpha)

        return

# -----------------------------------------------------------
# determine data limits

def RangeOfPlot(x,y,pad,origin):

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

    return

# -----------------------------------------------------------
# clean up x-axis of plot

def cleanx(time,logfile,verbose):

    status = 0
    try:
        time0 = float(int(time[0] / 100) * 100.0)
        if time0 < 2.4e6: time0 += 2.4e6
        timeout = time - time0
        label = 'BJD $-$ %d' % time0
    except:
        txt = 'ERROR -- KEPPLOT.CLEANX: cannot calculate plot scaling in x dimension'
        status = kepmsg.err(logfile,txt,verbose)
        label = ''

    return timeout, label, status

# -----------------------------------------------------------
# clean up y-axis of plot

def cleany(signal,cadence,logfile,verbose):

    status = 0
    try:
        signal /= cadence
        nrm = math.ceil(math.log10(numpy.nanmax(signal))) - 1.0
	signal = signal / 10**nrm
        if nrm == 0:
            label = 'Flux (e$^-$ s$^{-1}$)'
        else:
            label = 'Flux (10$^{%d}$ e$^-$ s$^{-1}$)' % nrm
    except:
        txt = 'ERROR -- KEPPLOT.CLEANY: cannot calculate plot scaling in y dimension'
        status = kepmsg.err(logfile,txt,verbose)
        label = ''

    return signal, label, status

# -----------------------------------------------------------
# plot limits

def limits(x,y,logfile,verbose):

    status = 0
    try:
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xr = xmax - xmin
        yr = ymax - ymin
        x = insert(x,[0],[x[0]])
        x = append(x,[x[-1]])
        y = insert(y,[0],[0.0])
        y = append(y,0.0)
    except:
        txt = 'ERROR -- KEPPLOT.LIMITS: cannot calculate plot limits'
        status = kepmsg.err(logfile,txt,verbose)

    return x, y, xmin,  xmax, xr, ymin, ymax, yr, status

# -----------------------------------------------------------
# plot labels

def labels(xlab,ylab,labcol,fs):

    plt.xlabel(xlab, fontsize=fs, color=labcol)
    plt.ylabel(ylab, fontsize=fs, color=labcol)

    return

# -----------------------------------------------------------
# intensity scale limits of 1d array

def intScale1D(image,imscale):

    seterr(all="ignore")
    nstat = 2; work2 = []
    work1 = array(sort(image),dtype=float32)
    for i in range(len(work1)):
        if 'nan' not in str(work1[i]).lower():
            work2.append(work1[i])
    work2 = array(work2,dtype=float32)
    if int(float(len(work2)) / 10 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 10 + 0.5)
    zmin = median(work2[:nstat])
    zmax = median(work2[-nstat:])
    if imscale == 'logarithmic':
        if zmin < 0.0: zmin = 100.0
        image = log10(image)
        zmin = log10(zmin)
        zmax = log10(zmax)
    if (imscale == 'squareroot'):
        if zmin < 0.0: zmin = 100.0
        image = sqrt(image)
        zmin = sqrt(zmin)
        zmax = sqrt(zmax)

    return image, zmin, zmax

# -----------------------------------------------------------
# intensity scale limits of 2d array

def intScale2D(image,imscale):

    seterr(all="ignore")
    nstat = 2
    work1 = numpy.array([],dtype='float32')
    (ysiz,xsiz) = numpy.shape(image)
    for i in range(ysiz):
        for j in range(xsiz):
            if numpy.isfinite(image[i,j]) and image[i,j] > 0.0:
                work1 = numpy.append(work1,image[i,j])
    work2 = array(sort(work1))
    if int(float(len(work2)) / 1000 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 1000 + 0.5)
    zmin = median(work2[:nstat])
    zmax = median(work2[-nstat:])
    if imscale == 'logarithmic':
        image = log10(image)
        zmin = log10(zmin)
        zmax = log10(zmax)
    if (imscale == 'squareroot'):
        image = sqrt(image)
        zmin = sqrt(zmin)
        zmax = sqrt(zmax)

    return image, zmin, zmax


# ------------------------------------------
# plot mask borders in CCD coordinates

def borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,bit,lcolor,lstyle,lwidth):

    for i in range(1,ydim):
        for j in range(1,xdim):
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

# corner cases

    for j in range(ydim):
        try:
            if kepstat.bitInBitmap(maskimg[j,0],bit) and not kepstat.bitInBitmap(maskimg[j-1,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[j+1,0],bit) and kepstat.bitInBitmap(maskimg[j,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) + 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j,0],bit):
            x = array([pixcoord1[0,j],pixcoord1[0,j]]) - 0.5
            try:
                y = array([pixcoord2[0,j],pixcoord2[0,j+1]]) - 0.5
            except:
                y = array([pixcoord2[0,j-1],pixcoord2[0,j]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j,xdim-1],bit):
            x = array([pixcoord1[xdim-1,j],pixcoord1[xdim-1,j]]) + 0.5
            try:
                y = array([pixcoord2[xdim-1,j],pixcoord2[xdim-1,j+1]]) - 0.5
            except:
                y = array([pixcoord2[xdim-1,j-1],pixcoord2[xdim-1,j]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
    for i in range(xdim):
        try:
            if kepstat.bitInBitmap(maskimg[0,i],bit) and not kepstat.bitInBitmap(maskimg[0,i-1],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) - 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[0,i+1],bit) and kepstat.bitInBitmap(maskimg[0,i],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) + 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0,i],bit):
            try:
                x = array([pixcoord1[i,0],pixcoord1[i+1,0]]) - 0.5
            except:
                x = array([pixcoord1[i-1,0],pixcoord1[i,0]]) + 0.5
            y = array([pixcoord2[i,0],pixcoord2[i,0]]) - 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim-1,i],bit):
            try:
                x = array([pixcoord1[i,ydim-1],pixcoord1[i+1,ydim-1]]) - 0.5
            except:
                x = array([pixcoord1[i-1,ydim-1],pixcoord1[i,ydim-1]]) - 0.5
            y = array([pixcoord2[i,ydim-1],pixcoord2[i,ydim-1]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[ydim-1,xdim-1],bit):
        x = array([pixcoord1[xdim-2,ydim-1],pixcoord1[xdim-1,ydim-1]]) + 0.5
        y = array([pixcoord2[xdim-1,ydim-1],pixcoord2[xdim-1,ydim-1]]) + 0.5
        plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0,xdim-1],bit):
        x = array([pixcoord1[xdim-1,0],pixcoord1[xdim-1,0]]) + 0.5
        y = array([pixcoord2[xdim-1,0],pixcoord2[xdim-1,1]]) - 0.5
        plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    return


# ------------------------------------------
# plot mask borders in CCD coordinates

def PrfBorders(maskimg,xdim,ydim,pixcoord1,pixcoord2,bit,lcolor,lstyle,lwidth):

    for i in range(1,ydim):
        for j in range(1,xdim):
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                plt.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                plt.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                plt.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                plt.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)

# corner cases

    for j in range(ydim):
        try:
            if kepstat.bitInBitmap(maskimg[j,0],bit) and not kepstat.bitInBitmap(maskimg[j-1,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[j+1,0],bit) and kepstat.bitInBitmap(maskimg[j,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) + 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j,0],bit):
            x = array([pixcoord1[0,j],pixcoord1[0,j]]) - 0.5
            try:
                y = array([pixcoord2[0,j],pixcoord2[0,j+1]]) - 0.5
            except:
                y = array([pixcoord2[0,j-1],pixcoord2[0,j]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j,xdim-1],bit):
            x = array([pixcoord1[xdim-1,j],pixcoord1[xdim-1,j]]) + 0.5
            try:
                y = array([pixcoord2[xdim-1,j],pixcoord2[xdim-1,j+1]]) - 0.5
            except:
                y = array([pixcoord2[xdim-1,j-1],pixcoord2[xdim-1,j]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
    for i in range(xdim):
        try:
            if kepstat.bitInBitmap(maskimg[0,i],bit) and not kepstat.bitInBitmap(maskimg[0,i-1],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) - 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[0,i+1],bit) and kepstat.bitInBitmap(maskimg[0,i],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) + 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0,i],bit):
            try:
                x = array([pixcoord1[i,0],pixcoord1[i+1,0]]) - 0.5
            except:
                x = array([pixcoord1[i-1,0],pixcoord1[i,0]]) + 0.5
            y = array([pixcoord2[i,0],pixcoord2[i,0]]) - 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim-1,i],bit):
            try:
                x = array([pixcoord1[i,ydim-1],pixcoord1[i+1,ydim-1]]) - 0.5
            except:
                x = array([pixcoord1[i-1,ydim-1],pixcoord1[i,ydim-1]]) - 0.5
            y = array([pixcoord2[i,ydim-1],pixcoord2[i,ydim-1]]) + 0.5
            plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[ydim-1,xdim-1],bit):
        x = array([pixcoord1[xdim-2,ydim-1],pixcoord1[xdim-1,ydim-1]]) + 0.5
        y = array([pixcoord2[xdim-1,ydim-1],pixcoord2[xdim-1,ydim-1]]) + 0.5
        plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0,xdim-1],bit):
        x = array([pixcoord1[xdim-1,0],pixcoord1[xdim-1,0]]) + 0.5
        y = array([pixcoord2[xdim-1,0],pixcoord2[xdim-1,1]]) - 0.5
        plt.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    return
