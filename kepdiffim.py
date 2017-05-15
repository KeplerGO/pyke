import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey, kepplot, kepstat
import sys, time, re, math

# -----------------------------------------------------------
# core code

def kepdiffim(infile,outfile,plotfile,imscale,colmap,filter,function,cutoff,
              clobber,verbose,logfile,status,cmdLine=False):
# input arguments

    status = 0
    np.seterr(all="ignore")

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDIFFIM -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'plotfile='+plotfile+' '
    call += 'imscale='+imscale+' '
    call += 'colmap='+colmap+' '
    filt = 'n'
    if (filter): filt = 'y'
    call += 'filter='+filt+ ' '
    call += 'function='+function+' '
    call += 'cutoff='+str(cutoff)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPDIFFIM started at: ',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPDIFFIM: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# reference color map

    if colmap == 'browse':
        status = cmap_plot()

# open TPF FITS file

    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
            kepio.readTPF(infile,'TIME',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, tcorr, status = \
            kepio.readTPF(infile,'TIMECORR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, cadno, status = \
            kepio.readTPF(infile,'CADENCENO',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, fluxpixels, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, errpixels, status = \
            kepio.readTPF(infile,'FLUX_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)

# read mask defintion data from TPF file

    if status == 0:
        maskimg, pixcoord1, pixcoord2, status = kepio.readMaskDefinition(infile,logfile,verbose)

# print target data

    if status == 0:
        print ''
        print '      KepID:  %s' % kepid
        print ' RA (J2000):  %s' % ra
        print 'Dec (J2000): %s' % dec
        print '     KepMag:  %s' % kepmag
        print '   SkyGroup:    %2s' % skygroup
        print '     Season:    %2s' % str(season)
        print '    Channel:    %2s' % channel
        print '     Module:    %2s' % module
        print '     Output:     %1s' % output
        print ''

# how many quality = 0 rows?

    if status == 0:
        npts = 0
        nrows = len(fluxpixels)
        for i in range(nrows):
            if qual[i] == 0 and np.isfinite(barytime[i]) \
                    and np.isfinite(fluxpixels[i,int(ydim*xdim/2)]):
                npts += 1
        time = np.empty((npts))
        timecorr = np.empty((npts))
        cadenceno = np.empty((npts))
        quality = np.empty((npts))
        pixseries = np.empty((ydim*xdim,npts))
        errseries = np.empty((ydim*xdim,npts))

# construct output light curves

    if status == 0:
        nptsx = 0
        for i in range(ydim*xdim):
            npts = 0
            for k in range(nrows):
                if (qual[k] == 0 and
                    np.isfinite(barytime[k]) and
                    np.isfinite(fluxpixels[k,int(ydim*xdim/2)])):
                    time[npts] = barytime[k]
                    timecorr[npts] = tcorr[k]
                    cadenceno[npts] = cadno[k]
                    quality[npts] = qual[k]
                    pixseries[i,npts] = fluxpixels[k,nptsx]
                    errseries[i,npts] = errpixels[k,nptsx]
                    npts += 1
            nptsx += 1

# define data sampling

    if status == 0 and filter:
        tpf, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0 and filter:
        cadence, status = kepkey.cadence(tpf[1],infile,logfile,verbose)
        tr = 1.0 / (cadence / 86400)
        timescale = 1.0 / (cutoff / tr)

# define convolution function

    if status == 0 and filter:
        if function == 'boxcar':
            filtfunc = np.ones(int(np.ceil(timescale)))
        elif function == 'gauss':
            timescale /= 2
            dx = np.ceil(timescale * 10 + 1)
            filtfunc = kepfunc.gauss()
            filtfunc = filtfunc([1.0,dx/2-1.0,timescale],linspace(0,dx-1,dx))
        elif function == 'sinc':
            dx = np.ceil(timescale * 12 + 1)
            fx = linspace(0,dx-1,dx)
            fx = fx - dx / 2 + 0.5
            fx /= timescale
            filtfunc = np.sinc(fx)
        filtfunc /= np.sum(filtfunc)

# pad time series at both ends with noise model

    if status == 0 and filter:
        for i in range(ydim * xdim):
            ave, sigma  = (np.mean(pixseries[i, :len(filtfunc)]),
                           np.std(pixseries[i, :len(filtfunc)]))
            padded = np.append(kepstat.randarray(np.ones(len(filtfunc)) * ave,
                    np.ones(len(filtfunc)) * sigma), pixseries[i,:])
            ave, sigma  = (np.mean(pixseries[i,-len(filtfunc):]),
                           np.std(pixseries[i,-len(filtfunc):]))
            padded = np.append(padded, kepstat.randarray(np.ones(len(filtfunc))
                    * ave, np.ones(len(filtfunc)) * sigma))

# convolve data
            if status == 0:
                convolved = np.convolve(padded,filtfunc,'same')
# remove padding from the output array
            if status == 0:
                outdata = convolved[len(filtfunc):-len(filtfunc)]
# subtract low frequencies
            if status == 0:
                outmedian = np.median(outdata)
                pixseries[i,:] = pixseries[i,:] - outdata + outmedian

# sum pixels over cadence

    if status == 0:
        nptsx = 0
        nrows = len(fluxpixels)
        pixsum = np.zeros((ydim*xdim))
        errsum = np.zeros((ydim*xdim))
        for i in range(npts):
            if quality[i] == 0:
                pixsum += pixseries[:,i]
                errsum += errseries[:,i]**2
                nptsx += 1
        pixsum /= nptsx
        errsum = np.sqrt(errsum) / nptsx

# calculate standard deviation pixels

    if status == 0:
        pixvar = np.zeros((ydim*xdim))
        for i in range(npts):
            if quality[i] == 0:
                pixvar += (pixsum - pixseries[:,i] / errseries[:,i])**2
        pixvar = np.sqrt(pixvar)

# median pixel errors

    if status == 0:
        errmed = np.empty((ydim*xdim))
        for i in range(ydim*xdim):
            errmed[i] = np.median(errseries[:,i])

# calculate chi distribution pixels

    if status == 0:
        pixdev = np.zeros((ydim*xdim))
        for i in range(npts):
            if quality[i] == 0:
                pixdev += ((pixsum - pixseries[:,i]) / pixsum)**2
        pixdev = np.sqrt(pixdev)

# image scale and intensity limits

    if status == 0:
        pixsum_pl, zminsum, zmaxsum = kepplot.intScale1D(pixsum,imscale)
        pixvar_pl, zminvar, zmaxvar = kepplot.intScale1D(pixvar,imscale)
        pixdev_pl, zmindev, zmaxdev = kepplot.intScale1D(pixdev,imscale)

# construct output summed image

    if status == 0:
        imgsum = np.empty((ydim,xdim))
        imgvar = np.empty((ydim,xdim))
        imgdev = np.empty((ydim,xdim))
        imgsum_pl = np.empty((ydim,xdim))
        imgvar_pl = np.empty((ydim,xdim))
        imgdev_pl = np.empty((ydim,xdim))
        n = 0
        for i in range(ydim):
            for j in range(xdim):
                imgsum[i,j] = pixsum[n]
                imgvar[i,j] = pixvar[n]
                imgdev[i,j] = pixdev[n]
                imgsum_pl[i,j] = pixsum_pl[n]
                imgvar_pl[i,j] = pixvar_pl[n]
                imgdev_pl[i,j] = pixdev_pl[n]
                n += 1

# construct output file

    if status == 0:
        instruct, status = kepio.openfits(infile,'readonly',logfile,verbose)
        status = kepkey.history(call,instruct[0],outfile,logfile,verbose)
        hdulist = pyfits.HDUList(instruct[0])
        hdulist.writeto(outfile)
        status = kepkey.new('EXTNAME','FLUX','name of extension',instruct[2],outfile,logfile,verbose)
        pyfits.append(outfile,imgsum,instruct[2].header)
        status = kepkey.new('EXTNAME','CHI','name of extension',instruct[2],outfile,logfile,verbose)
        pyfits.append(outfile,imgvar,instruct[2].header)
        status = kepkey.new('EXTNAME','STDDEV','name of extension',instruct[2],outfile,logfile,verbose)
        pyfits.append(outfile,imgdev,instruct[2].header)
        status = kepkey.new('EXTNAME','APERTURE','name of extension',instruct[2],outfile,logfile,verbose)
        pyfits.append(outfile,instruct[2].data,instruct[2].header)
        status = kepio.closefits(instruct,logfile,verbose)

# pixel limits of the subimage

    if status == 0:
        ymin = row
        ymax = ymin + ydim
        xmin = column
        xmax = xmin + xdim

# plot limits for summed image

        ymin = float(ymin) - 0.5
        ymax = float(ymax) - 0.5
        xmin = float(xmin) - 0.5
        xmax = float(xmax) - 0.5

# plot style

    if status == 0:
        plotimage(imgsum_pl,imgvar_pl,imgdev_pl,zminsum,zminvar,zmindev,
                  zmaxsum,zmaxvar,zmaxdev,xmin,xmax,ymin,ymax,colmap,plotfile,cmdLine)

# stop time

    kepmsg.clock('KEPDIFFIM ended at: ',logfile,verbose)

    return

# -----------------------------------------------------------
# plot channel image

def plotimage(imgsum_pl,imgvar_pl,imgdev_pl,zminsum,zminvar,zmindev,
              zmaxsum,zmaxvar,zmaxdev,xmin,xmax,ymin,ymax,colmap,plotfile,cmdLine):

    plt.figure(figsize=[15,6])
    #ion()
    plt.clf()

# plot the image window

    ax = plt.axes([0.04,0.11,0.31,0.78])
    plt.imshow(imgsum_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zminsum,vmax=zmaxsum,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    plt.gca().set_autoscale_on(False)
    labels = ax.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    plt.ylabel('Pixel Row Number', {'color' : 'k'})
    plt.title('Flux', {'color' : 'k', 'fontsize' : '24'})

# plot the variance window

    plt.axes([0.36,0.11,0.31,0.78])
    plt.imshow(imgvar_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zminvar,vmax=zmaxvar,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    plt.gca().set_autoscale_on(False)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.setp(plt.gca(),yticklabels=[])
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    try:
        plt.title(r'$\chi$ Distribution', {'color' : 'k', 'fontsize' : '28'})
    except:
        plt.title('Chi Distribution', {'color' : 'k', 'fontsize' : '24'})

# plot the normalized standard deviation window

    plt.axes([0.68,0.11,0.31,0.78])
    plt.imshow(imgdev_pl,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zmindev,vmax=zmaxdev,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    plt.gca().set_autoscale_on(False)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.setp(plt.gca(),yticklabels=[])
    plt.xlabel('Pixel Column Number', {'color' : 'k'})
    plt.title('Normalized Standard Deviation', {'color' : 'k', 'fontsize' : '24'})

# render plot
    plt.ion()
    plt.show()

    if plotfile.lower() != 'none':
        plt.savefig(plotfile)

    return

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot():

    plt.figure(1,figsize=[5,10])
    #ion()
    a=outer(ones(10),arange(0,1,0.01))
    plt.subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        print m
        plt.subplot(l,1,i+1)
        plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        plt.imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        plt.text(100.85,0.5,m,fontsize=10)
    #ioff()
    status = 1
    return status

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Difference imaging of pixels within a target mask')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--plotfile', default='None', help='name of output PNG plot file', type=str)
    parser.add_argument('--imscale', default='logarithmic', help='type of image intensity scale', type=str, choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap', type=str)
    parser.add_argument('--filter', action='store_true', help='High-pass Filter data?')

    parser.add_argument('--function', help='filter function', default='boxcar', type=str, choices=['boxcar','gauss','sinc'])
    parser.add_argument('--cutoff', help='Characteristic frequency cutoff of filter [1/days]', type=int, default=1.0)

    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepdiffim.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()

    cmdLine=True
    kepdiffim(args.infile, args.outfile, args.plotfile, args.imscale, args.cmap, args.filter, args.function, args.cutoff,
            args.clobber, args.verbose, args.logfile, args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepdiffim.par")
    t = iraf.IrafTaskFactory(taskname="kepdiffim", value=parfile, function=kepdiffim)
