import sys, time, math, re
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import numpy as np
import kepio, kepmsg, kepkey, kepfit, kepstat, kepfourier, keplab

def kepdynamic(infile,outfile,fcol,pmin,pmax,nfreq,deltat,nslice,
               plot,plotscale,cmap,clobber,verbose,logfile,status,cmdLine=False):

# startup parameters

    status = 0
    labelsize = 24
    ticksize = 16
    xsize = 12
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2
    np.seterr(all="ignore")

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDYNAMIC -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'fcol='+fcol+' '
    call += 'pmin='+str(pmin)+' '
    call += 'pmax='+str(pmax)+' '
    call += 'nfreq='+str(nfreq)+' '
    call += 'deltat='+str(deltat)+' '
    call += 'nslice='+str(nslice)+' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    call += 'plotscale='+plotscale+ ' '
    call += 'cmap='+str(cmap)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('Start time is',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# error checking

    if status == 0 and pmin >= pmax:
        message = 'ERROR -- KEPDYNAMIC: PMIN must be less than PMAX'
        status = kepmsg.err(logfile,message,verbose)


# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPDYNAMIC: ' + outfile + ' exists. Use clobber'
        status = kepmsg.err(logfile,message,verbose)

# plot color map

    if status == 0 and cmap == 'browse':
        status = keplab.cmap_plot()

# open input file

    if status == 0:
	    instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence

# fudge non-compliant FITS keywords with no values

    if status == 0:
	    instr = kepkey.emptykeys(instr,file,logfile,verbose)

# read table columns

    if status == 0:
        barytime, status = kepio.readtimecol(infile,instr[1].data,logfile,verbose)
    if status == 0:
        signal, status = kepio.readfitscol(infile,instr[1].data,fcol,logfile,verbose)
    if status == 0:
        barytime = barytime + bjdref
        signal = signal / cadenom

# remove infinite data from time series

    if status == 0:
	    incols = [barytime, signal]
	    outcols = kepstat.removeinfinlc(signal, incols)
	    barytime = outcols[0]
	    signal = outcols[1]

# period to frequency conversion

    if status == 0:
        fmin = 1.0 / pmax
        fmax = 1.0 / pmin
        deltaf = (fmax - fmin) / nfreq

# determine bounds of time slices

    if status == 0:
        t1 = []; t2 = []
        dt = barytime[-1] - barytime[0]
        dt -= deltat
        if dt < 0:
            message = 'ERROR -- KEPDYNAMIC: time slices are larger than data range'
            status = kepmsg.err(logfile,message,verbose)
        ds = dt / (nslice - 1)
        for i in range(nslice):
            t1.append(barytime[0] + ds * float(i))
            t2.append(barytime[0] + deltat + ds * float(i))

# loop through time slices

    if status == 0:
        dynam = []
        for i in range(nslice):
            x = []; y = []
            for j in range(len(barytime)):
                if (barytime[j] >= t1[i] and barytime[j] <= t2[i]):
                    x.append(barytime[j])
                    y.append(signal[j])
            x = np.array(x,dtype='float64')
            y = np.array(y,dtype='float')
            y = y - np.median(y)

# determine FT power

	    fr, power = kepfourier.ft(x,y,fmin,fmax,deltaf,False)
            for j in range(len(power)):
                dynam.append(power[j])
            print 'Timeslice: %.4f  Pmax: %.2E' % ((t2[i] + t1[i]) / 2, power.max())

# define shape of results array

        dynam = np.array(dynam,dtype='float64')
        dynam.shape = len(t1),len(power)

# write output file

    if status == 0:
        instr.append(pyfits.ImageHDU())
        instr[-1].data = dynam.transpose()
        instr[-1].header['EXTNAME'] = ('DYNAMIC FT','extension name')
        instr[-1].header['WCSAXES'] = (2,'number of WCS axes')
        instr[-1].header['CRPIX1' ] = (0.5,'reference pixel along axis 1')
        instr[-1].header['CRPIX2' ] = (0.5,'reference pixel along axis 2')
        instr[-1].header['CRVAL1' ] = (t1[0],'time at reference pixel (BJD)')
        instr[-1].header['CRVAL2' ] = (fmin,'frequency at reference pixel (1/day)')
        instr[-1].header['CDELT1' ] = ((barytime[-1] - barytime[0]) / nslice,
                        'pixel scale in dimension 1 (days)')
        instr[-1].header['CDELT2'] = (deltaf,'pixel scale in dimension 2 (1/day)')
        instr[-1].header['CTYPE1'] = ('BJD','data type of dimension 1')
        instr[-1].header['CTYPE2'] = ('FREQUENCY','data type of dimension 2')
        instr.writeto(outfile)

# history keyword in output file

    if status == 0:
	    status = kepkey.history(call,instr[0],outfile,logfile,verbose)

# close input file

    if status == 0:
	    status = kepio.closefits(instr,logfile,verbose)

# clean up x-axis unit

    if status == 0:
	time0 = float(int(barytime[0] / 100) * 100.0)
	barytime = barytime - time0
	xlab = 'BJD $-$ %d' % time0

# image intensity min and max

    if status == 0:
        if 'rithmic' in plotscale:
            dynam = np.log10(dynam)
        elif 'sq' in plotscale:
            dynam = np.sqrt(dynam)
        elif 'logoflog' in plotscale:
            dynam = np.log10(np.abs(np.log10(dynam)))
#        dynam = -dynam
        nstat = 2; pixels = []
        for i in range(dynam.shape[0]):
            for j in range(dynam.shape[1]):
                pixels.append(dynam[i,j])
        pixels = np.array(np.sort(pixels),dtype='float')
        if int(float(len(pixels)) * 0.1 + 0.5) > nstat:
            nstat = int(float(len(pixels)) * 0.1 + 0.5)
        zmin = np.median(pixels[:nstat])
        zmax = np.median(pixels[-1:])
        if np.isnan(zmax):
            zmax = np.median(pixels[-nstat/2:])
        if np.isnan(zmax):
            zmax = np.nanmax(pixels)

# plot power spectrum

    if status == 0 and plot:
        plt.figure(1,figsize=[xsize,ysize])
        plt.clf()
        plt.axes([0.08,0.113,0.91,0.86])
        dynam = dynam.transpose()
        plt.imshow(dynam,origin='lower',aspect='auto',cmap=cmap,vmin=zmin,vmax=zmax,
                     extent=[barytime[0],barytime[-1],fmin,fmax],interpolation='bilinear')
        plt.xlabel(xlab, {'color' : 'k'})
        plt.ylabel(r'Frequency (d$^{-1}$)', {'color' : 'k'})
        plt.grid()
        plt.savefig(re.sub('\.\S+','.png',outfile),dpi=100)

# render plot
        plt.ion()
        plt.show()

    return status

## end time

    if (status == 0):
	    message = 'KEPDYNAMIC completed at'
    else:
	    message = '\nKEPDYNAMIC aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Construct a dynamic (time-dependent) power spectrum from Kepler time series data')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column to plot', type=str, dest='fcol')

    parser.add_argument('--pmin', default=0.1, help='Minimum search period [days]', type=float)
    parser.add_argument('--pmax', default=10., help='Maximum search period [days]', type=float)
    parser.add_argument('--nfreq', default=100, help='Number of frequency intervals', type=int)
    parser.add_argument('--deltat', default=10., help='Length of time slice [days]',type=float)
    parser.add_argument('--nslice', default=10., help='Number of time slices', type=int)

    parser.add_argument('--plot', action='store_true', help='Plot result?')
    parser.add_argument('--plotscale', default='logarithmic', help='type of image intensity scale',
                        type=str, choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap', type=str)


    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log',
                        dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()

    cmdLine=True

    kepynamic(args.infile, args.outfile, args.fcol, args.pmin, args.pmax, args.nfreq, args.deltat, args.nslice,
          args.plot,args.plotscale,args.cmap,args.clobber,args.verbose,args.logfile,args.status,cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepdynamic.par")
    t = iraf.IrafTaskFactory(taskname="kepdynamic", value=parfile, function=kepdynamic)
