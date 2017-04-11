import numpy as np
import sys
from matplotlib import pyplot as plt
import kepio, kepmsg, kepkey
import re

__svnid__ = "$Id: kepdraw.py 6165 2014-03-26 21:16:27Z mstill $"
__url__ = "$URL: svn+ssh://mstill@murzim.amn.nasa.gov/data-repo/trunk/data/flight/go/PyKE/kepler/kepdraw.py $"


def kepdraw(infile,outfile,datacol,ploterr,errcol,quality,
	    lcolor,lwidth,fcolor,falpha,labelsize,ticksize,
	    xsize,ysize,fullrange,chooserange,y1,y2,plotgrid,
            ylabel,plottype,verbose,logfile,status,cmdLine=False):

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDRAW -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+datacol+' '
    perr = 'n'
    if (ploterr): perr = 'y'
    call += 'ploterr='+perr+ ' '
    call += 'errcol='+errcol+' '
    qual = 'n'
    if (quality): qual = 'y'
    call += 'quality='+qual+ ' '
    call += 'lcolor='+str(lcolor)+' '
    call += 'lwidth='+str(lwidth)+' '
    call += 'fcolor='+str(fcolor)+' '
    call += 'falpha='+str(falpha)+' '
    call += 'labelsize='+str(labelsize)+' '
    call += 'ticksize='+str(ticksize)+' '
    call += 'xsize='+str(xsize)+' '
    call += 'ysize='+str(ysize)+' '
    frange = 'n'
    if (fullrange): frange = 'y'
    call += 'fullrange='+frange+ ' '
    crange = 'n'
    if (chooserange): crange = 'y'
    call += 'chooserange='+crange+ ' '
    call += 'ymin='+str(y1)+' '
    call += 'ymax='+str(y2)+' '
    pgrid = 'n'
    if (plotgrid): pgrid = 'y'
    call += 'plotgrid='+pgrid+ ' '
    call += 'ylabel='+str(ylabel)+' '
    call += 'plottype='+plottype+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPDRAW started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# open input file

    if status == 0:
        struct, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(struct,infile,logfile,verbose,status)

# read table structure

    if status == 0:
	table, status = kepio.readfitstab(infile,struct[1],logfile,verbose)

# read table columns

    if status == 0:
        intime, status = kepio.readtimecol(infile,table,logfile,verbose)
        intime += bjdref
	indata, status = kepio.readfitscol(infile,table,datacol,logfile,verbose)
        indataerr, status = kepio.readfitscol(infile,table,errcol,logfile,verbose)
        qualty, status = kepio.readfitscol(infile,table,'SAP_QUALITY',logfile,verbose)

# close infile

    if status == 0:
	status = kepio.closefits(struct,logfile,verbose)

# remove infinities and bad data

    if status == 0:
        if np.isnan(np.nansum(indataerr)):
            indataerr[:] = 1.0e-5
        work1 = np.array([intime, indata, indataerr, qualty],dtype='float64')
        work1 = np.rot90(work1,3)
        work1 = work1[~np.isnan(work1).any(1)]
        work1 = work1[~np.isinf(work1).any(1)]
        if quality:
            work1 = work1[work1[:,0] == 0.0]
        barytime = np.array(work1[:,3],dtype='float64')
        data = np.array(work1[:,2],dtype='float32')
        dataerr = np.array(work1[:,1],dtype='float32')
        if len(barytime) == 0:
            message = 'ERROR -- KEPDRAW: Plotting arrays are full of NaN'
            status = kepmsg.err(logfile,message,verbose)

# clean up x-axis unit

    if status == 0:
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

#        ax = plt.axes([0.1,0.11,0.89,0.87])
	ax = plt.subplot(111)
	plt.subplots_adjust(0.06,0.15,0.92,0.83)

# force tick labels to be absolute rather than relative

        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        plt.setp(labels, 'rotation', 90, fontsize=ticksize)

# if plot type is 'fast' plot data time series as points

        if plottype == 'fast':
            plt.plot(barytime,data,'o',color=lcolor)

# if plot type is 'pretty' plot data time series as an unbroken line, retaining data gaps

        else:
            ltime = np.array([],dtype='float64')
            ldata = np.array([],dtype='float32')
            dt = 0
            work1 = 2.0 * cadence / 86400
            for i in range(1,len(data)-1):
                dt = barytime[i] - barytime[i-1]
                if dt < work1:
                    ltime = np.append(ltime,barytime[i])
                    ldata = np.append(ldata,data[i])
                else:
                    plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)
                    ltime = np.array([],dtype='float64')
                    ldata = np.array([],dtype='float32')
            plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)

# plot the fill color below data time series, with no data gaps

	plt.fill(barytime,data,fc=fcolor,linewidth=0.0,alpha=falpha)

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

        ax.minorticks_on()
        ax.tick_params('both', length=20, width=2, which='major')
        ax.tick_params('both', length=10, width=1, which='minor')

# save plot to file

    if status == 0 and outfile.lower() != 'none':
	plt.savefig(outfile)

# render plot
        plt.ion()
        plt.show()
# end time

    if (status == 0):
        message = 'KEPDRAW completed at'
    else:
        message = '\nKEPDRAW aborted at'
    kepmsg.clock(message,logfile,verbose)

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse

    parser = argparse.ArgumentParser(description='Interactive plotting of Kepler time series data')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='name of output PNG file', type=str)
    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column to plot', type=str)
    parser.add_argument('--ploterr', action='store_true', help='Plot data error bars?')
    parser.add_argument('--errcol', default='SAP_FLUX_ERR', help='Name of data error column', type=str)
    parser.add_argument('--quality', action='store_true', help='Ignore cadences where the data quality is questionable?')
    parser.add_argument('--lcolor', default='#0000ff', help='HTML color of data line within plot', type=str)
    parser.add_argument('--lwidth', default=1.0, help='type of image intensity scale', type=float)
    parser.add_argument('--fcolor', default='#ffff00', help='HTML color of data line within plot', type=str)
    parser.add_argument('--falpha', default=0.2, help='type of image intensity scale', type=float)
    parser.add_argument('--labelsize', default=24., help='Fontsize of axis labels', type=float)
    parser.add_argument('--ticksize', default=16., help='Fontsize of numeric tick labels', type=float)
    parser.add_argument('--xsize', default=18., help='X-dimension size of plot', type=float)
    parser.add_argument('--ysize', default=6., help='Y-dimension size of plot', type=float)
    parser.add_argument('--fullrange', action='store_true', help='Plot flux range from 0.0 e-/sec?')
    parser.add_argument('--chooserange', action='store_true', help='Choose Y-axis range?')
    parser.add_argument('--ymin', default=0., help='Low limit of the Y-axis range to plot [e-/s]', type=float)
    parser.add_argument('--ymax', default=100000., help='High limit of the Y-axis range to plot [e-/s]', type=float)
    parser.add_argument('--plotgrid', action='store_true', help='Plot axis grid?')
    parser.add_argument('--ylabel', default='e$^-$ s$^{-1}$', help='Plot axis label', type=str)
    parser.add_argument('--plottype', default='fast', help='plot type', type=str, choices=['fast','pretty'])
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    cmdLine=True
    kepdraw(args.infile,args.outfile,args.datacol,args.ploterr,args.errcol,args.quality,
            args.lcolor,args.lwidth,args.fcolor,args.falpha,args.labelsize,args.ticksize,
            args.xsize,args.ysize,args.fullrange,args.chooserange,args.ymin,args.ymax,
            args.plotgrid,args.ylabel,args.plottype,args.verbose,args.logfile,args.status,
            cmdLine)

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepdraw.par")
    t = iraf.IrafTaskFactory(taskname="kepdraw", value=parfile, function=kepdraw)
