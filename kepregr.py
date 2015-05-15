"""
Regression based on k-nearest neighbor algorithm
"""

from pyraf import iraf
import numpy, sys, time, pyfits, pylab
from pyfits import *
from pylab import *
from matplotlib import *
from math import *
import kepio, kepmsg, kepkey, kepfunc
from scipy.spatial import KDTree

def knn_predict(X, y, method, k_neighb):

    X_ = numpy.array(X, dtype='float64')
    X = X_.reshape(len(X), 1)
    y = numpy.array(y, dtype='float64')

    if method == "kdtree":
        tree = KDTree(X)
        # Index =  1 , for the second array in the
        # KDTree query return. This corresponds to the indices
        # of the k nearest neighbors.
        nb_ind = tree.query(X, k=k_neighb)[1]

    else:
        # A naive euclidean distance calculation.
        # as we have n_features = 1, 
        eucl_d = []
        for i in X:
            eucl_d.append(abs(i - X_))
        eucl_d = np.array(eucl_d)

        nb_ind = eucl_d.argsort(axis=1)[:, :k_neighb]

    return numpy.mean(y[nb_ind], axis=1)


def kepregr(infile,outfile,datacol,kmethod,kneighb,plot,plotlab,
              clobber,verbose,logfile,status): 
    """
    Perform a k-nearest neighbor regression analysis.
    """

## startup parameters

    status = 0
    labelsize = 24
    ticksize = 16
    xsize = 16
    ysize = 6
    lcolor = '#47AE10'
    lwidth = 1.0
    fcolor = '#9AFF9A'
    falpha = 0.3

## log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPREGR -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'kmethod='+str(kmethod)+' '
    call += 'kneighb='+str(kneighb)+' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    call += 'plotlab='+str(plotlab)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

## start time

    kepmsg.clock('KEPREGR started at',logfile,verbose)

## test log file

    logfile = kepmsg.test(logfile)

## clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPREGR: ' + outfile + ' exists. Use clobber=yes'
	    status = kepmsg.err(logfile,message,verbose)

## open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)
        if cadence == 0.0: 
            tstart, tstop, ncad, cadence, status = kepio.cadence(instr,infile,logfile,verbose,status) 
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence

## fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

## read table structure

    if status == 0:
	table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# read time and flux columns

    if status == 0:
        barytime, status = kepio.readtimecol(infile,table,logfile,verbose)
    if status == 0:
        flux, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)

# filter input data table

    if status == 0:
        try:
            nanclean = instr[1].header['NANCLEAN']
        except:
            naxis2 = 0
            for i in range(len(table.field(0))):
                if (numpy.isfinite(barytime[i]) and numpy.isfinite(flux[i]) and flux[i] != 0.0):
                    table[naxis2] = table[i]
                    naxis2 += 1
            instr[1].data = table[:naxis2]
            comment = 'NaN cadences removed from data'
            status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)

## read table columns

    if status == 0:
	try:
            intime = instr[1].data.field('barytime')
	except:
            intime, status = kepio.readfitscol(infile,instr[1].data,'time',logfile,verbose)
	indata, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)
    if status == 0:
        intime = intime + bjdref
        indata = indata / cadenom

    if status == 0:
   	outdata = knn_predict(intime, indata, kmethod, kneighb)

## comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

## clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
        if intime0 < 2.4e6: intime0 += 2.4e6
	ptime = intime - intime0
        # print ptime,intime,intime0
	xlab = 'BJD $-$ %d' % intime0

## clean up y-axis units

    if status == 0:
        pout = indata * 1.0
        pout2 = outdata * 1.0 
	nrm = len(str(int(numpy.nanmax(pout))))-1
	pout = pout / 10**nrm
	pout2 = pout2 / 10**nrm
	ylab = '10$^%d$ %s' % (nrm, plotlab)

## data limits

	xmin = numpy.nanmin(ptime)
	xmax = numpy.nanmax(ptime)
	ymin = numpy.min(pout)
	ymax = numpy.nanmax(pout)
	xr = xmax - xmin
	yr = ymax - ymin
        ptime = insert(ptime,[0],[ptime[0]]) 
        ptime = append(ptime,[ptime[-1]])
        pout = insert(pout,[0],[0.0]) 
        pout = append(pout,0.0)
        pout2 = insert(pout2,[0],[0.0]) 
        pout2 = append(pout2,0.0)

## plot light curve

    if status == 0 and plot:
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': labelsize,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': ticksize,
                      'ytick.labelsize': ticksize}
            rcParams.update(params)
        except:
            print 'ERROR -- KEPREGR: install latex for scientific plotting'
            status = 1
    if status == 0 and plot:
        pylab.figure(1,figsize=[xsize,ysize])

## plot regression data

        ax = pylab.axes([0.06,0.1,0.93,0.87])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        # pylab.plot(ptime,pout,color='#ff9900',linestyle='-',linewidth=lwidth)
	pylab.scatter(ptime, pout, color='#214CAE', s=5)
        fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
        pylab.plot(ptime[kneighb:-kneighb],pout2[kneighb:-kneighb],color=lcolor,linestyle='-',linewidth=lwidth*2.0)
        xlabel(xlab, {'color' : 'k'})
        ylabel(ylab, {'color' : 'k'})
        xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0: 
            ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            ylim(1.0e-10,ymax+yr*0.01)
        pylab.grid()
        pylab.draw()
        pylab.savefig(re.sub('\.\S+','.png',outfile),dpi=100)

## write output file

    if status == 0:
        for i in range(len(outdata)):
            instr[1].data.field(datacol)[i] = outdata[i]
        instr.writeto(outfile)
    
## close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

## end time

    if (status == 0):
	    message = 'KEPREGR completed at'
    else:
	    message = '\nKEPREGR aborted at'
    kepmsg.clock(message,logfile,verbose)

## main

parfile = iraf.osfn("kepler$kepregr.par")
t = iraf.IrafTaskFactory(taskname="kepregr", value=parfile, function=kepregr)
