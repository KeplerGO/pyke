"""
Simple dip detection using on of 
	- k-nearest neighbour mean algorithm or
	- outlier detection using chebyshev inequality
"""

from pyraf import iraf
import numpy, sys, time, pyfits, pylab
from pyfits import *
from pylab import *
from matplotlib import *
from math import *
import kepio, kepmsg, kepkey, kepfunc

def dip_kneighb_mean(k, i, y):
    y_i = y[i]
    max_left = max(y[i-k:i] - y_i)
    max_right = max(y[i+1:i+k+1] - y_i)
    min_left = min(y[i-k:i] - y_i)
    min_right = min(y[i+1:i+k+1] - y_i)

    # Compute the 'distance' maximum.
    # Here we compute the maximum absolute distance, then 
    #      return with their sign. 
    if (abs(max_left) < abs(min_left)):
        max_left = min_left
    if (abs(max_right) < abs(min_right)):
        max_right = min_right

    return (max_left + max_right)/2.0


def dip_chebyshev(k, i, y):
    y_i = y[i]
    # List of negibhors without y[i]
    N_wo_i = numpy.concatenate((y[i-k:i], y[i+1:i+k+1]))

    N_mean = numpy.mean(N_wo_i)
    N_std = numpy.std(N_wo_i)

    # One-tailed variant of Chebyshev inequality test, 
    # with a preset std multiplier. 
    # This can be made configurable if required.
    if (y_i <= N_mean) and (abs(y_i - N_mean) >= 0.025 * N_std):
        # We caught an outlier.
        return abs(y_i - N_mean)
    else:
        return 0.0

def _find_dips(X, y, method, k, h):

    dip = []; y_out = []; x_out = []
    X = numpy.array(X)
    y = numpy.array(y)

    X_len = len(X)

    if method == "kneigh":
        _dip_func = dip_kneighb_mean
    elif method == "chebyshev":
        _dip_func = dip_chebyshev

    for i in range(X_len):
        # Dip detection is based on sliding windows. For consistency, 
        # we skip the first k and the last k.
        if k > i or (k + i) >= len(y):
            dip.append('0.0')
            continue
        dip.append(_dip_func(k, i, y))

    dip = numpy.array(dip, dtype='float64')

    # We need to calculate mean/SD only for the
    #   +ve elements in the list.
    dip_p = [i for i in dip if i > 0]
    dp_mean = numpy.mean(dip_p)
    dp_std = numpy.std(dip_p)
    # print "Mean", dp_mean, "Std", dp_std, "h * dp_std", h * dp_std

    # Remove local dips, that are 'small' in global context
    for i in range(len(dip)):
        if(dip[i] > 0.0 and ((dip[i] - dp_mean) > h * dp_std)):
            y_out.append((y[i], i))
            # print i, ",", y[i], ",", X[i]

    # Retain only one dip within the reach of k.
    # We introduce a 'u' to take care of the proper indexing
    # after we pop elements from the list.
    u = 0
    for t in range(len(y_out) - 1):
        dip_i, i = y_out[t-u]
        dip_j, j = y_out[t-u + 1]
        if abs(i - j) <= k:
            if dip_i == max(dip_i, dip_j):
                y_out.pop(t-u)
            else:
                y_out.pop(t-u + 1)

            # Make sure, next read will include one of the current 
            # elements
            u += 1

    y_fmt = [0.0] * X_len
    for i in range(len(y_out)):
        x_out.append(X[y_out[i][1]])
        y_fmt[y_out[i][1]] = y_out[i][0]

    y_out = (numpy.array(y_out))[:,0]
    return numpy.array(x_out), y_out, numpy.array(y_fmt)


def kepdip(infile,outfile,datacol,dmethod,kneighb,hstd,plot,plotlab,
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
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#9AFF9A'
    falpha = 0.3

## log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPDIP -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'dmethod='+dmethod+' '
    call += 'hstd='+str(hstd)+' '
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

    kepmsg.clock('KEPDIP started at',logfile,verbose)

## test log file

    logfile = kepmsg.test(logfile)

## clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPDIP: ' + outfile + ' exists. Use clobber=yes'
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

## smooth data

    if status == 0:
        # outdata = knn_predict(intime, indata, kmethod, kneighb)
	outdata_t, outdata_l, outdata_fmt = _find_dips(intime, indata, dmethod, kneighb, hstd)

## comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

## clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
        if intime0 < 2.4e6: intime0 += 2.4e6
	ptime = intime - intime0
	ptime2 = outdata_t - intime0
        # print ptime,intime,intime0
	xlab = 'BJD $-$ %d' % intime0

## clean up y-axis units

    if status == 0:
        pout = indata * 1.0
        pout2 = outdata_l * 1.0 
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
	if (len(ptime2) > 0):
	        ptime2 = insert(ptime2,[0],[ptime2[0]]) 
        	ptime2 = append(ptime2,[ptime2[-1]])
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
            print 'ERROR -- KEPDIP: install latex for scientific plotting'
            status = 1
    if status == 0 and plot:
        pylab.figure(1,figsize=[xsize,ysize])

## plot regression data

        ax = pylab.axes([0.06,0.1,0.93,0.87])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
	pylab.scatter(ptime, pout, color='#214CAE', s=2)

	if (len(ptime2) > 0):
	        pylab.scatter(ptime2, pout2, color='#47AE10', s=35, marker='o', linewidths=2, alpha=0.4)
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
        for i in range(len(outdata_fmt)):
            instr[1].data.field(datacol)[i] = outdata_fmt[i]
        instr.writeto(outfile)
    
## close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

## end time

    if (status == 0):
	    message = 'KEPDIP completed at'
    else:
	    message = '\nKEPDIP aborted at'
    kepmsg.clock(message,logfile,verbose)

## main

parfile = iraf.osfn("kepler$kepdip.par")
t = iraf.IrafTaskFactory(taskname="kepdip", value=parfile, function=kepdip)
