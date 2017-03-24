#!/usr/bin/env python

from pyraf import iraf
import numpy, sys
from scipy.interpolate import interp1d
from scipy import integrate
from math import *
import kepmsg

# -----------------------------------------------------------
# rebin cadence data

# method -- the method of interpolation, options are:
#       'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic' or
#	 an integer (i), which interpolates using a spline of order (i)

# interpm -- method of performing the integration
#	 can either be 'quad' or 'romberg'

def rebin(date,flux,logfile,verbose,nbins=None,binwidth=None,ownbins=None,method='linear',interpm='quad'):

    status = 0
    bflux = []

# catch wrong sized input arrays

    if numpy.shape(date) != numpy.shape(flux):
        txt = 'ERROR -- KEPREBIN.REBIN: Time and flux arrays are different lengths'
        status = kepmsg.err(logfile,txt,verbose)

# catch multipe rebinning methods

    if status == 0:
	i = 0
	for binMethod in [nbins,binwidth,ownbins]:
            if binMethod is None:
                i += 1
        if i != 2:
            txt = 'ERROR -- KEPREBIN.REBIN: Supply one and only one of nbin, binwidth, ownbins'
            status = kepmsg.err(logfile,txt,verbose)

# catch bad choice of interpolation method

    if status == 0 and not method in [i for i in ['linear','nearest', 'zero', 'slinear',
                                                  'quadratic' or 'cubic']]:
        try:
            testcase =  "Integer %d" % int(method)
        except:
            txt = 'ERROR -- KEPREBIN.REBIN: Method needs to be one of: linear, nearest, zero, slinear, '
            txt += 'quadratic, cubic or an integer (i), which interpolates using a  spline of order (i)'
        status = kepmsg.err(logfile,txt,verbose)

# catch bad choice of integration method

    if status == 0 and not interpm in [i for i in ['quad','romberg']]:
        txt = 'ERROR -- KEPREBIN.REBIN: Integration method must either be quad or romberg'
        status = kepmsg.err(logfile,txt,verbose)

# perform the interpolation on the data

    if status == 0:
        try:
            intpl = scipy.interpolate.interp1d(date,flux,kind=method)
        except:
            txt = 'ERROR -- KEPREBIN.REBIN: Cannot define interpolation function with scipy.interpolate.interp1d'
            status = kepmsg.err(logfile,txt,verbose)

# calculate bin bounds, should be of length nbins+1

    if status == 0:
	if binwidth is not None:
            bounds = numpy.arange(date[0],date[-1] + 1.0e-10,binwidth)
            bdate = numpy.arange(date[0] + 0.5 * binwidth,date[-1] - 0.5 * binwidth,binwidth)
	elif nbins is not None:
            bounds = numpy.arange(date[0],date[-1] + 1.0e-10,(date[-1] - date[0]) / nbins)
            bdate = numpy.arange(date[0] + 0.5 * ((date[-1] - date[0]) / nbins), \
                                     date[-1],(date[-1] - date[0]) / nbins)
	elif ownbins is not None:
            bounds = ownbins
            bdate = numpy.zeros(len(bounds) - 1)
            for i in range(len(bdate)):
                bdate[i] = bounds[i] + (0.5 * (bounds[i+1] - bounds[i]))
        minbin = []
	for i in range(1,len(bounds)):
            minbin.append(bounds[i] - bounds[i-1])
        mincad = []
	for i in range(1,len(date)):
            mincad.append(date[i] - date[i-1])
	
# iterate over the bins starting with bounds[1]

    if status == 0:
	for i in range(1,len(bounds)):
            bin = []
            t = []
            extrabitHigh = 0.0
            extrabitLow = 0.0
            eb1 = None
            eb2 = None
		
# iterate over the length of the flux array

            for j in range(len(flux)):

# catch all date bins within the bounds of the bin

                if date[j] >= bounds[i-1] and date[j] < bounds[i]:
                    bin.append(flux[j])
                    t.append(j)
		try:
                    tmin = min(t)
                    tmax = max(t)
		
# deal with the bits either side of date[tmin] and date[tmax], first the bit between date[tmax] and bounds[i]
		
                    if bounds[i] > date[tmax] + 1.0e-6:
                        if interpm == 'quad':
                            extrabitHigh = scipy.integrate.quad(intpl,date[tmax],bounds[i])[0]
                        elif interpm == 'romberg':
                            extrabitHigh = scipy.integrate.romberg(intpl,date[tmax],bounds[i])

# now the bit between bound[i-1] and t[min]

                    if bounds[i-1] < date[tmin] - 1.0e-6:
                        if interpm == 'quad':
                            extrabitLow = scipy.integrate.quad(intpl,bounds[i-1],date[tmin])[0]
                        elif interpm == 'romberg':
                            extrabitLow = scipy.integrate.romberg(intpl,bounds[i-1],date[tmin])
					
# scale the extrabits as if they were the size of a full bin, catch error if there is no extrabit

                    try:
                        eb1 = (extrabitHigh / (float(bounds[i] - date[tmax])))
                        bin.append(eb1)
                    except ZeroDivisionError:
                        pass
                    try:
                        eb2 = (extrabitLow / (float(date[tmin] - bounds[i-1])))
                        bin.append(eb2)
                    except ZeroDivisionError:
                        pass
		
                except ValueError:
                    if interpm == 'quad':
                        totbin = scipy.integrate.quad(intpl,bounds[i-1],bounds[i])[0]
                    elif interpm == 'romberg':
                        totbin =  integrate.romberg(intpl,bounds[i-1],bounds[i])
                    eb = (totbin / (float(bounds[i] - bounds[i-1])))
                    bin.append(eb)
		
                binflux = numpy.mean(bin)
		bflux.append(binflux)

    return bdate, bflux
