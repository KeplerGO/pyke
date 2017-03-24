from pyraf import iraf
import numpy as n
import sys
sys.path.append('/usr/stsci/kepler/' )
import kepio, kepmsg, kepkey

from astropy.io import fits as pyfits
import pylab as plt
from scipy.interpolate import interp1d
from scipy import integrate
from pylab import rc, rcParams

def bin_funct(date,flux,nbins=None,binwidth=None,ownbins=None,method='linear',
	interpm='quad'):
	"""
	A function to bin light curve data

        Arguements:
        date -- the array containing the time the observations were made
	flux -- the flux values corresponding to the date array

	flux and date must be the same length and should have been cleaned of
	nan values

        nbin -- if specifying the number of bins then give the number here
        binwidth -- the width of each bin
        ownbins -- an array which has been read from a user supplied file

	only ONE of nbin, binwidth and ownbins should be defined
	No bin can be larger than the largest gap in data <- should fix

	method -- the method of interpolation, options are:
	'linear', 'nearest', 'zero', 'slinear', 'quadratic' or 'cubic'
	or
	an integer (i), which interpolates using a spline of order (i)

	interpm -- method of performing the integration
	can either be 'quad' or 'romberg'
	"""

	# Try to catch any bad values.
	if n.shape(date) != n.shape(flux):
		raise IndexError('date and flux arrays should be the same length')

	i = 0
	for binMethod in [nbins,binwidth,ownbins]:
		if binMethod is None:
			i += 1
	if i != 2:
		raise TypeError('Supply one and only one of nbin, binwidth, ownbins')

	if not method in [i for i in ['linear','nearest', 'zero', 'slinear',
	                              'quadratic' or 'cubic']]:
		try:
			testcase =  "Integer %d" % int(method)
		except ValueError:
			raise ValueError("method needs to be one of: 'linear', 'nearest', 'zero', 'slinear', 'quadratic' or 'cubic'or an integer (i), which interpolates using a  spline of order (i)")

	if not interpm in [i for i in ['quad','romberg']]:
		raise ValueError("Integration method must either be 'quad' or 'romberg'")



	#perform the interpolation on the data
	intpl = interp1d(date,flux,kind=method)

	#caculate bin bounds using the correct method, should be of length nbins+1
	if binwidth is not None:
		bounds = n.arange(date[0],date[-1]+1e-10,binwidth)
		bdate = n.arange(date[0]+0.5*binwidth,date[-1]-0.5*binwidth,binwidth)
	elif nbins is not None:
		bounds = n.arange(date[0],date[-1]+1e-10,(date[-1]-date[0])/nbins)
		bdate = n.arange(date[0]+0.5*((date[-1]-date[0])/nbins),
			date[-1],(date[-1]-date[0])/nbins)
	elif ownbins is not None:
		bounds = ownbins
		bdate = n.zeros(len(bounds)-1)
		for i in range(len(bdate)):
			bdate[i] = bounds[i]+(0.5*(bounds[i+1] - bounds[i]))


	minbin = []
	for i in range(1,len(bounds)):
		minbin.append(bounds[i] - bounds[i-1])

	mincad = []
	for i in range(1,len(date)):
		mincad.append(date[i] - date[i-1])

#	if min(minbin) < max(mincad):
#		raise ValueError('Smallest bin must be larger than largest gap between cadences')

	bflux =[]

	#Iterate over the bins starting with bounds[1]
	for i in range(1,len(bounds)):
		#Initialise variables that need to be reet after each iteration
		bin = []
		t = []
		extrabitHigh = 0.
		extrabitLow = 0.
		eb1 = None
		eb2 = None

		#Iterate over the length of the flux array
		for j in range(len(flux)):
			#Catch all date bins within the bounds of the bin
			if date[j] >= bounds[i-1] and date[j] < bounds[i]:
				bin.append(flux[j])
				t.append(j)
		try:
			tmin = min(t)
			tmax = max(t)

			#We now have to deal with the bits either side of date[tmin]
			#and date[tmax]

		#First the bit between date[tmax] and bounds[i]
			if bounds[i] > date[tmax]+1e-6:
				if interpm == 'quad':
					extrabitHigh = integrate.quad(intpl,date[tmax],bounds[i])[0]
				elif interpm == 'romberg':
					extrabitHigh = integrate.romberg(intpl,date[tmax],bounds[i])

			#Now the bit between bound[i-1] and t[min]
			if bounds[i-1] < date[tmin]-1e-6:
				if interpm == 'quad':
					extrabitLow = integrate.quad(intpl,bounds[i-1],date[tmin])[0]
				elif interpm == 'romberg':
					extrabitLow = integrate.romberg(intpl,bounds[i-1],date[tmin])

			#Scale the extrabits as if they were the size of a full bin
			#Catch error if there is no extrabit
			try:
				eb1 = (extrabitHigh/(float(bounds[i]-date[tmax])))
				bin.append(eb1)
			except ZeroDivisionError:
				pass
			try:
				eb2 = (extrabitLow/(float(date[tmin]-bounds[i-1])))
				bin.append(eb2)
			except ZeroDivisionError:
				pass
		except ValueError:
			if interpm == 'quad':
				totbin = integrate.quad(intpl,bounds[i-1],bounds[i])[0]
			elif interpm == 'romberg':
				totbin =  integrate.romberg(intpl,bounds[i-1],bounds[i])
			eb = (totbin/(float(bounds[i]-bounds[i-1])))
			bin.append(eb)
		binflux = n.mean(bin)
		bflux.append(binflux)
	return bdate,bflux

def cutBadData(date,flux):
	date2 = date[n.logical_and(n.logical_and(n.isfinite(date),
	n.isfinite(flux)),flux != 0.0)]
	flux2 = flux[n.logical_and(n.logical_and(n.isfinite(date),
	n.isfinite(flux)),flux != 0.0)]
	return date2,flux2


def do_plot(date,flux,status=0):
	xmin = min(date)
	xmax = max(date)
	ymin = min(flux)
	ymax = max(flux)
	xr = xmax - xmin
	yr = ymax - ymin
	try:
		params = {'backend': 'png','axes.linewidth': 2.5,
			'axes.labelsize': 24,
			'axes.font': 'sans-serif',
			'axes.fontweight' : 'bold',
			'text.fontsize': 12,
			'legend.fontsize': 12,
			'xtick.labelsize': 16,
			'ytick.labelsize': 16}
		rcParams.update(params)
	except:
		print 'ERROR -- KEPCLIP: install latex for scientific plotting'
		status = 1

	if status == 0:
		plt.clf()

        plt.plot(date,flux,color='#0000ff',linestyle='-',linewidth=1.0)
	date = n.insert(date,[0],[date[0]]) 
	date = n.append(date,[date[-1]])
	flux = n.insert(flux,[0],[0.0]) 
	flux = n.append(flux,[0.0])
	plt.fill(date,flux,fc='#ffff00',linewidth=0.0,alpha=0.2)
	plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
	if ymin-yr*0.01 <= 0.0:
		plt.ylim(1.0e-10,ymax+yr*0.01)
	else:
		plt.ylim(ymin-yr*0.01,ymax+yr*0.01)
	xlab = 'BJD'
	ylab = 'e- / cadence'
	plt.xlabel(xlab, {'color' : 'k'})
	plt.ylabel(ylab, {'color' : 'k'})
        plt.ion()
	plt.grid()
	plt.ioff()



def kepbin(infile,outfile,fluxcol,do_nbin,nbins,do_binwidth,binwidth,
	do_ownbins,binfile,method,interpm,plot,clobber,verbose,logfile,status):
	"""
	Setup the kepbin environment
	"""
	# log the call 
	hashline = '----------------------------------------------------------------------------'
	kepmsg.log(logfile,hashline,verbose)
	call = 'KEPBIN -- '
	call += 'infile='+infile+' '
	call += 'outfile='+outfile+' '
	call += 'fluxcol='+fluxcol+ ' '
	donbin = 'n'
	if (do_nbin): donbin = 'y'
	call += 'donbin='+donbin+ ' '
	dobinwidth = 'n'
	if (do_binwidth): dobinwidth = 'y'
	call += 'dbinwidth='+dobinwidth+ ' '
	doownbin = 'n'
	if (do_ownbins): doownbin = 'y'
	call += 'doownbin='+doownbin+ ' '
	call += 'method='+method+' '
	call += 'interpm='+interpm+' '
	plotit = 'n'
	if (plot): plotit = 'y'
	call += 'plot='+plotit+ ' '
	overwrite = 'n'
	if (clobber): overwrite = 'y'
	call += 'clobber='+overwrite+ ' '
	chatter = 'n'
	if (verbose): chatter = 'y'
	call += 'verbose='+chatter+' '
	call += 'logfile='+logfile
	kepmsg.log(logfile,call+'\n',verbose)


	# start time
	kepmsg.clock('KEPCLIP started at',logfile,verbose)

	# test log file
	logfile = kepmsg.test(logfile)
	# clobber output file
	if clobber:
		status = kepio.clobber(outfile,logfile,verbose)
	if kepio.fileexists(outfile):
		message = 'ERROR -- KEPCLIP: ' + outfile + ' exists. Use --clobber'
		status = kepmsg.err(logfile,message,verbose)

	# open input file
	if status == 0:
		instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
		tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,
			infile,logfile,verbose,status)

	# fudge non-compliant FITS keywords with no values
	if status == 0:
		instr = kepkey.emptykeys(instr,file,logfile,verbose)

	# input data
	if status == 0:
		table = instr[1].data

	# read time and flux columns
	date = table.field('barytime')
	flux = table.field(fluxcol)
	#cut out infinites and zero flux columns
	date,flux = cutBadData(date,flux)

        if do_nbin:
		bdate,bflux = bin_funct(date,flux,nbins=nbins
			,method=method,interpm=interpm)
	elif do_binwidth:
		bdate,bflux = bin_funct(date,flux,binwidth=binwidth
			,method=method,interpm=interpm)
	elif do_ownbins:
		filepointer = open(binfile,'r')
		ownbins = []
		for line in filepointer:
			splitted = line.split()
			ownbins.append(float(splitted[0]))
		ownbins = n.array(ownbins)
		bdate,bflux = bin_funct(date,flux,ownbins=ownbins
			,method=method,interpm=interpm)
	if plot:
		do_plot(bdate,bflux)
	if status == 0:
		col1 = pyfits.Column(name='bdate',format='E',unit='day',array=bdate)
		col2 = pyfits.Column(name='bflux',format='E',unit='e-/cadence',array=bflux)
		cols = pyfits.ColDefs([col1,col2])
		instr.append(pyfits.new_table(cols))
		instr[-1].header.update('EXTNAME','BINNED DATA','extension name')
		instr.writeto(outfile)
	# close input file
	if status == 0:
		status = kepio.closefits(instr,logfile,verbose)

	# end time
	if status == 0:
		message = 'KEPBIN completed at'
	else:
		message = '\nKEPBIN aborted at'
	kepmsg.clock(message,logfile,verbose)


#def main():
parfile = iraf.osfn("kepler$kepbin.par")
t = iraf.IrafTaskFactory(taskname="kepbin", value=parfile, function=kepbin)







#if __name__ == "__main__":
#    main()
