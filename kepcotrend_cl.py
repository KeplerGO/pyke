
"""
Name: kepcotrend.py
Written by: Tom Barclay
Date released: July 27 2011

Changelog:
1.0 released
1.0.1 caught not having a new version of numpy with the in1d function and changed the simplex algorithm so that it fits the first two BV then fits
		the rest of them, now much more stable (Sep 5 2011)
2.0 code will now work on short cadence data by interpolating the basis vectors down to short cadence, several interpolation methods are available
3.0 made it possible to run from the command line if there are arguements given to the call to kepcotrend
"""

__svnid__ = "$Id: kepcotrend_cl.py 4005 2012-03-14 21:51:22Z tsbarcl2 $"
__url__ = "$URL: svn+ssh://mstill@murzim.amn.nasa.gov/data-repo/trunk/data/flight/go/PyKE/kepler/kepcotrend_cl.py $"

#some of the module imports are at the bottom because because I 
#import different modules depending on the way the program is called


import sys
import matplotlib.pyplot as plt
from matplotlib.cbook import is_numlike
from scipy.optimize import leastsq
from scipy.optimize import fmin as effmin
from scipy.interpolate import interp1d
from astropy.io import fits as pyfits
import sys
from numpy.linalg import lstsq, inv
from numpy import interp as interpolat
from numpy import *
import kepmsg, kepio, kepkey



def cutBadData(cad,date,flux,err):
	"""
	this function finds cadences with bad data in and removes them returning only cadences which contain good data
	"""
	date2 = date[logical_and(logical_and(isfinite(date),
	isfinite(flux)),flux != 0.0)]
	cad2 = cad[logical_and(logical_and(isfinite(date),
	isfinite(flux)),flux != 0.0)]
	flux2 = flux[logical_and(logical_and(isfinite(date),
	isfinite(flux)),flux != 0.0)]
	err2 = err[logical_and(logical_and(isfinite(date),
	isfinite(flux)),flux != 0.0)]
	bad_data = logical_and(logical_and(isfinite(date),
	isfinite(flux)),flux != 0.0)
	return cad2,date2,flux2,err2,bad_data

def putInNans(bad_data,flux):
	"""
	Function finds the cadences where the data has been removed using cutBadData() and puts data back in. The flux data put back in is nan.
	
	This function is used when writing data to a FITS files.
	
	bad_data == True means the datapoint is good!!
	"""
	newflux = zeros(len(bad_data))
	j = 0
	for i in range(len(bad_data)):
		if bad_data[i] == True:
			newflux[i] = flux[j]
			j += 1
		elif bad_data[i] == False:
			newflux[i] = nan
	return newflux

# def get_pcomp_list_newformat(pcompdata,pcomplist,newcad):
# 	pcomp = zeros((len(pcomplist),len(newcad)))
# 	for i in range(len(array(pcomplist))):
# 		j = int(array(pcomplist)[i])-1
# 		dat = pcompdata[...,j+4]
# 		pcomp[i] = dat[in1d(pcompdata[...,0],newcad)]
# 	return pcomp

def get_pcomp_list_newformat(bvdat,pcomplist,newcad,short,scinterp):
	"""
	Finds cotrending basis vectors which have been requested to be used by the user and adds them to an array.
	"""
	status = False
	pcomp = zeros((len(pcomplist),len(newcad)))
	for i in range(len(array(pcomplist))):
		j = int(array(pcomplist)[i])
		dat = bvdat.field('VECTOR_%s' %j)[isnan(bvdat.field('CADENCENO')) == False]
		bvcadfull = bvdat.field('CADENCENO')[isnan(bvdat.field('CADENCENO')) == False]
		#try:
		if short:
			#if the data is short cadence the interpolate the basis vectors
			bv_data = dat[in1d(bvdat.field('CADENCENO'),newcad)]
			bv_cad = bvcadfull[in1d(bvdat.field('CADENCENO'),newcad)]
			#funny things happen why I use interp1d for linear interpolation
			#so I have opted to use the numpy interp function for linear
			if scinterp == 'linear':
				intpl = interpolat(newcad,bv_cad,bv_data,left=bv_data[0],right=bv_data[-1])
				pcomp[i] = intpl
			else:
				intpl = interp1d(bv_cad,bv_data,kind=scinterp,bounds_error=False,fill_value=None)
				pcomp[i] = where(isnan(intpl(newcad)),0,intpl(newcad))
				mid_pt = floor(median(arange(len(pcomp[i]))))

				p_len = len(pcomp[i])

				lower = [logical_and(arange(p_len) < mid_pt,pcomp[i] == 0)]
				upper = [logical_and(arange(p_len) > mid_pt,pcomp[i] == 0)]
				pcomp[i][lower] = bv_data[0]
				pcomp[i][upper] = bv_data[-1]

			
		else:
			pcomp[i] = dat[in1d(bvdat.field('CADENCENO'),newcad)]
		#except NameError:
		#	status = True
		#	break
	return pcomp,status

def make_sc_lc(obs_cad,bv_cad,flux):
	"""
	make short cadence data look like long cadence data
	"""
	newflux = zeros(len(bv_cad))
	for i in range(len(bv_cad)):
		mask = logical_and(obs_cad > bv_cad[i] - 15,obs_cad < bv_cad[i] + 15)
		newflux[i] = obs_cad[mask]
	return newflux
		
	


def near_intpl(xout,xin,yin): 
	""" 
	Interpolate the curve defined by (xin, yin) at points xout. The array 
	xin must be monotonically increasing. The output has the same data type as 
	the input yin. 
 
	:param yin: y values of input curve 
	:param xin: x values of input curve 
	:param xout: x values of output interpolated curve 
	:param method: interpolation method ('linear' | 'nearest') 
 
  	@:rtype: numpy array with interpolated curve 
  	"""  
  	lenxin = len(xin)  
  
  	i1 = searchsorted(xin, xout)  
  	i1[ i1==0 ] = 1  
  	i1[ i1==lenxin ] = lenxin-1  
  
  	x0 = xin[i1-1]  
  	x1 = xin[i1]  
  	y0 = yin[i1-1]  
  	y1 = yin[i1]  
  
	return where(abs(xout - x0) < abs(xout - x1), y0, y1)  

def get_pcomp_list(pcompdata,pcomplist,newcad):
	pcomp = zeros((len(pcomplist),len(newcad)))
	for i in range(len(array(pcomplist))):
		j = int(array(pcomplist)[i])-1
		dat = pcompdata[...,j+2]
		pcomp[i] = dat[in1d(pcompdata[...,1],newcad)]
	return pcomp


def do_lsq_uhat(pcomps,cad,flux,orthog=True):
	"""
	does a linear least squares fit of the basis vectors to the light curve using the 'matrix' method - U(transpose) * y = coeffs
	In my implimentation y is a horizontal 1D array and U is also a long thin array of length correpsonding to the number of basis vectors use. In effect what I have is my U is already transposed and y need to be transposed to used. First I convert to what is expected in leasts squares fitting
	
	Orthog is a boolean saying if the basis vectors are orthogonal - they are not orthogonal if there has been masking or iterative fitting
	-> now changed to force the fitting to always be a generic least squares fit instead of relying on any orthogonality
	
	"""
	U_hat = matrix(pcomps).transpose()
	y_hat = matrix(flux).transpose()
	
	U_trans = U_hat.transpose()
	
	if orthog == True:
		coeffs = inv(U_trans * U_hat) * U_trans * y_hat

	elif orthog == False:
		coeffs = inv(U_trans * U_hat) * U_trans * y_hat

	coeffs = array(0. - coeffs)

	return coeffs

def do_lsq_nlin(pcomps,cad,flux):
	"""
	does a linear least squares fit of the basis vectors to the light curve using the 'lst_sq' method - this performs a Levenberg-Marquart least squares fit. The initial guess is an array of zeros
	"""
	guess = append(array([1.]),zeros(len(pcomps)-1))
	t = leastsq(fitfunct,guess,
		args=(pcomps,cad,flux),full_output=0)
	return -array(t[0])

def do_lsq_fmin(pcomps,cad,flux):
	"""
	performs a simplex fit of the basis vectors to the light curve.
	Initial guess is an array with 1. as the first element and zero as the value of all other elements in the array
	"""
	guess = append(array([1.]),zeros(len(pcomps)-1))
	t = effmin(fitfunct_fmin,guess,args=(pcomps,cad,flux))
	return -array(t)

def do_lsq_fmin_pow(pcomps,cad,flux,order):
	"""
	performs a simplex fit of the basis vectors to the light curve.
	Initial guess is an array with 1. as the first element and zero as the value of all other elements in the array
	"""
	guess = array([1,0])
	initial = effmin(fitfunct_fmin_pow,guess,args=(pcomps[0:2],cad,flux,order))
	#guess = append(array([1.]),zeros(len(pcomps)-1))
	guess = append(initial,zeros(len(pcomps)-2))
	t = effmin(fitfunct_fmin_pow,guess,args=(pcomps,cad,flux,order))
	return -array(t)

def fitfunct_fmin(scale,pcomp,date,zeroflux):
	"""
	the function called by the simplex fitting alogirthm
	"""
	outflux = copy(zeroflux)
	for i in range(shape(pcomp)[0]):
		outflux -= scale[i]*pcomp[i]
	sumsq = sum(abs(array(outflux)))
	#sumsq = sum(array(outflux)**2)
	return sumsq

def fitfunct_fmin_pow(scale,pcomp,date,zeroflux,order):
	"""
	the function called by the simplex fitting alogirthm
	"""
	outflux = copy(zeroflux)
	for i in range(shape(pcomp)[0]):
		outflux -= scale[i]*pcomp[i]
	sumsq = sum(power(abs(array(outflux)),order))
	#sumsq = sum(array(outflux)**2)
	return sumsq

	
def fitfunct(scale,pcomp,date,zeroflux):
	"""
	the function called by the least squares fitting algorithm
	"""
	outflux = copy(zeroflux)
	for i in range(shape(pcomp)[0]):
		outflux -= scale[i]*pcomp[i]
	return outflux

def get_newflux(oldflux,pcomps,s):
	"""
	uses the coefficients found by the fitting of the basis vectors to the light curve to correct the flux in the light curve
	Each basis vector is multiplid by a coefficient and then subtracted from the light curve
	"""
	newflux = copy(oldflux)
	for i in range(len(s)):
		newflux += s[i]*pcomps[i]
	return newflux

def get_pcompsum(pcomps,s):
	"""
	calculates the sum of basis vectors which are to be subtracted from the light curve to produce the corrected data.
	"""
	pcompsum = 0.
	for i in range(len(s)):
		pcompsum += s[i]*pcomps[i]
	return pcompsum

def chi2_gtf(obs,expect,err,dof):
	"""
	calculates a chi squared of the model fit to the data
	"""
	chisqu = 0.
	obs = obs 
	expect  = expect 
	err = err 
	for i in range(len(obs)):
		chisqu += ((1.0/(err[i]))*((obs[i] - expect[i])))**2
	chisqu = chisqu * (1.0/float(dof))
	return chisqu

def rms(O,E):
	"""
	calculates a root mean square of the model fit to the data
	"""
	rms= sqrt(sum((O-E)**2)/len(O))
	return rms
	
def do_lst_iter(bvs,cad,flux,nsigma,niter,method,order):
	"""
	performs an iterative fit of the basis vectors to the light curve
	after each fit outliers further than nsigma from the fit are removed and the fit recalculated.
	
	The sigma actually means a median absolute deviation from the median.
	"""
	iiter = 1
	fluxnew = copy(flux)
	lcnew = copy(cad)
	bvsnew = copy(bvs)
	if method == 'matrix':
		t = do_lsq_uhat(bvsnew,lcnew,fluxnew,False)
	elif method == 'lst_sq':
		t = do_lsq_nlin(bvsnew,lcnew,fluxnew)
	elif method == 'simplex':
		t = do_lsq_fmin_pow(bvsnew,lcnew,fluxnew,order)
	elif method == 'simplex_abs':
		t = do_lsq_fmin_pow(bvsnew,lcnew,fluxnew)
	elif method == 'llsq':
		t = do_lsq_uhat(bvsnew,lcnew,fluxnew,False)
	bvsum = get_pcompsum(bvsnew,t)
	while (iiter < niter):
		iiter += 1
		matchrange = 1.4826 * nsigma*MAD_model(subtract(fluxnew,bvsum))
		mask = abs(fluxnew - bvsum) < matchrange
		fluxnew = fluxnew[mask]
		lcnew = lcnew[mask]
		try:
			bvsnew = copy(bvsnew2)
		except:
			pass
		bvsnew2 = newpcompsarray(bvsnew,mask)
		#print shape(bvsnew),shape(bvsnew2),shape(mask[mask])
		for i in range(shape(bvsnew)[0]):
			bvsnew2[i] = bvsnew[i][mask]
		if method == 'matrix':
			t = do_lsq_uhat(bvsnew2,lcnew,fluxnew,False)
		elif method == 'lst_sq':
			t = do_lsq_nlin(bvsnew2,lcnew,fluxnew)
		elif method == 'simplex':
			t = do_lsq_fmin_pow(bvsnew2,lcnew,fluxnew,order)
		elif method == 'simplex_abs':
			t = do_lsq_fmin_pow(bvsnew2,lcnew,fluxnew)
		bvsum = get_pcompsum(bvsnew2,t)

	return t,mask

def newpcompsarray(pcomp,mask):
	pcompnew = zeros((shape(pcomp)[0],len(mask[mask])))
	return pcompnew


def MAD_model(xx,minSd=1E-16): 
      """Median Absolute Deviation 
      """ 
      absdev=abs(xx) 
      mad=median(absdev,0) 
      mad=maximum(mad,multiply(ones(mad.shape,float32), 
                                           (minSd/1.48))) 
      return mad 


def make_outfile(fitsfile,outfile,flux_new,bvsum,version):
	"""
	creates a fits file identical to the input fits file save from containing two extra columns - CBVSAP_MODL and CBVSAP_FLUX which are the sum of basis vectors fit to the data and the resulting corrected flux after the basis vector fit has been subtracted
	"""
	
	if version == 1:
		unit = 'e-/cadence'
		flux_new = flux_new * 1625.3514 #convert to e-/cadence
	elif version == 2:
		unit = 'e-/s'
	col1 = pyfits.Column(name='CBVSAP_MODL',
		format='E13.7   ',unit=unit,array=bvsum)
	col2 = pyfits.Column(name='CBVSAP_FLUX',format='E13.7   ',
		unit=unit,array=flux_new)
	cols = fitsfile[1].columns + col1 + col2
	fitsfile[1] = pyfits.new_table(cols,header=fitsfile[1].header)
	fitsfile.writeto(outfile)
	
def do_plot(date,flux_old,flux_new,bvsum,cad,bad_data,cad_nans,version):
	try:
		params = {'backend': 'eps',
			'axes.linewidth': 2.5,
			'axes.labelsize': 24,
			'axes.font': 'sans-serif',
			'axes.fontweight' : 'bold',
			'text.fontsize': 12,
			'legend.fontsize': 12,
			'xtick.labelsize': 16,
			'ytick.labelsize': 16}
		plt.rcParams.update(params)
	except:
		#print 'ERROR -- KEPCLIP: install latex for scientific plotting'
		pass
	plt.figure(figsize=[15,8])
	plt.clf()
	
	
	if version == 1:
		barytime0 = float(int(date[0] / 100) * 100.0)
		date_sub = date - barytime0
		xlab = r'BJD $-$ %d' %(barytime0+2400000.)
	elif version == 2:
		barytime0 = float(int((date[0]+54833.) / 100) * 100.0)
		date_sub = date+54833. - barytime0
		xlab = r'BJD $-$ %d' %(barytime0+2400000.)
	
	try:
		nrm1 = len(str(int(flux_old.max())))-1
	except:
		nrm1 = 0
	flux_old_sub = flux_old / 10**nrm1
	bvsum_sub = bvsum / 10**nrm1
	ylab1 = r'10$^%d$ e$^-$ s$^{-1}$' % nrm1
	
	try:
		nrm2 = len(str(int(flux_new.max())))-1
	except:
		nrm2 = 0
	flux_new_sub = flux_new / 10**nrm2
	ylab2 = r'10$^%d$ e$^-$ s$^{-1}$' % nrm2
	
	xmin = min(date_sub)
	xmax = max(date_sub)
	ymin1 = min(min(flux_old_sub),min(bvsum_sub))
	ymax1 = max(max(flux_old_sub),max(bvsum_sub))
	ymin2 = min(flux_new_sub)
	ymax2 = max(flux_new_sub)
	xr = xmax - xmin
	yr1 = ymax1 - ymin1
	yr2 = ymax2 - ymin2
	
	#plt.axes([0.06,0.1,0.93,0.88])
	ax1 = plt.subplot(211)
	
	blocks = split_on_nans(bad_data,cad_nans)
	for i in range(len(blocks)):
		if i == 0:
			block = [blocks[0],blocks[i]]
		else:
			block = [blocks[i-1],blocks[i]]
		mask = logical_and(cad >= block[0],cad <= block[1])
		plot_x = date_sub[mask]
		plot_y = flux_old_sub[mask]
		if nan in plot_y:
			break
		plt.plot(plot_x,plot_y,color='#0000ff',linestyle='-',linewidth=1.0)
		plot_y = bvsum_sub[mask]
		plt.plot(plot_x,plot_y,color='red',linestyle='-',
		linewidth=2.0)
	date2 = insert(date_sub,[0],[date_sub[0]]) 
	date2 = append(date2,[date_sub[-1]])
	flux2 = insert(flux_old_sub,[0],[0.0]) 
	flux2 = append(flux2,[0.0])
	plt.fill(date2,flux2,fc='#FFFACD',linewidth=0.0)
	plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
	if ymin1-yr1*0.01 <= 0.0:
		plt.ylim(1.0e-10,ymax1+yr1*0.01)
	else:
		plt.ylim(ymin1-yr1*0.01,ymax1+yr1*0.01)
	plt.xlabel(xlab, {'color' : 'k'})
	plt.ylabel(ylab1, {'color' : 'k'})
	#plt.ion()
	plt.grid()
	
	
	ax2 = plt.subplot(212,sharex=ax1)
	
	for i in range(len(blocks)):
		if i == 0:
			block = [blocks[0],blocks[i]]
		else:
			block = [blocks[i-1],blocks[i]]
		mask = logical_and(cad >= block[0],cad <= block[1])
		plot_x = date_sub[mask]
		plot_y = flux_new_sub[mask]
		if nan in plot_y:
			break
		plt.plot(plot_x,plot_y,color='#0000ff',linestyle='-',linewidth=1.0)
		plot_y = bvsum_sub[mask]
	date2 = insert(date_sub,[0],[date_sub[0]]) 
	date2 = append(date2,[date_sub[-1]])
	flux2 = insert(flux_new_sub,[0],[0.0]) 
	flux2 = append(flux2,[0.0])
	plt.fill(date2,flux2,fc='#FFFACD',linewidth=0.0)
	plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
	if ymin2-yr2*0.01 <= 0.0:
		plt.ylim(1.0e-10,ymax2+yr2*0.01)
	else:
		plt.ylim(ymin2-yr2*0.01,ymax2+yr2*0.01)
	plt.xlabel(xlab, {'color' : 'k'})
	plt.ylabel(ylab2, {'color' : 'k'})
	plt.grid()
	
	plt.subplots_adjust(0.1,0.1,0.94,0.94,0.0,0.0)
	plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
	plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
	
	plt.draw()


def split_on_nans(bad_data,cad):
	blocks = []
	time_of_nans = cad[bad_data == False]
	if bad_data[0] == True:
		blocks.append(cad[0])
	for i in range(1,len(time_of_nans)):
		if time_of_nans[i] - time_of_nans[i-1] > 1:
			blocks.append(time_of_nans[i])
	if bad_data[-1] == True:
		blocks.append(cad[-1])
	return blocks
		

def kepcotrendsc(infile,outfile,bvfile,listbv,fitmethod,fitpower,iterate,sigma,maskfile,scinterp,plot,clobber,verbose,logfile,status):
	"""
	Setup the kepcotrend environment
	
	infile: 
	the input file in the FITS format obtained from MAST
	
	outfile:
	The output will be a fits file in the same style as the input file but with two additional columns: CBVSAP_MODL and CBVSAP_FLUX. The first of these is the best fitting linear combination of basis vectors. The second is the new flux with the basis vector sum subtracted. This is the new flux value. 
	
	plot:
	either True or False if you want to see a plot of the light curve
	The top plot shows the original light curve in blue and the sum of basis vectors in red
	The bottom plot has had the basis vector sum subracted
	
	bvfile:
	the name of the FITS file containing the basis vectors

	listbv:
	the basis vectors to fit to the data
	
	fitmethod:
	fit using either the 'llsq' or the 'simplex' method. 'llsq' is usually the correct one to use because as the basis vectors are orthogonal. Simplex gives you option of using a different merit function - ie. you can minimise the least absolute residual instead of the least squares which weights outliers less
	
	fitpower:
	if using a simplex you can chose your own power in the metir function - i.e. the merit function minimises abs(Obs - Mod)^P. P=2 is least squares, P = 1 minimises least absolutes
	
	iterate:
	should the program fit the basis vectors to the light curve data then remove data points further than 'sigma' from the fit and then refit
	
	maskfile:
	this is the name of a mask file which can be used to define regions of the flux time series to exclude from the fit. The easiest way to create this is by using keprange from the PyKE set of tools. You can also make this yourself with two BJDs on each line in the file specifying the beginning and ending date of the region to exclude.
	
	scinterp:
	the basis vectors are only calculated for long cadence data, therefore if you want to use short cadence data you have to interpolate the basis vectors. There are several methods to do this, the best of these probably being nearest which picks the value of the nearest long cadence data point.
	The options available are None|linear|nearest|zero|slinear|quadratic|cubic
	If you are using short cadence data don't choose none
	"""
	# log the call 
	hashline = '----------------------------------------------------------------------------'
	kepmsg.log(logfile,hashline,verbose)
	call = 'KEPCOTREND -- '
	call += 'infile='+infile+' '
	call += 'outfile='+outfile+' '
	call += 'bvfile='+bvfile+' '
#	call += 'numpcomp= '+str(numpcomp)+' '
	call += 'listbv= '+str(listbv)+' '
	call += 'fitmethod=' +str(fitmethod)+ ' '
	call += 'fitpower=' + str(fitpower)+ ' '
	iterateit = 'n'
	if (iterate): iterateit = 'y'
	call += 'iterate='+iterateit+ ' '
	call += 'sigma_clip='+str(sigma)+' '
	call += 'mask_file='+maskfile+' '
	call += 'scinterp=' + str(scinterp)+ ' '
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
	kepmsg.clock('KEPCOTREND started at',logfile,verbose)

	# test log file
	logfile = kepmsg.test(logfile)
    
	# clobber output file
	if clobber:
		status = kepio.clobber(outfile,logfile,verbose)
	if kepio.fileexists(outfile): 
		message = 'ERROR -- KEPCOTREND: ' + outfile + ' exists. Use --clobber'
		status = kepmsg.err(logfile,message,verbose)
	
	# open input file
	if status == 0:
		instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
		tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,
			infile,logfile,verbose,status)

	# fudge non-compliant FITS keywords with no values
	if status == 0:
		instr = kepkey.emptykeys(instr,file,logfile,verbose)

	if status == 0:
		if not kepio.fileexists(bvfile):
			message = 'ERROR -- KEPCOTREND: ' + bvfile + ' does not exist.'
			status = kepmsg.err(logfile,message,verbose)
	
	#lsq_sq - nonlinear least squares fitting and simplex_abs have been removed from the option in PyRAF but they are still in the code!
	if status == 0:
		if fitmethod not in ['llsq','matrix','lst_sq','simplex_abs','simplex']:
			message = 'Fit method must either: llsq, matrix, lst_sq or simplex'
			status = kepmsg.err(logfile,message,verbose)
	
	if status == 0:
		if not is_numlike(fitpower) and fitpower is not None:
			message = 'Fit power must be an real number or None'
			status = kepmsg.err(logfile,message,verbose)
	
	
		
	if status == 0:	
		if fitpower is None:
			fitpower = 1.

	# input data
	if status == 0:
		short = False
		try:
			test = str(instr[0].header['FILEVER'])
			version = 2
		except KeyError:
			version = 1

		table = instr[1].data
		if version == 1:
			if str(instr[1].header['DATATYPE']) == 'long cadence':
				#print 'Light curve was taken in Lond Cadence mode!'
				quarter = str(instr[1].header['QUARTER'])
				module = str(instr[1].header['MODULE'])
				output = str(instr[1].header['OUTPUT'])
				channel = str(instr[1].header['CHANNEL'])
	
				lc_cad_o = table.field('cadence_number')
				lc_date_o = table.field('barytime')
				lc_flux_o = table.field('ap_raw_flux') / 1625.3468 #convert to e-/s
				lc_err_o = table.field('ap_raw_err') / 1625.3468 #convert to e-/s
			elif str(instr[1].header['DATATYPE']) == 'short cadence':
				short = True
				#print 'Light curve was taken in Short Cadence mode!'
				quarter = str(instr[1].header['QUARTER'])
				module = str(instr[1].header['MODULE'])
				output = str(instr[1].header['OUTPUT'])
				channel = str(instr[1].header['CHANNEL'])
	
				lc_cad_o = table.field('cadence_number')
				lc_date_o = table.field('barytime')
				lc_flux_o = table.field('ap_raw_flux') / 54.178 #convert to e-/s
				lc_err_o = table.field('ap_raw_err') / 54.178 #convert to e-/s

		elif version == 2:
			if str(instr[0].header['OBSMODE']) == 'long cadence':
				#print 'Light curve was taken in Long Cadence mode!'

				quarter = str(instr[0].header['QUARTER'])
				module = str(instr[0].header['MODULE'])
				output = str(instr[0].header['OUTPUT'])
				channel = str(instr[0].header['CHANNEL'])
	
				lc_cad_o = table.field('CADENCENO')
				lc_date_o = table.field('TIME')
				lc_flux_o = table.field('SAP_FLUX')
				lc_err_o = table.field('SAP_FLUX_ERR')
			elif str(instr[0].header['OBSMODE']) == 'short cadence':
				#print 'Light curve was taken in Short Cadence mode!'
				short = True
				quarter = str(instr[0].header['QUARTER'])
				module = str(instr[0].header['MODULE'])
				output = str(instr[0].header['OUTPUT'])
				channel = str(instr[0].header['CHANNEL'])
	
				lc_cad_o = table.field('CADENCENO')
				lc_date_o = table.field('TIME')
				lc_flux_o = table.field('SAP_FLUX')
				lc_err_o = table.field('SAP_FLUX_ERR')
		

		if str(quarter) == str(4) and version == 1:
			lc_cad_o = lc_cad_o[lc_cad_o >= 11914]
			lc_date_o = lc_date_o[lc_cad_o >= 11914]
			lc_flux_o = lc_flux_o[lc_cad_o >= 11914]
			lc_err_o = lc_err_o[lc_cad_o >= 11914]

		# bvfilename = '%s/Q%s_%s_%s_map.txt' %(bvfile,quarter,module,output)
		# if str(quarter) == str(5):
		# 	bvdata = genfromtxt(bvfilename)
		# elif str(quarter) == str(3) or str(quarter) == str(4):
		# 	bvdata = genfromtxt(bvfilename,skip_header=22)
		# elif str(quarter) == str(1):
		# 	bvdata = genfromtxt(bvfilename,skip_header=10)
		# else:
		# 	bvdata = genfromtxt(bvfilename,skip_header=13)
		
		if short and scinterp == 'None':
			message = 'You cannot select None as the interpolation method because you are using short cadence data and therefore must use some form of interpolation. I reccommend nearest if you are unsure.'
			status = kepmsg.err(logfile,message,verbose)
		
		bvfiledata = pyfits.open(bvfile)
		bvdata = bvfiledata['MODOUT_%s_%s' %(module,output)].data 
		
		
		if int(bvfiledata[0].header['QUARTER']) != int(quarter):
			message = 'CBV file and light curve file are from different quarters. CBV file is from Q%s and light curve is from Q%s' %(int(bvfiledata[0].header['QUARTER']),int(quarter))
			status = kepmsg.err(logfile,message,verbose)
	
	if status == 0:
		if int(quarter) == 4 and int(module) == 3:
			message = 'Approximately twenty days into Q4 Module 3 failed. As a result, Q4 light curves contain these 20 day of data. However, we do not calculate CBVs for this section of data.'
			status = kepmsg.err(logfile,message,verbose)
	
	if status == 0:
		

		#cut out infinites and zero flux columns
		lc_cad,lc_date,lc_flux,lc_err,bad_data = cutBadData(lc_cad_o,
			lc_date_o,lc_flux_o,lc_err_o)

		#get a list of basis vectors to use from the list given
		#accept different seperators
		listbv = listbv.strip()
		if listbv[1] in [' ',',',':',';','|',', ']:
			separator = str(listbv)[1]
		else:
			message = 'You must separate your basis vector numbers to use with \' \' \',\' \':\' \';\' or \'|\' and the first basis vector to use must be between 1 and 9' 
			status = kepmsg.err(logfile,message,verbose)
	
	
	if status == 0:
		bvlist = fromstring(listbv,dtype=int,sep=separator)

		if bvlist[0] == 0:
			message = 'Must use at least one basis vector' 
			status = kepmsg.err(logfile,message,verbose)
	if status == 0:
		#pcomps = get_pcomp(pcompdata,n_comps,lc_cad)
		# if str(quarter) == str(5):
		# 	bvectors = get_pcomp_list(bvdata,bvlist,lc_cad)
		# else:
		#	bvectors = get_pcomp_list_newformat(bvdata,bvlist,lc_cad)
		
		if short:
			bvdata.field('CADENCENO')[:] = (((bvdata.field('CADENCENO')[:] + (7.5/15.) )* 30.) - 11540.).round()
		
		bvectors,in1derror = get_pcomp_list_newformat(bvdata,bvlist,lc_cad,short,scinterp)
		
		if in1derror:
			message = 'It seems that you have an old version of numpy which does not have the in1d function included. Please update your version of numpy to a version 1.4.0 or later'
			status = kepmsg.err(logfile,message,verbose)
	if status == 0:
		
		medflux = median(lc_flux)
		n_flux = (lc_flux /medflux)-1
		n_err = sqrt(pow(lc_err,2)/ pow(medflux,2))
		
		#plt.errorbar(lc_cad,n_flux,yerr=n_err)
		#plt.errorbar(lc_cad,lc_flux,yerr=lc_err)
		
		#n_err = median(lc_err/lc_flux) * n_flux
		#print n_err
		
		#does an iterative least squares fit
		#t1 = do_leastsq(pcomps,lc_cad,n_flux)
		#

		if maskfile != '':
			domasking = True
			if not kepio.fileexists(maskfile):
				message = 'Maskfile %s does not exist' %maskfile
				status = kepmsg.err(logfile,message,verbose)
		else:
			domasking = False
		
		
			
	if status == 0:	
		if domasking:

			lc_date_masked = copy(lc_date)
			n_flux_masked = copy(n_flux)
			lc_cad_masked = copy(lc_cad)
			n_err_masked = copy(n_err)
			maskdata = atleast_2d(genfromtxt(maskfile,delimiter=','))
			#make a mask of True values incase there are not regions in maskfile to exclude.
			mask = zeros(len(lc_date_masked)) == 0.
			for maskrange in maskdata:
				if version == 1:
					start = maskrange[0] - 2400000.0
					end = maskrange[1] - 2400000.0
				elif version == 2:
					start = maskrange[0] - 2454833.
					end = maskrange[1] - 2454833.
				masknew = logical_xor(lc_date < start,lc_date > end)
				mask = logical_and(mask,masknew)

			lc_date_masked = lc_date_masked[mask]
			n_flux_masked = n_flux_masked[mask]
			lc_cad_masked = lc_cad_masked[mask]
			n_err_masked = n_err_masked[mask]
		else:
			lc_date_masked = copy(lc_date)
			n_flux_masked = copy(n_flux)
			lc_cad_masked = copy(lc_cad)
			n_err_masked = copy(n_err)

	
		#pcomps = get_pcomp(pcompdata,n_comps,lc_cad)

		bvectors_masked,hasin1d = get_pcomp_list_newformat(bvdata,bvlist,lc_cad_masked,short,scinterp)
		

		if (iterate) and sigma is None:
			message = 'If fitting iteratively you must specify a clipping range'
			status = kepmsg.err(logfile,message,verbose)
			
	if status == 0:
		#uses Pvals = yhat * U_transpose
		if (iterate):
			coeffs,fittedmask = do_lst_iter(bvectors_masked,lc_cad_masked
				,n_flux_masked,sigma,50.,fitmethod,fitpower)
		else:
			if fitmethod == 'matrix' and domasking:
				coeffs = do_lsq_uhat(bvectors_masked,lc_cad_masked,n_flux_masked,False)
			if fitmethod == 'llsq' and domasking:
				coeffs = do_lsq_uhat(bvectors_masked,lc_cad_masked,n_flux_masked,False)
			elif fitmethod == 'lst_sq':
				coeffs = do_lsq_nlin(bvectors_masked,lc_cad_masked,n_flux_masked)
			elif fitmethod == 'simplex_abs':
				coeffs = do_lsq_fmin(bvectors_masked,lc_cad_masked,n_flux_masked)
			elif fitmethod == 'simplex':
				coeffs = do_lsq_fmin_pow(bvectors_masked,lc_cad_masked,n_flux_masked,fitpower)
			else:
				coeffs = do_lsq_uhat(bvectors_masked,lc_cad_masked,n_flux_masked)
		
		
		
		flux_after = (get_newflux(n_flux,bvectors,coeffs) +1) * medflux
		flux_after_masked = (get_newflux(n_flux_masked,bvectors_masked,coeffs) +1) * medflux
		bvsum = get_pcompsum(bvectors,coeffs)
		
		bvsum_masked =  get_pcompsum(bvectors_masked,coeffs)
		
		#print 'chi2: ' + str(chi2_gtf(n_flux,bvsum,n_err,2.*len(n_flux)-2))
		#print 'rms: ' + str(rms(n_flux,bvsum))
		

		bvsum_nans = putInNans(bad_data,bvsum)
		flux_after_nans = putInNans(bad_data,flux_after)
		

	if plot and status == 0:
		bvsum_un_norm = medflux*(1-bvsum)
		#bvsum_un_norm = 0-bvsum
		#lc_flux = n_flux
		do_plot(lc_date,lc_flux,flux_after,
			bvsum_un_norm,lc_cad,bad_data,lc_cad_o,version)
		
	if status== 0:
		make_outfile(instr,outfile,flux_after_nans,bvsum_nans,version)
	
	# close input file
	if status == 0:
		status = kepio.closefits(instr,logfile,verbose)	   
		
		#print some results to screen:
		print '      -----      '
		if iterate:
			flux_fit = n_flux_masked[fittedmask]
			sum_fit = bvsum_masked[fittedmask]
			err_fit = n_err_masked[fittedmask]
		else:
			flux_fit = n_flux_masked
			sum_fit = bvsum_masked
			err_fit = n_err_masked
		print 'reduced chi2: ' + str(chi2_gtf(flux_fit,sum_fit,err_fit,len(flux_fit)-len(coeffs)))
		print 'rms: ' + str(medflux*rms(flux_fit,sum_fit))
		for i in range(len(coeffs)):
			print 'Coefficient of CBV #%s: %s' %(i+1,coeffs[i])
		print '      -----      '
	

	# end time
	if (status == 0):
		message = 'KEPCOTREND completed at'
	else:
		message = '\nKEPCOTTREND aborted at'
	kepmsg.clock(message,logfile,verbose)
	
	return


if len(sys.argv[1:]) > 0:
	import argparse

	parser = argparse.ArgumentParser(description='Remove systematic trends in photometry using cotrending basis vectors')
	#parser.add_argument('--infile', '-i', help='Name of input file', dest='infile',type=str,required=True)
	parser.add_argument('infile', help='Name of input file', type=str)
	#parser.add_argument('--outfile', '-o', help='Name of FITS file to output', dest='outfile',type=str,required=True)
	parser.add_argument('outfile', help='Name of FITS file to output', type=str)
	#parser.add_argument('--cbvfile', '-c', help='Name of file containing the CBVs',dest='bvfile',type=str,required=True)
	parser.add_argument('cbvfile', help='Name of file containing the CBVs', type=str)
	#parser.add_argument('--vectors', '-v', help='The CBVs to use', dest='listbv',type=str,required=True)
	parser.add_argument('--vectors', dest='listbv', help='The CBVs to use', type=str)

	parser.add_argument('--method', '-m', help='Fitting method',default='llsq',dest='fitmethod',type=str,choices=['llsq','simplex','lst_sq'])
	parser.add_argument('--fitpower', '-f', help='The index of the merit function (simplex only)', default=1, type=float, dest='fitpower')
	parser.add_argument('--iterate', action='store_true', help='Fit iteratively ', dest='iterate')
	parser.add_argument('--sigmaclip', '-s', help='Sigma clip value when iteratively fitting', default=False, dest='sigma', type=float)
	parser.add_argument('--maskfile', '-q', help='Name of file containing a mask', default='', dest='maskfile', type=str)
	parser.add_argument('--scinterp', '-w', help='Short cadence interpolation method', default='None', dest='scinterp', type=str,choices=['linear','nearest','slinear','quadratic','cubic'])
	parser.add_argument('--plot', '-p', action='store_true', help='Plot result?', dest='plot')
	parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
	parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
	parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
	parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


	args = parser.parse_args()

	kepcotrendsc(args.infile,args.outfile,args.bvfile,args.listbv, args.fitmethod, args.fitpower, args.iterate, args.sigma, args.maskfile, args.scinterp, args.plot, args.clobber, args.verbose, args.logfile, args.status)

else:
	from pyraf import iraf

        parfile = iraf.osfn("kepler$kepcotrend.par")
	t = iraf.IrafTaskFactory(taskname="kepcotrend", value=parfile, function=kepcotrendsc)
