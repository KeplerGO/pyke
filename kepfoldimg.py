from pyraf import iraf
import numpy, sys, time, pylab
from numpy import *
from astropy.io import fits as pyfits
from pylab import *
from matplotlib import *
import kepio, kepmsg, kepkey, kepstat, kepfit

# global variables

def kepfoldimg(infile,outfile,datacol,period,phasezero,binmethod,threshold,niter,nbins,
               plot,plotlab,clobber,verbose,logfile,status):

# startup parameters

    status = 0
    labelsize = 24; ticksize = 16; xsize = 17; ysize = 7
    lcolor = '#0000ff'; lwidth = 1.0; fcolor = '#ffff00'; falpha = 0.2

# log the call

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFOLD -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+datacol+' '
    call += 'period='+str(period)+' '
    call += 'phasezero='+str(phasezero)+' '
    call += 'binmethod='+binmethod+' '
    call += 'threshold='+str(threshold)+' '
    call += 'niter='+str(niter)+' '
    call += 'nbins='+str(nbins)+' '
    plotres = 'n'
    if (plot): plotres = 'y'
    call += 'plot='+plotres+ ' '
    call += 'plotlab='+plotlab+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPFOLDIMG started at: ',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPFOLDIMG: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,infile,logfile,verbose)

# input data

    if status == 0:
        table = instr[1].data
        incards = instr[1].header.cards
        indata, status = kepio.readfitscol(infile,table,datacol,logfile,verbose)
        barytime, status = kepio.readtimecol(infile,table,logfile,verbose)

# filter out NaNs

    work1 = []; work2 = []
    if status == 0:
        for i in range(len(barytime)):
            if (numpy.isfinite(barytime[i]) and
                numpy.isfinite(indata[i]) and indata[i] != 0.0):
                work1.append(barytime[i])
                work2.append(indata[i])
        barytime = array(work1,dtype='float64')
        indata = array(work2,dtype='float32')

# calculate phase

    if status == 0:
        phase2 = []
        phase1 = (barytime - phasezero) / period
        for i in range(len(phase1)):
            phase2.append(phase1[i] - int(phase1[i]))
            if phase2[-1] < 0.0: phase2[-1] += 1.0
        phase2 = array(phase2,'float32')

# sort phases

    if status == 0:
        ptuple = []
        phase3 = []
        data3 = []
        for i in range(len(phase2)):
            ptuple.append([phase2[i], indata[i]])
        phsort = sorted(ptuple,key=lambda ph: ph[0])
        for i in range(len(phsort)):
            phase3.append(phsort[i][0])
            data3.append(phsort[i][1])
        phase3 = array(phase3,'float32')
        data3 = array(data3,'float32')

# bin phases

    if status == 0:
        work1 = array([data3[0]],'float32')
        phase4 = array([],'float32')
        data4 = array([],'float32')
        dt = (phase3[-1] - phase3[0]) / nbins
        nb = 0.0
        for i in range(len(phase3)):
            if phase3[i] < phase3[0] + nb * dt or phase3[i] >= phase3[0] + (nb + 1.0) * dt:
                if len(work1) > 0:
                    phase4 = append(phase4,phase3[0] + (nb + 0.5) * dt)
                    if (binmethod == 'mean'):
                        data4 = append(data4,kepstat.mean(work1))
                    elif (binmethod == 'median'):
                        data4 = append(data4,kepstat.median(work1,logfile))
                    else:
                        coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
                            kepfit.lsqclip('poly0',[1.0],arange(0.0,float(len(work1)),1.0),work1,None,
                                           threshold,threshold,niter,logfile,verbose)
                        data4 = append(data4,coeffs[0])
                work1 = array([],'float32')
                nb += 1.0
            else:
                work1 = append(work1,data3[i])

# update HDU1 for output file

    if status == 0:
        cols = (instr[1].columns + ColDefs([Column(name='PHASE',format='E',array=phase1)]))
        instr[1] = pyfits.new_table(cols)
        instr[1].header.cards['TTYPE20'].comment = 'column title: phase'
        instr[1].header.cards['TFORM20'].comment = 'data type: float32'
        for i in range(len(incards)):
            if incards[i].key not in instr[1].header.keys():
                instr[1].header.update(incards[i].key, incards[i].value, incards[i].comment)
            else:
                instr[1].header.cards[incards[i].key].comment = incards[i].comment
        instr[1].header.update('PERIOD',period,'period defining the phase [d]')
        instr[1].header.update('BJD0',phasezero,'time of phase zero [BJD]')

# write new phased data extension for output file

    if status == 0:
        col1 = Column(name='PHASE',format='E',array=phase4)
        col2 = Column(name=datacol,format='E',unit='e/s',array=data4/cadence)
        cols = ColDefs([col1,col2])
        instr.append(new_table(cols))
        instr[-1].header.cards['TTYPE1'].comment = 'column title: phase'
        instr[-1].header.cards['TTYPE2'].comment = 'column title: simple aperture photometry'
        instr[-1].header.cards['TFORM1'].comment = 'column type: float32'
        instr[-1].header.cards['TFORM2'].comment = 'column type: float32'
        instr[-1].header.cards['TUNIT2'].comment = 'column units: electrons per second'
        instr[-1].header.update('EXTNAME','FOLDED','extension name')
        instr[-1].header.update('PERIOD',period,'period defining the phase [d]')
        instr[-1].header.update('BJD0',phasezero,'time of phase zero [BJD]')
        instr[-1].header.update('BINMETHD',binmethod,'phase binning method')
        if binmethod =='sigclip':
            instr[-1].header.update('THRSHOLD',threshold,'sigma-clipping threshold [sigma]')
            instr[-1].header.update('NITER',niter,'max number of sigma-clipping iterations')
    
# history keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)
        instr.writeto(outfile)

# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

# clean up x-axis unit

    if status == 0:
        ptime = array([],'float32')
        pout = array([],'float32')
        work = data4
        for i in range(len(phase4)):
            if (phase4[i] > 0.5): 
                ptime = append(ptime,phase4[i] - 1.0)
                pout = append(pout,work[i] / cadence)
        ptime = append(ptime,phase4)
        pout = append(pout,work / cadence)
        for i in range(len(phase4)):
            if (phase4[i] <= 0.5): 
                ptime = append(ptime,phase4[i] + 1.0)
                pout = append(pout,work[i] / cadence)
	xlab = 'Phase ($\phi$)'

# clean up y-axis units

    if status == 0:
	nrm = len(str(int(pout.max())))-1
	pout = pout / 10**nrm
	ylab = '10$^%d$ %s' % (nrm, plotlab)

# data limits

	xmin = ptime.min()
	xmax = ptime.max()
	ymin = pout.min()
	ymax = pout.max()
	xr = xmax - xmin
	yr = ymax - ymin
        ptime = insert(ptime,[0],[ptime[0]])
        ptime = append(ptime,[ptime[-1]])
        pout = insert(pout,[0],[0.0])
        pout = append(pout,0.0)

# plot new light curve

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
           pylab.rcParams.update(params)
        except:
            print 'ERROR -- KEPFOLD: install latex for scientific plotting'
            status = 1
    if status == 0 and plot:
	pylab.figure(1,figsize=[17,7])
        pylab.clf()
        pylab.axes([0.06,0.1,0.93,0.87])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.plot(ptime,pout,color=lcolor,linestyle='-',linewidth=lwidth)
        fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
	xlabel(xlab, {'color' : 'k'})
	ylabel(ylab, {'color' : 'k'})
        xlim(-0.49999,1.49999)
        if ymin >= 0.0:
            ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            ylim(1.0e-10,ymax+yr*0.01)
        pylab.grid()
        pylab.draw()

# stop time

    kepmsg.clock('KEPFOLDIMG ended at: ',logfile,verbose)

# --------------------------------------------
# main

parfile = iraf.osfn("kepler$kepfoldimg.par")
t = iraf.IrafTaskFactory(taskname="kepfoldimg", value=parfile, function=kepfoldimg)

