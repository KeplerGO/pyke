import sys
import os
import scipy
import numpy as np
from matplotlib import pyplot as plt
from . import kepmsg, kepio, kepkey, kepplot, kepfit, kepfunc

def kepsff(infile, outfile, datacol, cenmethod, stepsize, npoly_cxcy,
           sigma_cxcy, npoly_ardx, npoly_dsdt, sigma_dsdt, npoly_arfl,
           sigma_arfl, plotres, clobber, verbose, logfile):

    # log the call
    hashline = '--------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = ('KEPSFF -- '
            'infile={}'.format(infile)
            'outfile={}'.format(outfile)
            'datacol={}'.format(datacol)
            'cenmethod={}'.format(cenmethod)
            'stepsize={}'.format(stepsize)
            'npoly_cxcy={}'.format(npoly_cxcy)
            'sigma_cxcy={}'.format(sigma_cxcy)
            'npoly_ardx={}'.format(npoly_ardx)
            'npoly_dsdt={}'.format(npoly_dsdt)
            'sigma_dsdt={}'.format(sigma_dsdt)
            'npoly_arfl={}'.format(npoly_arfl)
            'sigma_arfl={}'.format(sigma_arfl)
            'plotres={}'.format(savep)
            'clobber={}'.format(clobber)
            'verbose={}'.format(verbose)
            'logfile={}'.format(logfile))

    kepmsg.log(logfile, call+'\n', verbose)
    # start time
    kepmsg.clock('KEPSFF started at', logfile, verbose)

    # clobber output file
    if clobber:
        kepio.clobber(outfile, logfile, verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPSFF: {} exists. Use clobber=True'.format(outfile)
        kepmsg.err(logfile, message, verbose)
    # open input file
    instr = pyfits.open(infile, 'readonly')
    tstart, tstop, bjdref, cadence = kepio.timekeys(instr, infile,
                                                    logfile, verbose)
    try:
        work = instr[0].header['FILEVER']
        cadenom = 1.0
    except:
        cadenom = cadence
    # fudge non-compliant FITS keywords with no values
    instr = kepkey.emptykeys(instr, infile, logfile, verbose)
    # read table structure
    table = kepio.readfitstab(infile, instr[1], logfile, verbose)
    # determine sequence of windows in time
    frametim = instr[1].header['FRAMETIM']
    num_frm = instr[1].header['NUM_FRM']
    exptime = frametim * num_frm / 86400
    tstart = table.field('TIME')[0]
    tstop = table.field('TIME')[-1]
    winedge = np.arange(tstart, tstop, stepsize)
    if tstop > winedge[-1] + stepsize / 2:
        winedge = np.append(winedge, tstop)
    else:
        winedge[-1] = tstop
    winedge = (winedge - tstart) / exptime
    winedge = winedge.astype(int)
    if len(table.field('TIME')) > winedge[-1] + 1:
        winedge = append(winedge,len(table.field('TIME')))
    elif len(table.field('TIME')) < winedge[-1]:
        winedge[-1] = len(table.field('TIME'))
    # step through the time windows
    for iw in range(1,len(winedge)):
        t1 = winedge[iw-1]
        t2 = winedge[iw]
        # filter input data table
        work1 = numpy.array([table.field('TIME')[t1:t2],
                             table.field('CADENCENO')[t1:t2],
                             table.field(datacol)[t1:t2],
                             table.field('MOM_CENTR1')[t1:t2],
                             table.field('MOM_CENTR2')[t1:t2],
                             table.field('PSF_CENTR1')[t1:t2],
                             table.field('PSF_CENTR2')[t1:t2],
                             table.field('SAP_QUALITY')[t1:t2]],'float64')
        work1 = np.rot90(work1, 3)
        work2 = work1[~np.isnan(work1).any(1)]
        work2 = work2[(work2[:,0] == 0.0) | (work2[:, 0] > 1e5)]

        # assign table columns
        intime = work2[:,7] + bjdref
        cadenceno = work2[:, 6].astype(int)
        indata = work2[:,5]
        mom_centr1 = work2[:, 4]
        mom_centr2 = work2[:, 3]
        psf_centr1 = work2[:, 2]
        psf_centr2 = work2[:, 1]
        sap_quality = work2[:, 0]
        if cenmethod == 'moments':
            centr1 = np.copy(mom_centr1)
            centr2 = np.copy(mom_centr2)
        else:
            centr1 = np.copy(psf_centr1)
            centr2 = np.copy(psf_centr2)

        # fit centroid data with low-order polynomial
        cfit = np.zeros((len(centr2)))
        csig = np.zeros((len(centr2)))
        functype = 'poly' + str(npoly_cxcy)
        pinit = array([np.nanmean(centr2)])
        if npoly_cxcy > 0:
            for j in range(npoly_cxcy):
                pinit = np.append(pinit,0.0)
        try:
            coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty, status = \
                kepfit.lsqclip(functype,pinit,centr1,centr2,None,sigma_cxcy,sigma_cxcy,10,logfile,verbose)
            for j in range(len(coeffs)):
                cfit += coeffs[j] * np.power(centr1, j)
                csig[:] = sigma
        except:
            errmsg = ('ERROR -- KEPSFF: could not fit centroid data with'
                      ' polynomial. There are no data points within the'
                      ' range of input rows {} - {}. Either increase the'
                      ' stepsize (with an appreciation of the effects on'
                      ' light curve quality this will have!), or better yet'
                      ' - cut the timeseries up to remove large gaps in the'
                      'input light curve using kepclip.'.format(t1,t2))
            kepmsg.err(logfile, errmsg, verbose)

        # reject outliers
        time_good = array([], 'float64')
        centr1_good = array([], 'float32')
        centr2_good = array([], 'float32')
        flux_good = array([], 'float32')
        cad_good = array([], 'int')
        for i in range(len(cfit)):
            if abs(centr2[i] - cfit[i]) < sigma_cxcy * csig[i]:
                time_good = np.append(time_good, intime[i])
                centr1_good = np.append(centr1_good, centr1[i])
                centr2_good = np.append(centr2_good, centr2[i])
                flux_good = np.append(flux_good, indata[i])
                cad_good = np.append(cad_good, cadenceno[i])

        # covariance matrix for centroid time series
        centr = np.concatenate([[centr1_good] - np.mean(centr1_good),
                                [centr2_good] - np.mean(centr2_good)])
        covar = np.cov(centr)
        # eigenvector eigenvalues of covariance matrix
        [eval, evec] = np.linalg.eigh(covar)
        ex = arange(-10.0, 10.0, 0.1)
        epar = evec[1, 1] / evec[0, 1] * ex
        enor = evec[1, 0] / evec[0, 0] * ex
        ex = ex + np.mean(centr1)
        epar = epar + np.mean(centr2_good)
        enor = enor + np.mean(centr2_good)
        # rotate centroid data
        centr_rot = np.dot(evec.T, centr)
        # fit polynomial to rotated centroids
        rfit = np.zeros((len(centr2)))
        rsig = np.zeros((len(centr2)))
        functype = getattr(kepfunc, 'poly' + str(npoly_ardx))
        pinit = np.array([nanmean(centr_rot[0,:])])
        pinit = np.array([1.0])
        if npoly_ardx > 0:
            for j in range(npoly_ardx):
                pinit = append(pinit,0.0)
        try:
            coeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit, centr_rot[1, :],
                               centr_rot[0, :], None, 100.0, 100.0, 1, logfile,
                               verbose)
        except:
            errmsg = ('ERROR -- KEPSFF: could not fit rotated centroid data'
                      ' with polynomial')
            kepmsg.err(logfile, errmsg, verbose)
        rx = np.linspace(np.nanmin(centr_rot[1, :]),
                         np.nanmax(centr_rot[1, :]), 100)
        ry = np.zeros((len(rx)))
        for i in range(len(coeffs)):
            ry = ry + coeffs[i] * numpy.power(rx,i)

        # calculate arclength of centroids
        s = np.zeros((len(rx)))
        for i in range(1, len(s)):
            work3 = ((ry[i] - ry[i - 1]) / (rx[i] - rx[i-1])) ** 2
            s[i] = s[i - 1] + math.sqrt(1.0 + work3) * (rx[i] - rx[i-1])

        # fit arclength as a function of strongest eigenvector
        sfit = np.zeros((len(centr2)))
        ssig = np.zeros((len(centr2)))
        functype = getattr(kepfunc, 'poly' + str(npoly_ardx))
        pinit = np.array([nanmean(s)])
        if npoly_ardx > 0:
            for j in range(npoly_ardx):
                pinit = np.append(pinit, 0.0)
        try:
            acoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plotx, ploty = \
                kepfit.lsqclip(functype, pinit, rx, s, None, 100.0, 100.0, 100,
                               logfile, verbose)
        except:
            errmsg = ('ERROR -- KEPSFF: could not fit arclength data'
                      ' with polynomial')
            kepmsg.err(logfile, message, verbose)

        # correlate arclength with detrended flux
        t = np.copy(time_good)
        c = np.copy(cad_good)
        y = np.copy(flux_good)
        z = np.centr_rot[1, :]
        x = np.zeros((len(z)))
        for i in range(len(acoeffs)):
            x = x + acoeffs[i] * np.power(z, i)

        # calculate time derivative of arclength s
        dx = zeros((len(x)))
        for i in range(1,len(x)):
            dx[i] = (x[i] - x[i-1]) / (t[i] - t[i-1])
        dx[0] = dx[1]

        # fit polynomial to derivative and flag outliers (thruster firings)
        dfit = np.zeros((len(dx)))
        dsig = np.zeros((len(dx)))
        functype = 'poly' + str(npoly_dsdt)
        pinit = np.array([np.nanmean(dx)])
        if npoly_dsdt > 0:
            for j in range(npoly_dsdt):
                pinit = np.append(pinit, 0.0)
        try:
            dcoeffs, errors, covar, iiter, dsigma, chi2, dof, fit, dumx, dumy = \
                kepfit.lsqclip(functype, pinit, t, dx, None, 3.0, 3.0, 10,
                               logfile, verbose)
        except:
            errmsg = ('ERROR -- KEPSFF: could not fit arclength derivative'
                      ' with polynomial.')
            kepmsg.err(logfile, errmsg, verbose)
        for i in range(len(dcoeffs)):
            dfit = dfit + dcoeffs[i] * np.power(t,i)
        centr1_pnt = np.array([], 'float32')
        centr2_pnt = np.array([], 'float32')
        time_pnt = np.array([], 'float64')
        flux_pnt = np.array([], 'float32')
        dx_pnt = np.array([], 'float32')
        s_pnt = np.array([], 'float32')
        time_thr = np.array([], 'float64')
        flux_thr = np.array([], 'float32')
        dx_thr = np.array([], 'float32')
        thr_cadence = []
        for i in range(len(t)):
            if (dx[i] < dfit[i] + sigma_dsdt * dsigma
                and dx[i] > dfit[i] - sigma_dsdt * dsigma):
                time_pnt = np.append(time_pnt, time_good[i])
                flux_pnt = np.append(flux_pnt, flux_good[i])
                dx_pnt = np.append(dx_pnt, dx[i])
                s_pnt = np.append(s_pnt, x[i])
                centr1_pnt = np.append(centr1_pnt, centr1_good[i])
                centr2_pnt = np.append(centr2_pnt, centr2_good[i])
            else:
                time_thr = np.append(time_thr, time_good[i])
                flux_thr = np.append(flux_thr, flux_good[i])
                dx_thr = append(dx_thr, dx[i])
                thr_cadence.append(cad_good[i])

        # fit arclength-flux correlation
        cfit = np.zeros((len(time_pnt)))
        csig = np.zeros((len(time_pnt)))
        functype = getattr(kepfunc, 'poly' + str(npoly_arfl))
        pinit = np.array([np.nanmean(flux_pnt)])
        if npoly_arfl > 0:
            for j in range(npoly_arfl):
                pinit = np.append(pinit, 0.0)
        try:
            ccoeffs, errors, covar, iiter, sigma, chi2, dof, fit, plx, ply = \
                kepfit.lsqclip(functype, pinit, s_pnt, flux_pnt, None,
                               sigma_arfl, sigma_arfl, 100, logfile, verbose)
        except:
            errmsg = ('ERROR -- KEPSFF: could not fit arclength-flux'
                      ' correlation with polynomial')
            kepmsg.err(logfile, message, verbose)

        # correction factors for unfiltered data
        centr = np.concatenate([[centr1] - np.mean(centr1_good),
                                [centr2] - np.mean(centr2_good)])
        centr_rot = np.dot(evec.T, centr)
        yy = np.copy(indata)
        zz = centr_rot[1, :]
        xx = np.zeros((len(zz)))
        cfac = np.zeros((len(zz)))
        for i in range(len(acoeffs)):
            xx = xx + acoeffs[i] * np.power(zz,i)
        for i in range(len(ccoeffs)):
            cfac = cfac + ccoeffs[i] * np.power(xx,i)

        # apply correction to flux time-series
        out_detsap = indata / cfac
        # split time-series data for plotting
        tim_gd = array([], 'float32')
        flx_gd = array([], 'float32')
        tim_bd = array([], 'float32')
        flx_bd = array([], 'float32')
        for i in range(len(indata)):
            if intime[i] in time_pnt:
                tim_gd = np.append(tim_gd, intime[i])
                flx_gd = np.append(flx_gd, out_detsap[i])
            else:
                tim_bd = np.append(tim_bd, intime[i])
                flx_bd = np.append(flx_bd, out_detsap[i])
        # plot style and size
        kepplot.define(16, 14, logfile, verbose)
        plt.figure(figsize=[20, 8])
        plt.clf()

        # plot x-centroid vs y-centroid
        ax = kepplot.location([0.04, 0.57, 0.16, 0.41])
        px = np.copy(centr1)
        py = np.copy(centr2)
        pxmin = px.min()
        pxmax = px.max()
        pymin = py.min()
        pymax = py.max()
        pxr = pxmax - pxmin
        pyr = pymax - pymin
        pad = 0.05
        if pxr > pyr:
            dely = (pxr - pyr) / 2
            plt.xlim(pxmin - pxr * pad, pxmax + pxr * pad)
            plt.ylim(pymin - dely - pyr * pad, pymax + dely + pyr * pad)
        else:
            delx = (pyr - pxr) / 2
            plt.ylim(pymin - pyr * pad, pymax + pyr * pad)
            plt.xlim(pxmin - delx - pxr * pad, pxmax + delx + pxr * pad)
        pylab.plot(px, py, color='#980000', markersize=5, marker='D', ls='')
        pylab.plot(centr1_good,centr2_good,color='#009900',markersize=5, marker='D', ls='')
        pylab.plot(ex,epar,color='k',ls='-')
        pylab.plot(ex,enor,color='k',ls='-')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        kepplot.labels('CCD Column', 'CCD Row', 'k', 16)
        plt.grid()

        # plot arclength fits vs drift along strongest eigenvector
        ax = kepplot.location([0.24, 0.57, 0.16, 0.41])
        px = rx - rx[0]
        py = s - rx - (s[0] - rx[0])
        py, ylab, status = kepplot.cleany(py, 1.0, logfile, verbose)
        kepplot.RangeOfPlot(px, py, 0.05, False)
        pylab.plot(px,py,color='#009900',markersize=5,marker='D',ls='')
        px = plotx - rx[0]
        py = ploty - plotx - (s[0] - rx[0])
        py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
        pylab.plot(px, py, color='r', ls='-', lw=3)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        ylab = re.sub(' e\S+', ' pixels)', ylab)
        ylab = re.sub(' s\S+', '', ylab)
        ylab = re.sub('Flux', 's $-$ x\'', ylab)
        kepplot.labels('Linear Drift [x\'] (pixels)', ylab, 'k', 16)
        plt.grid()

        # plot time derivative of arclength s
        ax = kepplot.location([0.04,0.08,0.16,0.41])
        px = copy(time_pnt)
        py = copy(dx_pnt)
        px, xlab, status = kepplot.cleanx(px,logfile,verbose)
        kepplot.RangeOfPlot(px, dx, 0.05, False)
        pylab.plot(px, py, color='#009900', markersize=5, marker='D', ls='')
        try:
            px = copy(time_thr)
            py = copy(dx_thr)
            px, xlab, status = kepplot.cleanx(px, logfile, verbose)
            pylab.plot(px, py, color='#980000', markersize=5, marker='D', ls='')
        except:
            pass
        px = np.copy(t)
        py = np.copy(dfit)
        px, xlab = kepplot.cleanx(px, logfile, verbose)
        plt.plot(px, py, color='r', ls='-', lw=3)
        py = np.copy(dfit + sigma_dsdt * dsigma)
        pylab.plot(px, py, color='r', ls='--', lw=3)
        py = copy(dfit-sigma_dsdt*dsigma)
        pylab.plot(px,py,color='r',ls='--',lw=3)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        kepplot.labels(xlab, 'ds/dt (pixels day$^{-1}$)', 'k', 16)
        plt.grid()
        # plot relation of arclength vs detrended flux
        ax = kepplot.location([0.24, 0.08, 0.16, 0.41])
        px = np.copy(s_pnt)
        py = np.copy(flux_pnt)
        py, ylab = kepplot.cleany(py, 1.0, logfile, verbose)
        kepplot.RangeOfPlot(px, py, 0.05, False)
        pylab.plot(px, py, color='#009900', markersize=5, marker='D', ls='')
        pylab.plot(plx, ply, color='r', ls='-', lw=3)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        kepplot.labels('Arclength [s] (pixels)', ylab, 'k', 16)
        plt.grid()

        # plot aperture photometry
        kepplot.location([0.44, 0.53, 0.55, 0.45])
        px, xlab, status = kepplot.cleanx(intime, logfile, verbose)
        py, ylab, status = kepplot.cleany(indata, 1.0, logfile, verbose)
        kepplot.RangeOfPlot(px, py, 0.01, True)
        kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2, True)
        kepplot.labels(' ',ylab,'k',16)
        pylab.setp(pylab.gca(),xticklabels=[])
        kepplot.labels(xlab,re.sub('Flux','Aperture Flux',ylab),'k',16)
        pylab.grid()
        kepplot.location([0.44, 0.08, 0.55, 0.45])
        kepplot.RangeOfPlot(px, py, 0.01, True)
        px, xlab = kepplot.cleanx(tim_gd, logfile, verbose)
        py, ylab = kepplot.cleany(flx_gd, 1.0, logfile, verbose)
        kepplot.plot1d(px, py, cadence, '#0000ff', 1.0, '#ffff00', 0.2, True)
        try:
            px, xlab = kepplot.cleanx(tim_bd,logfile,verbose)
            py = np.copy(flx_bd)
            plt.plot(px, py,color='#980000',markersize=5,marker='D',ls='')
        except:
            pass
        kepplot.labels(xlab,re.sub('Flux','Corrected Flux',ylab),'k',16)   # labels
        pylab.grid()                                                    # grid lines
        # render plot
        if plotres:
            kepplot.render(cmdLine)
        # save plot to file
        if plotres:
            pylab.savefig(re.sub('.fits','_%d.png' % (iw + 1),outfile))

        # correct fluxes within the output file
        intime = work1[:,7] + bjdref
        cadenceno = work1[:,6].astype(int)
        indata = work1[:,5]
        mom_centr1 = work1[:,4]
        mom_centr2 = work1[:,3]
        psf_centr1 = work1[:,2]
        psf_centr2 = work1[:,1]
        centr1 = copy(mom_centr1)
        centr2 = copy(mom_centr2)
        centr = concatenate([[centr1] - mean(centr1_good), [centr2] - mean(centr2_good)])
        centr_rot = dot(evec.T,centr)
        yy = copy(indata)
        zz = centr_rot[1,:]
        xx = zeros((len(zz)))
        cfac = zeros((len(zz)))
        for i in range(len(acoeffs)):
            xx = xx + acoeffs[i] * numpy.power(zz,i)
        for i in range(len(ccoeffs)):
            cfac = cfac + ccoeffs[i] * numpy.power(xx,i)
        out_detsap = yy / cfac
        instr[1].data.field('SAP_FLUX')[t1:t2] /= cfac
        instr[1].data.field('PDCSAP_FLUX')[t1:t2] /= cfac
        try:
            instr[1].data.field('DETSAP_FLUX')[t1:t2] /= cfac
        except:
            pass

# add quality flag to output file for thruster firings

        for i in range(len(intime)):
            if cadenceno[i] in thr_cadence:
                instr[1].data.field('SAP_QUALITY')[t1+i] += 131072

# write output file

    if status == 0:
        instr.writeto(outfile)
    
# close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

# end time

    if (status == 0):
	    message = 'KEPSFF completed at'
    else:
	    message = '\nKEPSFF aborted at'
    kepmsg.clock(message,logfile,verbose)


# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Correct aperture photmetry using target motion')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input FITS file', type=str)
    parser.add_argument('outfile', help='Name of output FITS file', type=str)
    parser.add_argument('--datacol', default='DETSAP_FLUX', help='Name of data column', type=str)
    parser.add_argument('--cenmethod', default='moments', help='Use which centroiding method, center0f-light or PSF fit?', type=str, choices=['moments','psf'])
    parser.add_argument('--stepsize', default=4.0, help='Stepsize over which to calibrate data [days]', type=float)
    parser.add_argument('--npoly_cxcy', default=1, help='Order of ploynomial fit to target centroids', type=int)
    parser.add_argument('--sigma_cxcy', default=6.0, help='Sigma-clipping threshold for fit to target centroids [sigma]', type=float)
    parser.add_argument('--npoly_ardx', default=6, help='Order of ploynomial fit for thruster firing detection', type=int)
    parser.add_argument('--npoly_dsdt', default=2, help='Order of ploynomial fit for thruster firing detection', type=int)
    parser.add_argument('--sigma_dsdt', default=3.0, help='Sigma-clipping threshold for thruster firing detection [sigma]', type=float)
    parser.add_argument('--npoly_arfl', default=3, help='Order of ploynomial for for arclength-flux calibration', type=int)
    parser.add_argument('--sigma_arfl', default=3.0, help='Sigma-clipping threshold for arclength-flux calibration [sigma]', type=float)
    parser.add_argument('--plotres', action='store_true', help='Save hardcopies of the plots?')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepsff.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    kepsff(args.infile,args.outfile,args.datacol,args.cenmethod,args.stepsize,
           args.npoly_cxcy,args.sigma_cxcy,args.npoly_ardx,args.npoly_dsdt,
           args.sigma_dsdt,args.npoly_arfl,args.sigma_arfl,args.plotres,
           args.clobber,args.verbose,args.logfile,args.status,cmdLine)    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepsff.par")
    t = iraf.IrafTaskFactory(taskname="kepsff", value=parfile, function=kepsff)
