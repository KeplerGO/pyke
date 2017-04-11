"""
Name: keptransit.py
Written by: Tom Barclay
Date released: ??

Changelog:
1.0 released

__svnid__ = "$Id$"
__url__ = "$URL$"

This code is intended to fit a transit model to a Kepler light curve. We assume that the data has been cleaned
in some way to remove instrumental signals and stellar variability.

Reference:
The transit model a Mandel and Agol model
<http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>.

This code calls a a module called lightcurve
This code was created by Ian Crossfield <http://www.astro.ucla.edu/~ianc/>
and Susanne Aigraine.

The lighcurve module has been modified by TSB in order to sample the model
on a finer grid than the original data.

"""

import lightcurve as tmod
import matplotlib.pyplot as plt
import numpy as np
from numpy import logical_and, isfinite, nan
import pyfits
from scipy.optimize import leastsq, fmin
#remove this line
import sys
#sys.path.append('/Users/tom/svn_code/PyKE/kepler/')
import kepio, kepmsg, kepkey, kepfit, kepstat

np.seterr(all="ignore")

def cutBadData(date,flux,err,removeflaggeddata,qualflag):
    """
    this function finds cadences with bad data in and removes them returning only cadences which contain good data
    """
    mask = logical_and(logical_and(isfinite(date),
    isfinite(flux)),flux != 0.0)
    if removeflaggeddata:
        quality = np.where(qualflag == 0,True,False)
        finmask = logical_and(mask,quality)
    else:
        finmask = mask


    date2 = date[finmask]
    flux2 = flux[finmask]
    err2 = err[finmask]
    bad_data = finmask
    return date2,flux2,err2,bad_data

def get_chi2(obs,model,error):
    chi2 = (obs - model)**2 / error**2
    return chi2

def fit_tmod(params,args):
    period_d,rprs,T0, Ecc,ars, inc, omega, sec, fluxoffset = params
    LDparams,time,flux,error,fixed_dict,guess_params = args

    #fix things and stop params going crazy
    #here!
    if fixed_dict['period'] == True:
        period_d = guess_params[0]
    if fixed_dict['rprs'] == True:
        rprs = guess_params[1]
    if fixed_dict['T0'] == True:
        T0 = guess_params[2]
    if fixed_dict['Ecc'] == True:
        Ecc = guess_params[3]
    if fixed_dict['ars'] == True:
        ars = guess_params[4]
    if fixed_dict['inc'] == True:
        inc = guess_params[5]
    if fixed_dict['omega'] == True:
        inc = guess_params[6]
    if fixed_dict['sec'] == True:
        inc = guess_params[7]

    if fixed_dict['fluxoffset'] == False:
        flux = flux + fluxoffset

    if inc > np.pi / 2.:
        return 10**(12.)

    if omega > np.pi * 2.:
        return 10**(12.)

    mod_output = tmod.lightcurve(time,period_d,rprs,T0,Ecc,ars, inc, omega, LDparams, sec)
    model = mod_output
    chi2 = get_chi2(flux,model,error)

    return chi2

def fit_tmod2(params,LDparams,time,flux,error,fixed_dict,guess_params):
    period_d,rprs,T0, Ecc,ars, inc, omega, sec, fluxoffset = params

    #fix things and stop params going crazy
    #here!
    if fixed_dict['period'] == True:
        period_d = guess_params[0]
    if fixed_dict['rprs'] == True:
        rprs = guess_params[1]
    if fixed_dict['T0'] == True:
        T0 = guess_params[2]
    if fixed_dict['Ecc'] == True:
        Ecc = guess_params[3]
    if fixed_dict['ars'] == True:
        ars = guess_params[4]
    if fixed_dict['inc'] == True:
        inc = guess_params[5]
    if fixed_dict['omega'] == True:
        omega = guess_params[6]
    if fixed_dict['sec'] == True:
        sec = guess_params[7]

    if fixed_dict['fluxoffset'] == False:
        flux = flux + fluxoffset

    if inc > np.pi / 2.:
        return 10**(12.)

    if omega > np.pi * 2.:
        return 10**(12.)

    mod_output = tmod.lightcurve(time,period_d,rprs,T0,Ecc,ars, inc, omega, LDparams, sec)
    model = mod_output

    #try ignoring unity elements in the model
    #mask = model != 1.0

    #flux2 = flux[mask]
    #model2 = model[mask]
    #error2 = error[mask]


    chi2 = get_chi2(flux,model,error)

    return np.sum(chi2)

def fix_params(fixperiod,fixrprs,fixT0,
    fixEcc,fixars,fixinc,fixomega,fixsec, fixfluxoffset):
    fixed_dict = {'period' : fixperiod,
                    'rprs' : fixrprs,
                    'T0' : fixT0,
                    'Ecc' : fixEcc,
                    'ars' : fixars,
                    'inc' : fixinc,
                    'omega' : fixomega,
                    'sec' : fixsec,
                    'fluxoffset' : fixfluxoffset}
    return fixed_dict

def do_plot(time,model,flux,error,period,T0,cmdLine=False):
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

    plt.subplots_adjust(0.09,0.1,0.98,0.95,0.19,0.27)
    ax1 = plt.subplot(211)

    plt.plot(time,flux,color='#0000ff',linestyle='-',linewidth=1.0)
    plt.plot(time,model,color='red',linestyle='-',
        linewidth=2.0)
    time2 = np.insert(time,[0],[time[0]])
    time2 = np.append(time2,[time[-1]])
    flux2 = np.insert(flux,[0],[0.0])
    flux2 = np.append(flux2,[0.0])
    plt.fill(time2,flux2,fc='#FFFACD',linewidth=0.0)

    plt.xlim([min(time),max(time)])
    ymin = min(min(flux),min(model))
    ymax = max(max(flux),max(model))
    yr = (ymax - ymin) * 0.1
    plt.ylim([ymin-yr,ymax+yr])
    plt.xlabel('Time (BJD - 2454833)', {'color' : 'k'})
    plt.ylabel('Flux', {'color' : 'k'})
    plt.grid()

    #fold data
    phi, fluxfold, modelfold, errorfold, notused = fold_data(time,model,flux,
        error,period,T0)

    ax2 = plt.subplot(212)

    plt.plot(phi,fluxfold,color='#0000ff',linestyle='-',linewidth=1.0)
    plt.plot(phi,modelfold,color='red',linestyle='-',
        linewidth=2.0)
    time2 = np.insert(phi,[0],[phi[0]])
    time2 = np.append(time2,[phi[-1]])
    flux2 = np.insert(fluxfold,[0],[0.0])
    flux2 = np.append(flux2,[0.0])
    plt.fill(time2,flux2,fc='#FFFACD',linewidth=0.0)
    plt.xlim([-10.,10.])
    plt.ylim([ymin-yr,ymax+yr])
    plt.xlabel('Hours from mid-transit', {'color' : 'k'})
    plt.ylabel('Flux', {'color' : 'k'})
    plt.grid()
    if cmdLine:
        plt.show()
    else:
        plt.ion()
        plt.plot([])
        plt.ioff()

def fold_data(time,model,flux,error,period,T0):
    date1 = (time - T0) + period
    date1 = (time - T0) + 0.5*period
    phi1 = (((date1 / period) -
        np.floor(date1/period)) *24.*period) - 12.*period

    sort_mask = np.argsort(phi1)

    phi = phi1[sort_mask]
    fluxfold = flux[sort_mask]
    modelfold = model[sort_mask]
    errorfold = error[sort_mask]
    return phi, fluxfold, modelfold, errorfold, phi1


def make_outfile(fitsfile,outfile,phi,model,baddata):
    """
    creates a fits file identical to the input fits file save from containing two extra columns - TRANSIT_MODL and PHASE which are the sum of basis vectors fit to the data and the resulting corrected flux after the basis vector fit has been subtracted
    """
    newmodel = putInNans(baddata,model)
    newphi = putInNans(phi,model)
    col1 = pyfits.Column(name='TRANSIT_MODL',
        format='E13.7   ',unit='',array=newmodel)
    col2 = pyfits.Column(name='PHASE',format='E13.7   ',
        unit='days',array=newphi / 24.) # phi is in hours by default
    cols = fitsfile[1].columns + col1 + col2
    fitsfile[1] = pyfits.new_table(cols,header=fitsfile[1].header)
    fitsfile.writeto(outfile)

def putInNans(bad_data,flux):
    """
    Function finds the cadences where the data has been removed using cutBadData() and puts data back in. The flux data put back in is nan.
    This function is used when writing data to a FITS files.
    bad_data == True means the datapoint is good!!
    """
    newflux = np.zeros(len(bad_data))
    j = 0
    for i in range(len(bad_data)):
        if bad_data[i] == True:
            newflux[i] = flux[j]
            j += 1
        elif bad_data[i] == False:
            newflux[i] = nan
    return newflux

def keptransit(inputfile,outputfile,datacol,errorcol,periodini_d,rprsini,T0ini,
    Eccini,arsini,incini,omegaini,LDparams,secini,fixperiod,fixrprs,fixT0,
    fixEcc,fixars,fixinc,fixomega,fixsec,fixfluxoffset,removeflaggeddata,ftol=0.0001,fitter='nothing',norm=False,
    clobber=False, plot=True,verbose=0,logfile='logfile.dat',status=0,cmdLine=False):
    """
    tmod.lightcurve(xdata,period,rprs,T0,Ecc,ars, incl, omega, ld, sec)

    input transit parameters are
    Period in days
    T0
    rplanet / rstar
    a / rstar
    inclination

    limb darkening code number:
    0 = uniform
    1 = linear
    2 = quadratic
    3 = square root
    4 = non linear

    LDarr:
    u      -- linear limb-darkening (set NL=1)
    a, b   -- quadratic limb-darkening (set NL=2)
    c,  d  -- root-square limb-darkening (set NL= -2)
    a1, a2, a3, a4 -- nonlinear limb-darkening  (set NL=4)
    Nothing at all -- uniform limb-darkening (set NL=0)
    """

    np.seterr(all="ignore")

    #write to a logfile
    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPTRANSIT -- '
    call += 'inputfile='+inputfile+' '
    call += 'outputfile='+outputfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'errorcol='+str(errorcol)+' '
    call += 'periodini_d='+str(periodini_d)+' '
    call += 'rprsini='+str(rprsini)+' '
    call += 'T0ini='+str(T0ini)+' '
    call += 'Eccini='+str(Eccini)+' '
    call += 'arsini='+str(arsini)+' '
    call += 'incini='+str(incini)+' '
    call += 'omegaini='+str(omegaini)+' '
    call += 'LDparams='+str(LDparams)+' '
    call += 'secini='+str(secini)+' '
    call += 'fixperiod='+str(fixperiod)+' '
    call += 'fixrprs='+str(fixrprs)+' '
    call += 'fixT0='+str(fixT0)+' '
    call += 'fixEcc='+str(fixEcc)+' '
    call += 'fixars='+str(fixars)+' '
    call += 'fixinc='+str(fixinc)+' '
    call += 'fixomega='+str(fixomega)+' '
    call += 'fixsec='+str(fixsec)+' '
    call += 'fixfluxoffset='+str(fixfluxoffset)+' '
    call += 'removeflaggeddata='+str(removeflaggeddata)+' '
    call += 'ftol='+str(ftol)+' '
    call += 'fitter='+str(fitter)+' '
    call += 'norm='+str(norm)+' '

    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    #chatter = 'n'
    #if (verbose): chatter = 'y'
    #call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)


    kepmsg.clock('KEPTRANSIT started at',logfile,verbose)

    # test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber:
        status = kepio.clobber(outputfile,logfile,verbose)
    if kepio.fileexists(outputfile):
        message = 'ERROR -- KEPTRANSIT: ' + outputfile + ' exists. Use clobber=yes'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    if status == 0:
        instr, status = kepio.openfits(inputfile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,
            inputfile,logfile,verbose,status)
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# read table structure

    if status == 0:
        table, status = kepio.readfitstab(inputfile,instr[1],logfile,verbose)

    if status == 0:
        intime_o = table.field('time')
        influx_o = table.field(datacol)
        inerr_o = table.field(errorcol)
        try:
            qualflag = table.field('SAP_QUALITY')
        except:
            qualflag = np.zeros(len(intime_o))

    if status == 0:
        intime, indata, inerr, baddata = cutBadData(intime_o, influx_o, inerr_o,removeflaggeddata,qualflag)

    if status == 0 and norm:
        #first remove outliers before normalizing
        threesig = 3.* np.std(indata)
        mask = np.logical_and(indata< indata + threesig,indata > indata - threesig)
        #now normalize
        indata = indata / np.median(indata[mask])

    if status == 0:
        #need to check if LD params are sensible and in right format
        LDparams = [float(i) for i in LDparams.split()]

        incini = incini * np.pi / 180.

        omegaini = omegaini * np.pi / 180.

    if arsini*np.cos(incini) > 1.0 + rprsini:
        message = 'The guess inclination and a/r* values result in a non-transing planet'
        status = kepmsg.err(logfile,message,verbose)

    if status == 0:
        fixed_dict = fix_params(fixperiod,fixrprs,fixT0,
            fixEcc,fixars,fixinc,fixomega,fixsec,fixfluxoffset)

    #force flux offset to be guessed at zero
    fluxoffsetini = 0.0

    if status == 0:
        guess_params = [periodini_d,rprsini,T0ini,Eccini,arsini, incini, omegaini,
        secini,fluxoffsetini]

        print 'cleaning done: about to fit transit'

        if fitter == 'leastsq':
            fit_output = leastsq(fit_tmod,guess_params,
                args=[LDparams,intime,indata,inerr,fixed_dict,guess_params],
                full_output=True,ftol=ftol)
        elif fitter == 'fmin':

            fit_output = fmin(fit_tmod2,guess_params,
                args=[LDparams,intime,indata,inerr,fixed_dict,guess_params],
                full_output=True,ftol=ftol,xtol=ftol)


    if status == 0:
        if fixed_dict['period'] == True:
            newperiod = guess_params[0]
            print 'Fixed period (days) = ' + str(newperiod)
        else:
            newperiod = fit_output[0][0]
            print 'Fit period (days) = ' + str(newperiod)
        if fixed_dict['rprs'] == True:
            newrprs = guess_params[1]
            print 'Fixed R_planet / R_star = ' + str(newrprs)
        else:
            newrprs = fit_output[0][1]
            print 'Fit R_planet / R_star = ' + str(newrprs)
        if fixed_dict['T0'] == True:
            newT0 = guess_params[2]
            print 'Fixed T0 (BJD) = ' + str(newT0)
        else:
            newT0 = fit_output[0][2]
            print 'Fit T0 (BJD) = ' + str(newT0)
        if fixed_dict['Ecc'] == True:
            newEcc = guess_params[3]
            print 'Fixed eccentricity = ' + str(newEcc)
        else:
            newEcc = fit_output[0][3]
            print 'Fit eccentricity = ' + str(newEcc)
        if fixed_dict['ars'] == True:
            newars = guess_params[4]
            print 'Fixed a / R_star = ' + str(newars)
        else:
            newars = fit_output[0][4]
            print 'Fit a / R_star = ' + str(newars)
        if fixed_dict['inc'] == True:
            newinc = guess_params[5]
            print 'Fixed inclination (deg) = ' + str(newinc* 180. / np.pi)
        else:
            newinc = fit_output[0][5]
            print 'Fit inclination (deg) = ' + str(newinc* 180. / np.pi)
        if fixed_dict['omega'] == True:
            newomega = guess_params[6]
            print 'Fixed omega = ' + str(newomega)
        else:
            newomega = fit_output[0][6]
            print 'Fit omega = ' + str(newomega)
        if fixed_dict['sec'] == True:
            newsec = guess_params[7]
            print 'Fixed seconary eclipse depth = ' + str(newsec)
        else:
            newsec = fit_output[0][7]
            print 'Fit seconary eclipse depth = ' + str(newsec)
        if fixfluxoffset == False:
            newfluxoffset = fit_output[0][8]
            print 'Fit flux offset = ' + str(newfluxoffset)



        modelfit = tmod.lightcurve(intime,newperiod,newrprs,newT0,newEcc,
            newars,newinc,newomega,LDparams,newsec)

        if fixfluxoffset == False:
            modelfit += newfluxoffset


        #output to a file
        phi, fluxfold, modelfold, errorfold, phiNotFold = fold_data(intime,
            modelfit,indata,inerr,newperiod,newT0)


        make_outfile(instr,outputfile,phiNotFold,modelfit, baddata)


    # end time

    if (status == 0):
        message = 'KEPTRANSIT completed at'
    else:
        message = '\nKEPTRANSIT aborted at'
    kepmsg.clock(message,logfile,verbose)

    if plot and status == 0:
        do_plot(intime,modelfit,indata,inerr,newperiod,newT0,cmdLine)











if '--shell' in sys.argv:
    import argparse

    #inputfile,outputfile,datacol,errorcol,periodini_d,rprsini,T0ini,
    #Eccini,arsini,incini,omegaini,LDparams,secini,fixperiod,fixrprs,fixT0,
    #fixEcc,fixars,fixinc,fixomega,fixsec,ftol=0.0001,fitter='nothing'
    #clobber=False, plot=True,logfile='logfile.dat',verbose=0,status=0,cmdLine=False

    parser = argparse.ArgumentParser(description='Fit a exoplanet transit model to the light curve')
    #parser.add_argument('--infile', '-i', help='Name of input file', dest='infile',type=str,required=True)
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str, dest=inputfile)
    parser.add_argument('outfile', help='Name of output FITS file', type=str, dest=outputfile)
    parser.add_argument('--datacol', help='Column containing flux data to fit', type=str,
        default='DETSAP_FLUX')
    parser.add_argument('--errorcol', help='Column containing flux uncertainty', type=str,
        default='DETSAP_FLUX_ERR')
    parser.add_argument('period', help='Guess for planet orbital period',
        type=float, dest=periodini_d)
    parser.add_argument('rprs', help='Guess for radius of the planet / radius of the star',
        type=float, dest=rprsini)
    parser.add_argument('T0', help='Guess mid-time of first transit',
        type=float, dest=T0ini)
    parser.add_argument('--ecc', help='Guess eccentricity',
        type=float, dest=Eccini, default=0.)
    parser.add_argument('ars', help='Guess semi-major axis / radius of the star',
        type=float, dest=arsini)
    parser.add_argument('--inc', help='Guess inclination',
        type=float, dest=incini, default=90.)
    parser.add_argument('--omega', help='Guess periastron angle',
        type=float, dest=omegaini, default=0.)
    parser.add_argument('--LDparams', help='Limb darkening parameters, seperate by a space',
        type=str, default='')
    parser.add_argument('--sec', help='Guess secondary eclipse depth',
        type=float, dest=secini, default=0.)
    parser.add_argument('--fixperiod', action='store_true', help='Fix period?')
    parser.add_argument('--fixrprs', action='store_true', help='Fix rp/r*?')
    parser.add_argument('--fixT0', action='store_true', help='Fix T0?')
    parser.add_argument('--fixEcc', action='store_true', help='Fix eccentricity?')
    parser.add_argument('--fixars', action='store_true', help='Fix a/r*?')
    parser.add_argument('--fixinc', action='store_true', help='Fix inclination?')
    parser.add_argument('--fixomega', action='store_true', help='Fix period?')
    parser.add_argument('--fixsec', action='store_true', help='Fix secondary eclipse depth?')
    parser.add_argument('--fixfluxoffset', action='store_true', help='Fix the out of transit flux to unity?')
    parser.add_argument('--removeflaggeddata', action='store_true', help='Remove data with quality flag >0?')
    parser.add_argument('--ftol', help='Fix period?')
    parser.add_argument('--fitter', help='Fix period?')
    parser.add_argument('--norm', action='store_true', help='Normalize data to unity')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--plot', '-p', action='store_true', help='Plot result?', dest='plot')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='keptransit.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    cmdLine=True





    args = parser.parse_args()

    keptransit(args.inputfile, args.outputfile, args.datacol, args.errorcol,
        args.periodini_d, args.rprsini, args.T0ini,
        args.Eccini, args.arsini, args.incini, args.omegaini, args.LDparams, args.secini,
        args.fixperiod, args.fixrprs, args.fixT0,
        args.fixEcc, args.fixars, args.fixinc, args.fixomega, args.fixsec, args.fixfluxoffset,
        args.removeflaggeddata, args.ftol,
        args.fitter, args.norm,
        args.clobber, args.plot, args.verbose, args.logfile, args.status,
        cmdLine)



else:
    from pyraf import iraf

    parfile = iraf.osfn("kepler$keptransit.par")
    t = iraf.IrafTaskFactory(taskname="keptransit", value=parfile, function=keptransit)

