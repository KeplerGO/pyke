import lightcurve as tmod
import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import fits as pyfits
from scipy.optimize import leastsq, fmin
#remove this line
import sys
#sys.path.append('/Users/tom/svn_code/PyKE/kepler/')
import kepio, kepmsg, kepkey, kepfit, kepstat


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
    plt.xlabel('Time (BJD)', {'color' : 'k'})
    plt.ylabel('Flux', {'color' : 'k'})

    #fold data
    phi, fluxfold, modelfold, errorfold, notused = fold_data(time,model,flux,
        error,period,T0)

    ax2 = plt.subplot(212)

    plt.plot(phi,fluxfold,color='#0000ff',linestyle='-',linewidth=1.0)
    plt.plot(phi,modelfold,color='red',linestyle='-',
        linewidth=2.0)
    time2 = np.insert(phi,[0],[phi[0]]) 
    time2 = np.append(time2,[phi[-1]])
    flux2 = np.insert(modelfold,[0],[0.0]) 
    flux2 = np.append(flux2,[0.0])
    plt.fill(time2,flux2,fc='#FFFACD',linewidth=0.0)
    plt.xlim([-10.,10.])
    plt.ylim([ymin-yr,ymax+yr])
    plt.xlabel('Hours from mid-transit', {'color' : 'k'})
    plt.ylabel('Flux', {'color' : 'k'})
    
    if cmdLine: 
        plt.show()
    else: 
        plt.ion()
        plt.plot([])
        plt.ioff()


def keptransitmodel(inputfile,datacol,errorcol,period_d,rprs,T0,
    Ecc,ars,inc,omega,LDparams,sec,norm=False,
    verbose=0,logfile='logfile.dat',status=0,cmdLine=False):
    

    #write to a logfile
    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPTRANSIT -- '
    call += 'inputfile='+inputfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'errorcol='+str(errorcol)+' '
    call += 'period_d='+str(period_d)+' '
    call += 'rprs='+str(rprs)+' '
    call += 'T0='+str(T0)+' '
    call += 'Ecc='+str(Ecc)+' '
    call += 'ars='+str(ars)+' '
    call += 'inc='+str(inc)+' '
    call += 'omega='+str(omega)+' '
    call += 'LDparams='+str(LDparams)+' '
    call += 'sec='+str(sec)+' '
    #to finish


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

# filter input data table

    if status == 0:
        try:
            nanclean = instr[1].header['NANCLEAN']
        except:
            naxis2 = 0
            try:
                for i in range(len(table.field(0))):
                    if np.isfinite(table.field('barytime')[i]) and \
                            np.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
            except:
                for i in range(len(table.field(0))):
                    if np.isfinite(table.field('time')[i]) and \
                            np.isfinite(table.field(datacol)[i]):
                        table[naxis2] = table[i]
                        naxis2 += 1
                        instr[1].data = table[:naxis2]
#            comment = 'NaN cadences removed from data'
#            status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)
 
# read table columns

    if status == 0:
        try:
            intime = instr[1].data.field('barytime') + 2.4e6
        except:
            intime, status = kepio.readfitscol(inputfile,instr[1].data,'time',logfile,verbose)
        
        indata, status = kepio.readfitscol(inputfile,instr[1].data,datacol,logfile,verbose)
        inerr, status = kepio.readfitscol(inputfile,instr[1].data,errorcol,logfile,verbose)
    if status == 0:
        intime = intime + bjdref
        indata = indata / cadenom
        inerr = inerr / cadenom

    if status == 0 and norm:
        #first remove outliers before normalizing
        threesig = 3.* np.std(indata)
        mask = np.logical_and(indata< indata + threesig,indata > indata - threesig)
        #now normalize
        indata = indata / np.median(indata[mask])

    if status == 0:
        #need to check if LD params are sensible and in right format
        LDparams = [float(i) for i in LDparams.split()]

        inc = inc * np.pi / 180.


    if status == 0:
        modelfit = tmod.lightcurve(intime,period_d,rprs,T0,Ecc,
            ars,inc,omega,LDparams,sec)

    if status == 0:
        phi, fluxfold, modelfold, errorfold, phiNotFold = fold_data(intime, 
            modelfit,indata,inerr,period_d,T0)

    if status == 0:
        do_plot(intime,modelfit,indata,inerr,period_d,T0,cmdLine)


 

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
    parser.add_argument('--norm', action='store_true', help='Normalize data to unity')
#    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
#    parser.add_argument('--plot', '-p', action='store_true', help='Plot result?', dest='plot')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='keptransit.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    cmdLine=True


    


    args = parser.parse_args()
    
    keptransitmodel(args.inputfile, args.datacol, args.errorcol,
        args.periodini_d, args.rprsini, args.T0ini,
        args.Eccini, args.arsini, args.incini, args.omegaini, args.LDparams, args.secini,
        args.norm,
        args.verbose, args.logfile, args.status,
        cmdLine)
    
    

else:
    from pyraf import iraf
    
    parfile = iraf.osfn("kepler$keptransitmodel.par")
    t = iraf.IrafTaskFactory(taskname="keptransitmodel", value=parfile, function=keptransitmodel)
    







