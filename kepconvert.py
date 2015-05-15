
import numpy, sys, time, pyfits, re
from pyfits import *
from numpy import *
import kepio, kepmsg, kepkey

__svnid__ = "$Id: kepconvert.py 6165 2014-03-26 21:16:27Z mstill $"
__url__ = "$URL: svn+ssh://mstill@murzim.amn.nasa.gov/data-repo/trunk/data/flight/go/PyKE/kepler/kepconvert.py $"


def kepconvert(infile,outfile,conversion,columns,baddata,clobber,verbose,logfile,status): 

# startup parameters

    status = 0

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPCONVERT -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'conversion='+conversion+' '
    call += 'columns='+columns+ ' '
    writebad = 'n'
    if (baddata): writebad = 'y'
    call += 'baddata='+writebad+ ' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)
    
# start time

    kepmsg.clock('KEPCONVERT started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# data columns

    if status == 0:
        colnames = columns.strip().split(',')
        ncol = len(colnames)
        if ncol < 1:
            message = 'ERROR -- KEPCONVERT: no data columns specified'
            status = kepmsg.err(logfile,message,verbose)

# input file exists

    if status == 0 and not kepio.fileexists(infile):
        message = 'ERROR -- KEPCONVERT: input file '+infile+' does not exist'
        status = kepmsg.err(logfile,message,verbose)

# clobber output file

    if status == 0:
        if clobber: status = kepio.clobber(outfile,logfile,verbose)
        if kepio.fileexists(outfile): 
            message = 'ERROR -- KEPCONVERT: ' + outfile + ' exists. Use clobber=yes'
            status = kepmsg.err(logfile,message,verbose)

# open FITS input file


    if status == 0 and conversion == 'fits2asc':
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# read FITS table data

    if status == 0 and conversion == 'fits2asc':
        table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# check columns exist in FITS file
    if not baddata and status == 0 and conversion == 'fits2asc':
	    try:
	    	qualcol = table.field('SAP_QUALITY') == 0
	    except:
	    	message = 'No SAP_QUALITY column in data, are you using an old FITS file?'
	    	status = kepmsg.err(logfile,message,verbose)

    if status == 0 and conversion == 'fits2asc':
        work = []
        for colname in colnames:
            try:
                if colname.lower() == 'time':
                    work.append(table.field(colname) + bjdref)
                else:
                    work.append(table.field(colname))
            except:
                message = 'ERROR -- KEPCONVERT: no column ' + colname + ' in ' + infile
                status = kepmsg.err(logfile,message,verbose)
        if not baddata:
            for i in range(len(work)):
                work[i] = work[i][qualcol]
# close input file

    if status == 0 and conversion == 'fits2asc':
        status = kepio.closefits(instr,logfile,verbose)

## write output file

    if status == 0 and conversion == 'fits2asc':
        # table, status = kepio.openascii(outfile,'w',logfile,verbose)
        # for i in range(len(work[0])):
            # txt = ''
            # for j in range(len(work)):
                # if numpy.isfinite(work[j][i]):
                    # txt += str(work[j][i]) + ' '
            # txt = txt.strip()
            # if len(re.sub('\s+',',',txt).split(',')) == ncol:
                # table.write(txt + '\n')
        # status = kepio.closeascii(table,logfile,verbose)
        savetxt(outfile,array(work).T)
    	
## open and read ASCII input file

    if status == 0 and conversion == 'asc2fits':
        table, status = kepio.openascii(infile,'r',logfile,verbose)

## organize ASCII table into arrays

    if status == 0 and conversion == 'asc2fits':
        work = []
        for i in range(ncol):
            work.append([])
        nline = 0
        for line in table:
            line = line.strip()
            line = re.sub('\s+',',',line)
            line = re.sub('\|',',',line)
            line = re.sub(';',',',line)
            if '#' not in line:
                nline + 1
                line = line.split(',')
                if len(line) == ncol:
                    for i in range(len(line)):
                        try:
                            work[i].append(float(line[i]))
                        except:
                            message = 'ERROR --KEPCONVERT: ' + str(line[i]) + ' is not float'
                            status = kepmsg.err(logfile,message,verbose)
                            break
                else:
                    message  = 'ERROR --KEPCONVERT: ' + str(ncol) + ' columns required but '
                    message += str(len(line)) + ' columns supplied by ' + infile
                    message += ' at line' + str(nline)
                    status = kepmsg.err(logfile,message,verbose)
                    break
        for i in range(ncol):
            work[i] = numpy.array(work[i],dtype='float64')

## timing keywords for output file

    if status == 0 and conversion == 'asc2fits':
        for i in range(ncol):
            if 'time' in colnames[i].lower():
                if work[i][1] > 54000.0 and work[i][1] < 60000.0:
                    work[i] += 2.4e6
#                work[i] += 2.4553e6
                tstart = work[i].min()
                tstop = work[i].max()
                lc_start = tstart
                lc_end = tstop
                if lc_start > 2.4e6: lc_start -= 2.4e6
                if lc_end > 2.4e6: lc_end -= 2.4e6
                dts = []
                for j in range(1,len(work[i])):
                   dts.append(work[i][j] - work[i][j-1])
                dts = numpy.array(dts,dtype='float32')
                cadence = numpy.median(dts)
                if cadence * 86400.0 > 58.0 and cadence * 86400.0 < 61.0:
                    obsmode = 'short cadence'
                elif cadence * 86400.0 > 1600.0 and cadence * 86400.0 < 2000.0:
                    obsmode = 'long cadence'
                else:
                    obsmode = 'unknown'

## Create the outfile primary extension

    if status == 0 and conversion == 'asc2fits':
        hdu0 = PrimaryHDU()
        try:
            hdu0.header.update('EXTNAME','PRIMARY','name of extension')
            hdu0.header.update('EXTVER',1.0,'extension version number')
            hdu0.header.update('ORIGIN','NASA/Ames','organization that generated this file')
            hdu0.header.update('DATE',time.asctime(time.localtime()),'file creation date')
            hdu0.header.update('CREATOR','kepconvert','SW version used to create this file')
            hdu0.header.update('PROCVER','None','processing script version')
            hdu0.header.update('FILEVER','2.0','file format version')
            hdu0.header.update('TIMVERSN','OGIP/93-003','OGIP memo number for file format')
            hdu0.header.update('TELESCOP','Kepler','telescope')
            hdu0.header.update('INSTRUME','Kepler photometer','detector type')
            hdu0.header.update('OBJECT','Unknown','string version of kepID')
            hdu0.header.update('KEPLERID','Unknown','unique Kepler target identifier')
            hdu0.header.update('CHANNEL','Unknown','CCD channel')
            hdu0.header.update('SKYGROUP','Unknown','roll-independent location of channel')
            hdu0.header.update('MODULE','Unknown','CCD module')
            hdu0.header.update('OUTPUT','Unknown','CCD output')
            hdu0.header.update('QUARTER','Unknown','mission quarter during which data was collected')
            hdu0.header.update('SEASON','Unknown','mission season during which data was collected')
            hdu0.header.update('DATA_REL','Unknown','version of data release notes describing data')
            hdu0.header.update('OBSMODE',obsmode,'observing mode')
            hdu0.header.update('RADESYS','Unknown','reference frame of celestial coordinates')
            hdu0.header.update('RA_OBJ','Unknown','[deg] right ascension from KIC')
            hdu0.header.update('DEC_OBJ','Unknown','[deg] declination from KIC')
            hdu0.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
            hdu0.header.update('PMRA','Unknown','[arcsec/yr] RA proper motion')
            hdu0.header.update('PMDEC','Unknown','[arcsec/yr] Dec proper motion')
            hdu0.header.update('PMTOTAL','Unknown','[arcsec/yr] total proper motion')
            hdu0.header.update('PARALLAX','Unknown','[arcsec] parallax')
            hdu0.header.update('GLON','Unknown','[deg] galactic longitude')
            hdu0.header.update('GLAT','Unknown','[deg] galactic latitude')
            hdu0.header.update('GMAG','Unknown','[mag] SDSS g band magnitude from KIC')
            hdu0.header.update('RMAG','Unknown','[mag] SDSS r band magnitude from KIC')
            hdu0.header.update('IMAG','Unknown','[mag] SDSS i band magnitude from KIC')
            hdu0.header.update('ZMAG','Unknown','[mag] SDSS z band magnitude from KIC')
            hdu0.header.update('D51MAG','Unknown','[mag] D51 magnitude, from KIC')
            hdu0.header.update('JMAG','Unknown','[mag] J band magnitude from 2MASS')
            hdu0.header.update('HMAG','Unknown','[mag] H band magnitude from 2MASS')
            hdu0.header.update('KMAG','Unknown','[mag] K band magnitude from 2MASS')
            hdu0.header.update('KEPMAG','Unknown','[mag] Kepler magnitude (Kp) from KIC')
            hdu0.header.update('GRCOLOR','Unknown','[mag] (g-r) color, SDSS bands')
            hdu0.header.update('JKCOLOR','Unknown','[mag] (J-K) color, 2MASS bands')
            hdu0.header.update('GKCOLOR','Unknown','[mag] (g-K) color, SDSS g - 2MASS K')
            hdu0.header.update('TEFF','Unknown','[K] effective temperature from KIC')
            hdu0.header.update('LOGG','Unknown','[cm/s2] log10 surface gravity from KIC')
            hdu0.header.update('FEH','Unknown','[log10([Fe/H])] metallicity from KIC')
            hdu0.header.update('EBMINUSV','Unknown','[mag] E(B-V) redenning from KIC')
            hdu0.header.update('AV','Unknown','[mag] A_v extinction from KIC')
            hdu0.header.update('RADIUS','Unknown','[solar radii] stellar radius from KIC')
            hdu0.header.update('TMINDEX','Unknown','unique 2MASS catalog ID from KIC')
            hdu0.header.update('SCPID','Unknown','unique SCP processing ID from KIC') 
            hdulist = HDUList(hdu0)
        except:
            message = 'ERROR -- KEPCONVERT: cannot create primary extension in ' + outfile
            status = kepmsg.err(logfile,message,verbose)
            
## create the outfile HDU 1 extension

    if status == 0 and conversion == 'asc2fits':
        try:
            fitscol = []
            for i in range(ncol):
                fitscol.append(Column(name=colnames[i],format='D',array=work[i]))
            fitscols = ColDefs(fitscol)
            hdu1 = new_table(fitscols)
            hdulist.append(hdu1)
            hdu1.header.update('INHERIT',True,'inherit primary keywords')
            hdu1.header.update('EXTNAME','LIGHTCURVE','name of extension')
            hdu1.header.update('EXTVER',1,'extension version number')
            hdu1.header.update('TELESCOP','Kepler','telescope')
            hdu1.header.update('INSTRUME','Kepler photometer','detector type')
            hdu1.header.update('OBJECT','Unknown','string version of kepID')
            hdu1.header.update('KEPLERID','Unknown','unique Kepler target identifier')
            hdu1.header.update('RADESYS','Unknown','reference frame of celestial coordinates')
            hdu1.header.update('RA_OBJ','Unknown','[deg] right ascension from KIC')
            hdu1.header.update('DEC_OBJ','Unknown','[deg] declination from KIC')
            hdu1.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
            hdu1.header.update('TIMEREF','Unknown','barycentric correction applied to times')
            hdu1.header.update('TASSIGN','Unknown','where time is assigned')
            hdu1.header.update('TIMESYS','Unknown','time system is barycentric JD')
            hdu1.header.update('BJDREFI',0.0,'integer part of BJD reference date')
            hdu1.header.update('BJDREFF',0.0,'fraction of day in BJD reference date')
            hdu1.header.update('TIMEUNIT','Unknown','time unit for TIME, TSTART and TSTOP')
            hdu1.header.update('TSTART',tstart,'observation start time in JD - BJDREF')
            hdu1.header.update('TSTOP',tstop,'observation stop time in JD - BJDREF')
            hdu1.header.update('LC_START',lc_start,'observation start time in MJD')
            hdu1.header.update('LC_END',lc_end,'observation stop time in MJD')
            hdu1.header.update('TELAPSE',tstop-tstart,'[d] TSTOP - TSTART')
            hdu1.header.update('LIVETIME','Unknown','[d] TELAPSE multiplied by DEADC')
            hdu1.header.update('EXPOSURE','Unknown','[d] time on source')
            hdu1.header.update('DEADC','Unknown','deadtime correction')
            hdu1.header.update('TIMEPIXR','Unknown','bin time beginning=0 middle=0.5 end=1')
            hdu1.header.update('TIERRELA','Unknown','[d] relative time error')
            hdu1.header.update('TIERABSO','Unknown','[d] absolute time error')
            hdu1.header.update('INT_TIME','Unknown','[s] photon accumulation time per frame')
            hdu1.header.update('READTIME','Unknown','[s] readout time per frame')
            hdu1.header.update('FRAMETIM','Unknown','[s] frame time (INT_TIME + READTIME)')
            hdu1.header.update('NUM_FRM','Unknown','number of frames per time stamp')
            hdu1.header.update('TIMEDEL','Unknown','[d] time resolution of data')
            hdu1.header.update('DATE-OBS','Unknown','TSTART as UT calendar date')
            hdu1.header.update('DATE-END','Unknown','TSTOP as UT calendar date')
            hdu1.header.update('BACKAPP','Unknown','background is subtracted')
            hdu1.header.update('DEADAPP','Unknown','deadtime applied')
            hdu1.header.update('VIGNAPP','Unknown','vignetting or collimator correction applied')
            hdu1.header.update('GAIN','Unknown','channel gain [electrons/count]')
            hdu1.header.update('READNOIS','Unknown','read noise [electrons]')
            hdu1.header.update('NREADOUT','Unknown','number of reads per cadence')
            hdu1.header.update('TIMSLICE','Unknown','time-slice readout sequence section')
            hdu1.header.update('MEANBLCK','Unknown','FSW mean black level [count]')
            hdu1.header.update('PDCSAPFL','Unknown','SAP PDC processing flags (bit code)')
            hdu1.header.update('PDCDIAFL','Unknown','DIA PDC processing flags (bit code)')
            hdu1.header.update('MISPXSAP','Unknown','no of optimal aperture pixels missing from SAP')
            hdu1.header.update('MISPXDIA','Unknown','no of optimal aperture pixels missing from DIA')
            hdu1.header.update('CROWDSAP','Unknown','crowding metric evaluated over SAP opt. ap.')
            hdu1.header.update('CROWDDIA','Unknown','crowding metric evaluated over DIA aperture')
        except:
            message = 'ERROR -- KEPCONVERT: cannot create light curve extension in ' + outfile
            status = kepmsg.err(logfile,message,verbose)

## history keyword in output file

    if status == 0 and conversion == 'asc2fits':
        status = kepkey.history(call,hdu0,outfile,logfile,verbose)

## filter data table

    if status == 0 and conversion == 'asc2fits':
        instr, status = kepio.filterNaN(hdulist,colnames[min(array([1,len(colnames)-1],dtype='int'))],
                                        outfile,logfile,verbose)

## write output FITS file

    if status == 0 and conversion == 'asc2fits':
        hdulist.writeto(outfile,checksum=True)

## end time

    if (status == 0):
	    message = 'KEPCONVERT completed at'
    else:
	    message = '\nKEPCONVERT aborted at'
    kepmsg.clock(message,logfile,verbose)

## main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert Kepler FITS time series to or from a different file format')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of output file', type=str)

    parser.add_argument('conversion', help='Type of data conversion', type=str,choices=['fits2asc','asc2fits'])
    parser.add_argument('--columns', '-c', default='TIME,SAP_FLUX,SAP_FLUX_ERR', dest='columns', help='Comma-delimited list of data columns', type=str)
    
    parser.add_argument('--baddata', action='store_false', help='Output rows which have been flagged as questionable')

    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    kepconvert(args.infile, args.outfile, args.conversion, args.columns, args.baddata, args.clobber, args.verbose, args.logfile, args.status)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepconvert.par")
    t = iraf.IrafTaskFactory(taskname="kepconvert", value=parfile, function=kepconvert)
