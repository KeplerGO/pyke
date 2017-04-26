import sys, time, re
import numpy as np
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey


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
        np.savetxt(outfile,np.array(work).T)

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
            work[i] = np.array(work[i],dtype='float64')

## timing keywords for output file

    if status == 0 and conversion == 'asc2fits':
        for i in range(ncol):
            if 'time' in colnames[i].lower():
                if work[i][1] > 54000.0 and work[i][1] < 60000.0:
                    work[i] += 2.4e6
                tstart = work[i].min()
                tstop = work[i].max()
                lc_start = tstart
                lc_end = tstop
                if lc_start > 2.4e6: lc_start -= 2.4e6
                if lc_end > 2.4e6: lc_end -= 2.4e6
                dts = []
                for j in range(1,len(work[i])):
                   dts.append(work[i][j] - work[i][j-1])
                dts = np.array(dts,dtype='float32')
                cadence = np.median(dts)
                if cadence * 86400.0 > 58.0 and cadence * 86400.0 < 61.0:
                    obsmode = 'short cadence'
                elif cadence * 86400.0 > 1600.0 and cadence * 86400.0 < 2000.0:
                    obsmode = 'long cadence'
                else:
                    obsmode = 'unknown'

## Create the outfile primary extension

    if status == 0 and conversion == 'asc2fits':
        hdu0 = pyfits.PrimaryHDU()
        try:
            hdu0.header['EXTNAME'] = ('PRIMARY','name of extension')
            hdu0.header['EXTVER'] = (1.0,'extension version number')
            hdu0.header['ORIGIN'] = ('NASA/Ames','organization that generated this file')
            hdu0.header['DATE'] = (time.asctime(time.localtime()),'file creation date')
            hdu0.header['CREATOR'] = ('kepconvert','SW version used to create this file')
            hdu0.header['PROCVER'] = ('None','processing script version')
            hdu0.header['FILEVER'] = ('2.0','file format version')
            hdu0.header['TIMVERSN'] = ('OGIP/93-003','OGIP memo number for file format')
            hdu0.header['TELESCOP'] = ('Kepler','telescope')
            hdu0.header['INSTRUME'] = ('Kepler photometer','detector type')
            hdu0.header['OBJECT'] = ('Unknown','string version of kepID')
            hdu0.header['KEPLERID'] = ('Unknown','unique Kepler target identifier')
            hdu0.header['CHANNEL'] = ('Unknown','CCD channel')
            hdu0.header['SKYGROUP'] = ('Unknown','roll-independent location of channel')
            hdu0.header['MODULE'] = ('Unknown','CCD module')
            hdu0.header['OUTPUT'] = ('Unknown','CCD output')
            hdu0.header['QUARTER'] = ('Unknown','mission quarter during which data was collected')
            hdu0.header['SEASON'] = ('Unknown','mission season during which data was collected')
            hdu0.header['DATA_REL'] = ('Unknown','version of data release notes describing data')
            hdu0.header['OBSMODE'] = (obsmode,'observing mode')
            hdu0.header['RADESYS'] = ('Unknown','reference frame of celestial coordinates')
            hdu0.header['RA_OBJ'] = ('Unknown','[deg] right ascension from KIC')
            hdu0.header['DEC_OBJ'] = ('Unknown','[deg] declination from KIC')
            hdu0.header['EQUINOX'] = (2000.0,'equinox of celestial coordinate system')
            hdu0.header['PMRA'] = ('Unknown','[arcsec/yr] RA proper motion')
            hdu0.header['PMDEC'] = ('Unknown','[arcsec/yr] Dec proper motion')
            hdu0.header['PMTOTAL'] = ('Unknown','[arcsec/yr] total proper motion')
            hdu0.header['PARALLAX'] = ('Unknown','[arcsec] parallax')
            hdu0.header['GLON'] = ('Unknown','[deg] galactic longitude')
            hdu0.header['GLAT'] = ('Unknown','[deg] galactic latitude')
            hdu0.header['GMAG'] = ('Unknown','[mag] SDSS g band magnitude from KIC')
            hdu0.header['RMAG'] = ('Unknown','[mag] SDSS r band magnitude from KIC')
            hdu0.header['IMAG'] = ('Unknown','[mag] SDSS i band magnitude from KIC')
            hdu0.header['ZMAG'] = ('Unknown','[mag] SDSS z band magnitude from KIC')
            hdu0.header['D51MAG'] = ('Unknown','[mag] D51 magnitude], from KIC')
            hdu0.header['JMAG'] = ('Unknown','[mag] J band magnitude from 2MASS')
            hdu0.header['HMAG'] = ('Unknown','[mag] H band magnitude from 2MASS')
            hdu0.header['KMAG'] = ('Unknown','[mag] K band magnitude from 2MASS')
            hdu0.header['KEPMAG'] = ('Unknown','[mag] Kepler magnitude (Kp) from KIC')
            hdu0.header['GRCOLOR'] = ('Unknown','[mag] (g-r) color, SDSS bands')
            hdu0.header['JKCOLOR'] = ('Unknown','[mag] (J-K) color, 2MASS bands')
            hdu0.header['GKCOLOR'] = ('Unknown','[mag] (g-K) color, SDSS g - 2MASS K')
            hdu0.header['TEFF'] = ('Unknown','[K] effective temperature from KIC')
            hdu0.header['LOGG'] = ('Unknown','[cm/s2] log10 surface gravity from KIC')
            hdu0.header['FEH'] = ('Unknown','[log10([Fe/H])] metallicity from KIC')
            hdu0.header['EBMINUSV'] = ('Unknown','[mag] E(B-V) redenning from KIC')
            hdu0.header['AV'] = ('Unknown','[mag] A_v extinction from KIC')
            hdu0.header['RADIUS'] = ('Unknown','[solar radii] stellar radius from KIC')
            hdu0.header['TMINDEX'] = ('Unknown','unique 2MASS catalog ID from KIC')
            hdu0.header['SCPID'] = ('Unknown','unique SCP processing ID from KIC')
            hdulist = pyfits.HDUList(hdu0)
        except:
            message = 'ERROR -- KEPCONVERT: cannot create primary extension in ' + outfile
            status = kepmsg.err(logfile,message,verbose)
## create the outfile HDU 1 extension

    if status == 0 and conversion == 'asc2fits':
        try:
            fitscol = []
            for i in range(ncol):
                fitscol.append(pytfits.Column(name=colnames[i],format='D',array=work[i]))
            fitscols = pyfits.ColDefs(fitscol)
            hdu1 = pyfits.BinTableHDU.from_columns(fitscols)
            hdulist.append(hdu1)
            hdu1.header['INHERIT'] = (True,'inherit primary keywords')
            hdu1.header['EXTNAME'] = ('LIGHTCURVE','name of extension')
            hdu1.header['EXTVER'] = (1,'extension version number')
            hdu1.header['TELESCOP'] = ('Kepler','telescope')
            hdu1.header['INSTRUME'] = ('Kepler photometer','detector type')
            hdu1.header['OBJECT'] = ('Unknown','string version of kepID')
            hdu1.header['KEPLERID'] = ('Unknown','unique Kepler target identifier')
            hdu1.header['RADESYS'] = ('Unknown','reference frame of celestial coordinates')
            hdu1.header['RA_OBJ'] = ('Unknown','[deg] right ascension from KIC')
            hdu1.header['DEC_OBJ'] = ('Unknown','[deg] declination from KIC')
            hdu1.header['EQUINOX'] = (2000.0,'equinox of celestial coordinate system')
            hdu1.header['TIMEREF'] = ('Unknown','barycentric correction applied to times')
            hdu1.header['TASSIGN'] = ('Unknown','where time is assigned')
            hdu1.header['TIMESYS'] = ('Unknown','time system is barycentric JD')
            hdu1.header['BJDREFI'] = (0.0,'integer part of BJD reference date')
            hdu1.header['BJDREFF'] = (0.0,'fraction of day in BJD reference date')
            hdu1.header['TIMEUNIT'] = ('Unknown','time unit for TIME, TSTART and TSTOP')
            hdu1.header['TSTART'] = (tstart,'observation start time in JD - BJDREF')
            hdu1.header['TSTOP'] = (tstop,'observation stop time in JD - BJDREF')
            hdu1.header['LC_START'] = (lc_start,'observation start time in MJD')
            hdu1.header['LC_END'] = (lc_end,'observation stop time in MJD')
            hdu1.header['TELAPSE'] = (tstop-tstart,'[d] TSTOP - TSTART')
            hdu1.header['LIVETIME'] = ('Unknown','[d] TELAPSE multiplied by DEADC')
            hdu1.header['EXPOSURE'] = ('Unknown','[d] time on source')
            hdu1.header['DEADC'] = ('Unknown','deadtime correction')
            hdu1.header['TIMEPIXR'] = ('Unknown','bin time beginning=0 middle=0.5 end=1')
            hdu1.header['TIERRELA'] = ('Unknown','[d] relative time error')
            hdu1.header['TIERABSO'] = ('Unknown','[d] absolute time error')
            hdu1.header['INT_TIME'] = ('Unknown','[s] photon accumulation time per frame')
            hdu1.header['READTIME'] = ('Unknown','[s] readout time per frame')
            hdu1.header['FRAMETIM'] = ('Unknown','[s] frame time (INT_TIME + READTIME)')
            hdu1.header['NUM_FRM'] = ('Unknown','number of frames per time stamp')
            hdu1.header['TIMEDEL'] = ('Unknown','[d] time resolution of data')
            hdu1.header['DATE-OBS'] = ('Unknown','TSTART as UT calendar date')
            hdu1.header['DATE-END'] = ('Unknown','TSTOP as UT calendar date')
            hdu1.header['BACKAPP'] = ('Unknown','background is subtracted')
            hdu1.header['DEADAPP'] = ('Unknown','deadtime applied')
            hdu1.header['VIGNAPP'] = ('Unknown','vignetting or collimator correction applied')
            hdu1.header['GAIN'] = ('Unknown','channel gain [electrons/count]')
            hdu1.header['READNOIS'] = ('Unknown','read noise [electrons]')
            hdu1.header['NREADOUT'] = ('Unknown','number of reads per cadence')
            hdu1.header['TIMSLICE'] = ('Unknown','time-slice readout sequence section')
            hdu1.header['MEANBLCK'] = ('Unknown','FSW mean black level [count]')
            hdu1.header['PDCSAPFL'] = ('Unknown','SAP PDC processing flags (bit code)')
            hdu1.header['PDCDIAFL'] = ('Unknown','DIA PDC processing flags (bit code)')
            hdu1.header['MISPXSAP'] = ('Unknown','no of optimal aperture pixels missing from SAP')
            hdu1.header['MISPXDIA'] = ('Unknown','no of optimal aperture pixels missing from DIA')
            hdu1.header['CROWDSAP'] = ('Unknown','crowding metric evaluated over SAP opt. ap.')
            hdu1.header['CROWDDIA'] = ('Unknown','crowding metric evaluated over DIA aperture')
        except:
            message = 'ERROR -- KEPCONVERT: cannot create light curve extension in ' + outfile
            status = kepmsg.err(logfile,message,verbose)

## history keyword in output file
    if status == 0 and conversion == 'asc2fits':
        status = kepkey.history(call,hdu0,outfile,logfile,verbose)

## filter data table

    if status == 0 and conversion == 'asc2fits':
        instr, status = kepio.filterNaN(hdulist,colnames[min(1,len(colnames)-1)],
                                        outfile,logfile,verbose)

## write output FITS file

    if status == 0 and conversion == 'asc2fits':
        hdulist.writeto(outfile,checksum=True)

## end time

    if status == 0:
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

    parser.add_argument('conversion', help='Type of data conversion', type=str,
                        choices=['fits2asc','asc2fits'])
    parser.add_argument('--columns', '-c', default='TIME,SAP_FLUX,SAP_FLUX_ERR',
                        dest='columns', help='Comma-delimited list of data columns', type=str)
    parser.add_argument('--baddata', action='store_false',
                        help='Output rows which have been flagged as questionable')
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    kepconvert(args.infile, args.outfile, args.conversion, args.columns,
               args.baddata, args.clobber, args.verbose, args.logfile, args.status)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepconvert.par")
    t = iraf.IrafTaskFactory(taskname="kepconvert", value=parfile, function=kepconvert)
