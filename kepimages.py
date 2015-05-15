import time, numpy, sys, pyfits
from numpy import array
from pyfits import *
import kepio, kepmsg, kepkey, kepstat

def kepimages(infile,outfix,imtype,ranges,clobber,verbose,logfile,status): 

# startup parameters

    status = 0

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPIMAGES -- '
    call += 'infile='+infile+' '
    call += 'outfix='+outfix+' '
    call += 'imtype='+imtype+' '
    call += 'ranges='+str(ranges)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPIMAGES started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# open input file

    status = 0
    print ' '
    instr = pyfits.open(infile,mode='readonly',memmap=True)
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# ingest time series data

    if status == 0:
        time = instr[1].data.field('TIME')[:] + 2454833.0
        timecorr = instr[1].data.field('TIMECORR')[:]
        cadenceno = instr[1].data.field('CADENCENO')[:]
        raw_cnts = instr[1].data.field('RAW_CNTS')[:]
        flux = instr[1].data.field('FLUX')[:]
        flux_err = instr[1].data.field('FLUX_ERR')[:]
        flux_bkg = instr[1].data.field('FLUX_BKG')[:]
        flux_bkg_err = instr[1].data.field('FLUX_BKG_ERR')[:]
        cosmic_rays = instr[1].data.field('COSMIC_RAYS')[:]
        quality = instr[1].data.field('QUALITY')[:]
        pos_corr1 = instr[1].data.field('POS_CORR1')[:]
        pos_corr2 = instr[1].data.field('POS_CORR2')[:]

# choose output image

    if status == 0:
        if imtype.lower() == 'raw_cnts':
            outim = raw_cnts
        elif imtype.lower() == 'flux_err':
            outim = flux_err
        elif imtype.lower() == 'flux_bkg':
            outim = flux_bkg
        elif imtype.lower() == 'flux_bkg_err':
            outim = flux_bkg_err
        elif imtype.lower() == 'cosmic_rays':
            outim = cosmic_rays
        else:
            outim = flux

# identify images to be exported

    if status == 0:
        tim = array([]); dat = array([]); err = array([])
        tstart, tstop, status = kepio.timeranges(ranges,logfile,verbose)
    if status == 0:
        cadencelis, status = kepstat.filterOnRange(time,tstart,tstop)

# provide name for each output file and clobber if file exists

    if status == 0:
        for cadence in cadencelis:
            outfile = outfix + '_BJD%.4f' % time[cadence] + '.fits'
            if clobber and status == 0: status = kepio.clobber(outfile,logfile,verbose)
            if kepio.fileexists(outfile) and status == 0: 
                message = 'ERROR -- KEPIMAGES: ' + outfile + ' exists. Use --clobber'
                status = kepmsg.err(logfile,message,True)

# construct output primary extension

    if status == 0:
        ncad = 0
        for cadence in cadencelis:
            outfile = outfix + '_BJD%.4f' % time[cadence] + '.fits'
            hdu0 = pyfits.PrimaryHDU()
            for i in range(len(cards0)):
                try:
                    if cards0[i].key not in hdu0.header.keys():
                        hdu0.header.update(cards0[i].key, cards0[i].value, cards0[i].comment)
                    else:
                        hdu0.header.cards[cards0[i].key].comment = cards0[i].comment
                except:
                    pass
            status = kepkey.history(call,hdu0,outfile,logfile,verbose)
            outstr = HDUList(hdu0)

# construct output image extension

            hdu1 = ImageHDU(flux[cadence])
            for i in range(len(cards2)):
                try:
                    if cards2[i].key not in hdu1.header.keys():
                        hdu1.header.update(cards2[i].key, cards2[i].value, cards2[i].comment)
                except:
                    pass
            for i in range(len(cards1)):
                if (cards1[i].key not in hdu1.header.keys() and
                    cards1[i].key[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                              '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                              '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                              '12PC','21PC','22PC','WCSN','TFIE']):
                    hdu1.header.update(cards1[i].key, cards1[i].value, cards1[i].comment)
            try:
                int_time = cards1['INT_TIME'].value
            except:
                kepmsg.warn(logfile,'WARNING -- KEPIMAGES: cannot find INT_TIME keyword')
            try:
                frametim = cards1['FRAMETIM'].value
            except:
                kepmsg.warn(logfile,'WARNING -- KEPIMAGES: cannot find FRAMETIM keyword')
            try:
                num_frm = cards1['NUM_FRM'].value
            except:
                kepmsg.warn(logfile,'WARNING -- KEPIMAGES: cannot find NUM_FRM keyword')
            hdu1.header.update('EXTNAME','IMAGE','name of extension')
            try:
                hdu1.header.update('TELAPSE',frametim * num_frm,'[s] elapsed time for exposure')
            except:
                hdu1.header.update('TELAPSE',-999,'[s] elapsed time for exposure')
            try:
                hdu1.header.update('LIVETIME',int_time * num_frm,'[s] TELASPE multiplied by DEADC')
            except:
                hdu1.header.update('LIVETIME',-999,'[s] TELASPE multiplied by DEADC')
            try:
                hdu1.header.update('EXPOSURE',int_time * num_frm,'[s] time on source')
            except:
                hdu1.header.update('EXPOSURE',-999,'[s] time on source')
            try:
                hdu1.header.update('MIDTIME',time[cadence],'[BJD] mid-time of exposure')
            except:
                hdu1.header.update('MIDTIME',-999,'[BJD] mid-time of exposure')
            try:
                hdu1.header.update('TIMECORR',timecorr[cadence],'[d] barycenter - timeslice correction')
            except:
                hdu1.header.update('TIMECORR',-999,'[d] barycenter - timeslice correction')
            try:
                hdu1.header.update('CADENCEN',cadenceno[cadence],'unique cadence number')
            except:
                hdu1.header.update('CADENCEN',-999,'unique cadence number')
            try:
                hdu1.header.update('QUALITY',quality[cadence],'pixel quality flag')
            except:
                hdu1.header.update('QUALITY',-999,'pixel quality flag')
            try:
                if True in numpy.isfinite(cosmic_rays[cadence]):
                    hdu1.header.update('COSM_RAY',True,'cosmic ray detected?')
                else:
                    hdu1.header.update('COSM_RAY',False,'cosmic ray detected?')
            except:
                hdu1.header.update('COSM_RAY',-999,'cosmic ray detected?')
            try:
                pc1 = str(pos_corr1[cadence])
                pc2 = str(pos_corr2[cadence])
                hdu1.header.update('POSCORR1',pc1,'[pix] column position correction')
                hdu1.header.update('POSCORR2',pc2,'[pix] row position correction')
            except:
                hdu1.header.update('POSCORR1',-999,'[pix] column position correction')
                hdu1.header.update('POSCORR2',-999,'[pix] row position correction')
            outstr.append(hdu1)

# write output file

            if status == 0:
                outstr.writeto(outfile,checksum=True)
                ncad += 1
                txt  = '\r%3d%% ' % (float(ncad) / float(len(cadencelis)) * 100.0)
                txt += '%s ' % outfile
                sys.stdout.write(txt)
                sys.stdout.flush()

# close input structure

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    
        print '\n'

# end time

    kepmsg.clock('KEPIMAGES finished at',logfile,verbose)

# -------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Export images within a Target Pixel File to a series of FITS image files')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('outfix', help='Prefix of name for each output image file', type=str)
    parser.add_argument('--imtype', default='FLUX', dest='imtype', help='Image type', type=str,
                        choices=['RAW_CNTS','FLUX','FLUX_ERR','FLUX_BKG','FLUX_BKG_ERR','COSMIC_RAYS'])
    parser.add_argument('--ranges', help='Time ranges for output [BJD]', type=str)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', help='Name of ascii log file', default='kepimages.log', 
        dest='logfile', type=str)
    parser.add_argument('--status', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    
    kepimages(args.infile,args.outfix,args.imtype,args.ranges,args.clobber,args.verbose,args.logfile,args.status)

else:
    from pyraf import iraf
        
    parfile = iraf.osfn("kepler$kepimages.par")
    t = iraf.IrafTaskFactory(taskname="kepimages", value=parfile, function=kepimages)
