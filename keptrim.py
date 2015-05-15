import numpy, sys, time, pyfits, pylab, urllib
from numpy import *
from pyfits import *
from pylab import *
from matplotlib import *
import kepio, kepmsg, kepkey, kepstat

def keptrim(infile,outfile,kepid,column,row,imsize,clobber,verbose,logfile,status): 

# startup parameters

    status = 0

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPTRIM -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'kepid='+str(kepid)+' '
    call += 'column='+str(column)+' '
    call += 'row='+str(row)+' '
    call += 'imsize='+str(imsize)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPTRIM started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPTRIM: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    status = 0
    instr = pyfits.open(infile,mode='readonly',memmap=True)
    cards0 = instr[0].header.cards
    cards1 = instr[1].header.cards
    cards2 = instr[2].header.cards

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# identify the season of observation

    if status == 0:
        try:
            season = cards0['SEASON'].value
        except:
            season = 0

# retrieve column and row from KIC

        try:
            kic = FOVKepID(str(kepid))
            column = int(kic[98 + season * 5])
            row = int(kic[97 + season * 5])
        except:
            pass

# convert CCD column and row to image column and row

    if status == 0:
        if imsize % 2 == 0: imsize += 1
        crpix1p = cards2['CRPIX1P'].value
        crpix2p = cards2['CRPIX2P'].value
        crval1p = cards2['CRVAL1P'].value
        crval2p = cards2['CRVAL2P'].value
        cdelt1p = cards2['CDELT1P'].value
        cdelt2p = cards2['CDELT2P'].value        
        imcol = (column - crval1p) * cdelt1p + crpix1p - 1
        imrow = (row - crval2p) * cdelt2p + crpix2p - 1
        crval1p = column - imsize / 2 + 0.5
        crval2p = row - imsize / 2 + 0.5

# check subimage is contained inside the input image

    if status == 0:
        naxis1 = cards2['NAXIS1'].value
        naxis2 = cards2['NAXIS2'].value
        x1 = imcol - imsize / 2 + 0.5; x2 = x1 + imsize
        y1 = imrow - imsize / 2 + 0.5; y2 = y1 + imsize
        if x1 < 0 or y1 < 0 or x2 > naxis1 or y2 > naxis2:
            message =  'ERROR -- KEPTRIM: Requested pixel area falls outside of the pixel image in file ' + infile
            message += '. Make the pixel area smaller or relocate it''s center.'
            status = kepmsg.err(logfile,message,verbose)

# time series data

    if status == 0:
        time = instr[1].data.field('TIME')[:]
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

# resize time series

    if status == 0:
        raw_cnts = raw_cnts[:,y1:y2,x1:x2]
        flux = flux[:,y1:y2,x1:x2]
        flux_err = flux_err[:,y1:y2,x1:x2]
        flux_bkg = flux_bkg[:,y1:y2,x1:x2]
        flux_bkg_err = flux_bkg_err[:,y1:y2,x1:x2]
        cosmic_rays = cosmic_rays[:,y1:y2,x1:x2]

# reshape time series images

    if status == 0:
        isize = numpy.shape(flux)[0]
        jsize = numpy.shape(flux)[1]
        ksize = numpy.shape(flux)[2]
        raw_cnts = numpy.reshape(raw_cnts,(isize,jsize*ksize))
        flux = numpy.reshape(flux,(isize,jsize*ksize))
        flux_err = numpy.reshape(flux_err,(isize,jsize*ksize))
        flux_bkg = numpy.reshape(flux_bkg,(isize,jsize*ksize))
        flux_bkg_err = numpy.reshape(flux_bkg_err,(isize,jsize*ksize))
        cosmic_rays = numpy.reshape(cosmic_rays,(isize,jsize*ksize))
        
# pixel map data

    if status == 0:
        maskmap = array(instr[2].data[y1:y2,x1:x2])

# construct output primary extension

    if status == 0:
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

# construct output light curve extension

    if status == 0:
        coldim = '(' + str(imsize) + ',' + str(imsize) + ')'
        eformat = str(imsize*imsize) + 'E'
        jformat = str(imsize*imsize) + 'J'
        kformat = str(imsize*imsize) + 'K'
        col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
        col2 = Column(name='TIMECORR',format='E',unit='d',array=timecorr)
        col3 = Column(name='CADENCENO',format='J',array=cadenceno)
        col4 = Column(name='RAW_CNTS',format=jformat,unit='count',dim=coldim,array=raw_cnts)
        col5 = Column(name='FLUX',format=eformat,unit='e-/s',dim=coldim,array=flux)
        col6 = Column(name='FLUX_ERR',format=eformat,unit='e-/s',dim=coldim,array=flux_err)
        col7 = Column(name='FLUX_BKG',format=eformat,unit='e-/s',dim=coldim,array=flux_bkg)
        col8 = Column(name='FLUX_BKG_ERR',format=eformat,unit='e-/s',dim=coldim,array=flux_bkg_err)
        col9 = Column(name='COSMIC_RAYS',format=eformat,unit='e-/s',dim=coldim,array=cosmic_rays)
        col10 = Column(name='QUALITY',format='J',array=quality)
        col11 = Column(name='POS_CORR1',format='E',unit='pixel',array=pos_corr1)
        col12 = Column(name='POS_CORR2',format='E',unit='pixel',array=pos_corr2)
        cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
        hdu1 = new_table(cols)
        for i in range(len(cards1)):
            try:
                if cards1[i].key not in hdu1.header.keys():
                    hdu1.header.update(cards1[i].key, cards1[i].value, cards1[i].comment)
                else:
                    hdu1.header.cards[cards1[i].key].comment = cards1[i].comment
            except:
                pass
        hdu1.header.update('1CRV4P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV4P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX4',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX4',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV5P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV5P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX5',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX5',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV6P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV6P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX6',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX6',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV7P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV7P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX7',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX7',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV8P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV8P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX8',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX8',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        hdu1.header.update('1CRV9P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('2CRV9P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu1.header.update('1CRPX9',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu1.header.update('2CRPX9',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        outstr.append(hdu1)

# construct output mask bitmap extension

    if status == 0:
        hdu2 = ImageHDU(maskmap)
        for i in range(len(cards2)):
            try:
                if cards2[i].key not in hdu2.header.keys():
                    hdu2.header.update(cards2[i].key, cards2[i].value, cards2[i].comment)
                else:
                    hdu2.header.cards[cards2[i].key].comment = cards2[i].comment
            except:
                pass
        hdu2.header.update('NAXIS1',imsize,'')
        hdu2.header.update('NAXIS2',imsize,'')
        hdu2.header.update('CRVAL1P',crval1p,'[pixel] detector coordinate at reference pixel')
        hdu2.header.update('CRVAL2P',crval2p,'[pixel] detector coordinate at reference pixel')
        hdu2.header.update('CRPIX1',(imsize + 1) / 2,'[pixel] reference pixel along image axis 1')
        hdu2.header.update('CRPIX2',(imsize + 1) / 2,'[pixel] reference pixel along image axis 2')
        outstr.append(hdu2)

# write output file

    if status == 0:
        outstr.writeto(outfile,checksum=True)

# close input structure

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

# end time

    kepmsg.clock('KEPTRIM finished at',logfile,verbose)


# -------------------------------------
# KIC retrieval based upon KepID

def FOVKepID(id):

# build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&kic_kepler_id=' + id
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=CSV'

# retrieve results from MAST

    out = ''
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')

    return out


# -------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description=
        'Trim unwanted pixels from a Target Pixel File')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('outfile', help='Name of output target pixel file', type=str)
    parser.add_argument('--kepid', help='Kepler ID number from the Kepler Input Catalog', dest='kepid', type=int)
    parser.add_argument('--column', help='CCD column number of the target', dest='column', type=int)
    parser.add_argument('--row', help='CCD row number of the target', dest='row', type=int)
    parser.add_argument('--imsize', help='Number of pixels to extract in both row and column dimensions', dest='imsize', type=int)
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='keptrim.log', 
        dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    keptrim(args.infile,args.outfile,args.kepid,args.column,args.row,args.imsize,
            args.clobber,args.verbose,args.logfile,args.status)
    
    

else:
    from pyraf import iraf
    
    
    parfile = iraf.osfn("kepler$keptrim.par")
    t = iraf.IrafTaskFactory(taskname="keptrim", value=parfile, function=keptrim)

