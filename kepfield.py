import pylab, numpy, scipy
from pylab import *
from matplotlib import *
from numpy import *
from scipy import *
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey, kepplot, kepfit, keparray, kepstat
import sys, time, re, math, glob, urllib

# -----------------------------------------------------------
# core code

def kepfield(infile,plotfile,rownum,imscale,colmap,lcolor,srctab,verbose,logfile,status,cmdLine=False): 

# input arguments

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFIELD -- '
    call += 'infile='+infile+' '
    call += 'plotfile='+plotfile+' '
    call += 'rownum='+str(rownum)+' '
    call += 'imscale='+imscale+' '
    call += 'colmap='+colmap+' '
    call += 'lcolor='+lcolor+' '
    srct = 'n'
    if (srctab): srct = 'y'
    call += 'srctab='+srct+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPFIELD started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# reference color map

    if colmap == 'browse':
        status = cmap_plot(cmdLine)

# open TPF FITS file

    if status == 0:
        try:
            kepid, channel, skygroup, module, output, quarter, season, \
                ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
                kepio.readTPF(infile,'TIME',logfile,verbose)
        except:
            message = 'ERROR -- KEPFIELD: is %s a Target Pixel File? ' % infile
            status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, tcorr, status = \
            kepio.readTPF(infile,'TIMECORR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, cadno, status = \
            kepio.readTPF(infile,'CADENCENO',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, fluxpixels, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, errpixels, status = \
            kepio.readTPF(infile,'FLUX_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)

# read mask defintion data from TPF file

    if status == 0:
        maskimg, pixcoord1, pixcoord2, status = kepio.readMaskDefinition(infile,logfile,verbose)

# observed or simulated data?
        
    if status == 0:
        coa = False
        instr = pyfits.open(infile,mode='readonly',memmap=True)
        filever, status = kepkey.get(infile,instr[0],'FILEVER',logfile,verbose)
        if filever == 'COA': coa = True

# print target data

    if status == 0 and verbose:
        print ''
        print '      KepID: %s' % kepid
        print '        BJD: %.2f' % (barytime[rownum-1] + 2454833.0)
        print ' RA (J2000): %s' % ra
        print 'Dec (J2000):  %s' % dec
        print '     KepMag:  %s' % kepmag
        print '   SkyGroup:   %2s' % skygroup
        print '     Season:   %2s' % str(season)
        print '    Channel:   %2s' % channel
        print '     Module:   %2s' % module
        print '     Output:    %1s' % output
        print ''

# is this a good row with finite timestamp and pixels?

    if status == 0:
        if not numpy.isfinite(barytime[rownum-1]) or not numpy.nansum(fluxpixels[rownum-1,:]):
            message = 'ERROR -- KEPFIELD: Row ' + str(rownum) + ' is a bad quality timestamp'
            status = kepmsg.err(logfile,message,verbose)

# construct input pixel image

    if status == 0:
        flux = fluxpixels[rownum-1,:]

# image scale and intensity limits of pixel data

    if status == 0:
        flux_pl, zminfl, zmaxfl = kepplot.intScale1D(flux,imscale)
        n = 0
        imgflux_pl = empty((ydim+2,xdim+2))
        for i in range(ydim+2):
            for j in range(xdim+2):
                imgflux_pl[i,j] = numpy.nan
        for i in range(ydim):
            for j in range(xdim):
                imgflux_pl[i+1,j+1] = flux_pl[n]
                n += 1
        
# cone search around target coordinates using the MAST target search form 

    if status == 0:
        dr = max([ydim+2,xdim+2]) * 4.0
        kepid,ra,dec,kepmag = MASTRADec(float(ra),float(dec),dr,srctab)

# convert celestial coordinates to detector coordinates

    if status == 0:
        sx = numpy.array([])
        sy = numpy.array([])
        inf, status = kepio.openfits(infile,'readonly',logfile,verbose)
        try:
            crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, pc, status = \
                kepkey.getWCSs(infile,inf['APERTURE'],logfile,verbose) 
            crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status = \
                kepkey.getWCSp(infile,inf['APERTURE'],logfile,verbose)     
            for i in range(len(kepid)):
                dra = (ra[i] - crval1) * math.cos(math.radians(dec[i])) / cdelt1
                ddec = (dec[i] - crval2) / cdelt2
                if coa:
                    sx = numpy.append(sx,-(pc[0,0] * dra + pc[0,1] * ddec) + crpix1 + crval1p - 1.0)
                else:
                    sx = numpy.append(sx,pc[0,0] * dra + pc[0,1] * ddec + crpix1 + crval1p - 1.0) 
                sy = numpy.append(sy,pc[1,0] * dra + pc[1,1] * ddec + crpix2 + crval2p - 1.0)
        except:
            message = 'ERROR -- KEPFIELD: Non-compliant WCS information within file %s' % infile
            status = kepmsg.err(logfile,message,verbose)    

# plot style

    if status == 0:
        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': 48,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': 20,
                      'ytick.labelsize': 20}
            pylab.rcParams.update(params)
        except:
            pass
        pylab.figure(figsize=[10,10])
        pylab.clf()
            
# pixel limits of the subimage

    if status == 0:
        ymin = copy(float(row))
        ymax = ymin + ydim
        xmin = copy(float(column))
        xmax = xmin + xdim

# plot limits for flux image

    if status == 0:
        ymin = float(ymin) - 1.5
        ymax = float(ymax) + 0.5
        xmin = float(xmin) - 1.5
        xmax = float(xmax) + 0.5

# plot the image window
        
    if status == 0:
        ax = pylab.axes([0.1,0.11,0.88,0.88])
        pylab.imshow(imgflux_pl,aspect='auto',interpolation='nearest',origin='lower',
                     vmin=zminfl,vmax=zmaxfl,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
        pylab.gca().set_autoscale_on(False)
        labels = ax.get_yticklabels()
        setp(labels, 'rotation', 90)
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.xlabel('Pixel Column Number', {'color' : 'k'})
        pylab.ylabel('Pixel Row Number', {'color' : 'k'})

# plot mask borders

    if status == 0:
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,1,lcolor,'--',0.5)

# plot aperture borders

    if status == 0:
        kepplot.borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,2,lcolor,'-',4.0)

# list sources

    if status == 0:
        print 'Column    Row  RA J2000 Dec J2000    Kp    Kepler ID'
        print '----------------------------------------------------'
        for i in range(len(sx)-1,-1,-1):
            if sx[i] >= xmin and sx[i] < xmax and sy[i] >= ymin and sy[i] < ymax:
                if kepid[i] != 0 and kepmag[i] != 0.0:
                    print '%6.1f %6.1f %9.5f  %8.5f %5.2f KIC %d' % \
                        (float(sx[i]),float(sy[i]),float(ra[i]),float(dec[i]),float(kepmag[i]),int(kepid[i]))
                elif kepid[i] != 0 and kepmag[i] == 0.0:
                    print '%6.1f %6.1f %9.5f  %8.5f       KIC %d' % \
                        (float(sx[i]),float(sy[i]),float(ra[i]),float(dec[i]),int(kepid[i]))
                else:
                    print '%6.1f %6.1f %9.5f  %8.5f' % (float(sx[i]),float(sy[i]),float(ra[i]),float(dec[i]))

# plot sources

    if status == 0:
        for i in range(len(sx)-1,-1,-1):
            if kepid[i] != 0 and kepmag[i] != 0.0:
                size = max(array([80.0,80.0 + (2.5**(18.0 - max(12.0,float(kepmag[i])))) * 250.0]))
                pylab.scatter(sx[i],sy[i],s=size,facecolors='g',edgecolors='k',alpha=0.4)
            else:
                pylab.scatter(sx[i],sy[i],s=80,facecolors='r',edgecolors='k',alpha=0.4)

# render plot

    if status == 0 and len(plotfile) > 0 and plotfile.lower() != 'none':
        pylab.savefig(plotfile)
    if status == 0:
        if cmdLine: 
            pylab.show(block=True)
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
# stop time

    kepmsg.clock('\nKEPFIELD ended at',logfile,verbose)

    return

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot(cmdLine):

    pylab.figure(figsize=[5,10])
    a=outer(ones(10),arange(0,1,0.01))
    subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        print m
        subplot(l,1,i+1)
        pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        pylab.text(100.85,0.5,m,fontsize=10)

# render plot

    if cmdLine: 
        pylab.show(block=True)
    else: 
        pylab.ion()
        pylab.plot([])
        pylab.ioff()
	
    status = 1
    return status

# -------------------------------------
# detector location retrieval based upon RA and Dec

def MASTRADec(ra,dec,darcsec,srctab):

# coordinate limits

    darcsec /= 3600.0
    ra1 = ra - darcsec / cos(dec * pi / 180)
    ra2 = ra + darcsec / cos(dec * pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

# build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&masterRA=' + str(ra1) + '..' + str(ra2)
    url += '&masterDec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=10000'
    url += '&verb=3'
    url += '&outputformat=CSV'

# retrieve results from MAST

    if srctab:
        try:
            lines = urllib.urlopen(url)
        except:
            message = 'WARNING -- KEPFIELD: Cannot retrieve data from MAST'
            status = kepmsg.warn(logfile,message)
            lines = ''
    else:
        lines = ''

# collate nearby sources

    kepid = []
    kepmag = []
    ra = []
    dec = []
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')
            r,d = sex2dec(out[0],out[1])
            try:
                if out[-22] != 'Possible_artifact': kepid.append(int(out[2]))
            except:
                if out[-22] != 'Possible_artifact': kepid.append(0)
            try:
                if out[-22] != 'Possible_artifact': kepmag.append(float(out[42]))
            except:
                if out[-22] != 'Possible_artifact': kepmag.append(0.0)
            if out[-22] != 'Possible_artifact': ra.append(r)
            if out[-22] != 'Possible_artifact': dec.append(d)
    kepid = array(kepid)
    kepmag = array(kepmag)
    ra = array(ra)
    dec = array(dec)

    return kepid,ra,dec,kepmag

# -----------------------------------
# convert sexadecimal hours to decimal degrees

def sex2dec(ra,dec):

    ra = re.sub('\s+','|',ra.strip())
    ra = re.sub(':','|',ra.strip())
    ra = re.sub(';','|',ra.strip())
    ra = re.sub(',','|',ra.strip())
    ra = re.sub('-','|',ra.strip())
    ra = ra.split('|')
    outra = (float(ra[0]) + float(ra[1]) / 60 + float(ra[2]) / 3600) * 15.0

    dec = re.sub('\s+','|',dec.strip())
    dec = re.sub(':','|',dec.strip())
    dec = re.sub(';','|',dec.strip())
    dec = re.sub(',','|',dec.strip())
    dec = re.sub('-','|',dec.strip())
    dec = dec.split('|')
    if float(dec[0]) > 0.0:
        outdec = float(dec[0]) + float(dec[1]) / 60 + float(dec[2]) / 3600
    else:
        outdec = float(dec[0]) - float(dec[1]) / 60 - float(dec[2]) / 3600

    return outra, outdec


# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Displaying the target mask and local field sources')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input target pixel file', type=str)
    parser.add_argument('--plotfile', '-p', help='Name of output PNG plot file', default='None', dest='plotfile', type=str)
    parser.add_argument('--rownum', '-r', default=2200, help='Row number of image stored in infile', dest='rownum', type=int)
    parser.add_argument('--imscale', '-i', help='Type of image intensity scale', default='linear', dest='imscale', type=str,choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', '-c', help='Image colormap', default='YlOrBr', dest='cmap', type=str,choices=['Accent','Blues','BrBG','BuGn','BuPu','Dark2','GnBu','Greens','Greys','OrRd','Oranges','PRGn','Paired','Pastel1','Pastel2','PiYG','PuBu','PuBuGn','PuOr','PuRd','Purples','RdBu','RdGy','RdPu','RdYlBu','RdYlGn','Reds','Set1','Set2','Set3','Spectral','YlGn','YlGnBu','YlOrBr','YlOrRd','afmhot','autumn','binary','bone','brg','bwr','cool','copper','flag','gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_yarg','gnuplot','gnuplot2','gray','hot','hsv','jet','ocean','pink','prism','rainbow','seismic','spectral','spring','summer','terrain','winter','browse'])
    parser.add_argument('--lcolor', default='#0000ff', help='HTML color of data line within plot', type=str)
    parser.add_argument('--srctab', action='store_true', help='Extract data from target table at MAST?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepfieldphot.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepfield(args.infile,args.plotfile,args.rownum,args.imscale,args.cmap,args.lcolor,args.srctab,
             args.verbose,args.logfile,args.status,cmdLine)
    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepfield.par")
    t = iraf.IrafTaskFactory(taskname="kepfield", value=parfile, function=kepfield)
