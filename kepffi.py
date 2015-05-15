
import pylab, numpy, pyfits
from pylab import *
from matplotlib import *
from numpy import *
import kepio, kepmsg, kepkey
import sys, urllib, time, re

# global variables

ffifile = False; aperfile = False; maskfile = 'KeplerFFI.txt'
plotfile = 'KeplerFFI.png'; pimg = None; mask = []; zscale = False
xmin = 0.0; xmax = 1000.0; ymin = 0.0; ymax = 1000.0; zmin = False; zmax = False
kepid = ''; ra = ''; dec = ''; kepmag = ''; season = ''; quarter = -1
skygroup = ''; channel = ''; module = ''; output = ''; column = ''; row = ''
colmap='jet'; aid = None; cid = None; did = None; eid = None; fid = None
pkepmag = None; pkepid = None; pra = None; pdec = None; cmdLine = False

# -----------------------------------------------------------
# core code

def kepffi(ffifile,kepid,ra,dec,aperfile,imin,imax,iscale,cmap,npix,
            verbose,logfile,status,cmdLine=False): 

    global pimg, zscale, zmin, zmax, xmin, xmax, ymin, ymax, quarter
    global kepmag, skygroup, season, channel
    global module, output, row, column, maskfile, plotfile
    global pkepid, pkepmag, pra, pdec, colmap, mask

# input arguments

    status = 0
    seterr(all="ignore") 
    maskfile = 'kepffi-' + str(kepid) + '.txt'
    plotfile = 'kepffi-' + str(kepid) + '.png'
    zmin = imin; zmax = imax; zscale = iscale; colmap = cmap

# logg the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFFI -- '
    call += 'ffifile='+ffifile+' '
    call += 'kepid='+str(kepid)+' '
    call += 'ra='+ra+' '
    call += 'dec='+dec+' '
    call += 'aperfile='+aperfile+' '
    call += 'imin='+str(imin)+' '
    call += 'imax='+str(imax)+' '
    call += 'iscale='+str(iscale)+' '
    call += 'cmap'+str(cmap)+' '
    call += 'npix='+str(npix)+' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPFFI started at',logfile,verbose)

# reference color map

    if cmap == 'browse':
        status = cmap_plot(cmdLine)

# open existing mask file

    if kepio.fileexists(aperfile):
        lines, status = kepio.openascii(aperfile,'r',logfile,verbose)
        for line in lines:
            line = line.strip().split('|')
            y0 = int(line[3])
            x0 = int(line[4])
            pixels = line[5].split(';')
            for pixel in pixels:
                m = y0 + int(pixel.split(',')[0])
                n = x0 + int(pixel.split(',')[1])
                mask.append(str(m)+','+str(n))
        status = kepio.closeascii(lines,logfile,verbose)

# RA and Dec conversion

    if kepid == 'None' or kepid == 'none' or kepid.strip() == '':
        try:
            mra = float(ra)
            mdec = float(dec)
        except:
            try:
                mra,mdec = sex2dec(ra,dec)
            except:
                txt = 'ERROR -- no sensible RA and Dec coordinates provided'
                sys.exit(txt)

# open FFI FITS file

    if status == 0:
        ffi, status = openfits(ffifile,'readonly')
        try:
            quarter = ffi[0].header['QUARTER']
        except:
            try:
                dateobs = ffi[0].header['DATE-OBS']
                if dateobs == '2009-04-24': quarter = 0
                if dateobs == '2009-04-25': quarter = 0
                if dateobs == '2009-04-26': quarter = 0
                if dateobs == '2009-06-19': quarter = 2
                if dateobs == '2009-08-19': quarter = 2
                if dateobs == '2009-09-17': quarter = 2
                if dateobs == '2009-10-19': quarter = 3
                if dateobs == '2009-11-18': quarter = 3
                if dateobs == '2009-12-17': quarter = 3
            except:
                txt  = 'ERROR -- cannot determine quarter when FFI was taken. Either a\n'
                txt += 'QUARTER or DATE-OBS keyword is expected in the primary header'
                sys.exit(txt)
        if quarter == 0: quarter = 1
        if quarter < 0:
            txt  = 'ERROR -- cannot determine quarter from FFI. Try downloading a new\n'
            txt += 'version of KeplerFFI.py from http://keplergo.arc.nasa.gov'
            sys.exit()
    if int(quarter) == 0:
        season = 3
    else:
        season = (int(quarter) - 2) % 4

# locate target in MAST

        try:
            int(kepid)
            kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column \
                = MASTKepID(kepid,season)
            pkepmag = kepmag; pkepid = kepid
        except:
            kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column \
                = MASTRADec(mra,mdec,8.0,season)
            ra,dec = dec2sex(ra,dec)
        pra = ra; pdec = dec
        print kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column
       
# read and close FFI FITS file

        img, status = readimage(ffi,int(channel))
        status = closefits(ffi)

# print target data

        print ''
        print '      KepID:  %s' % kepid
        print ' RA (J2000):  %s' % ra
        print 'Dec (J2000): %s' % dec
        print '     KepMag:  %s' % kepmag
        print '   SkyGroup:    %2s' % skygroup
        print '     Season:    %2s' % str(season)
        print '    Channel:    %2s' % channel
        print '     Module:    %2s' % module
        print '     Output:     %1s' % output
        print '     Column:  %4s' % column
        print '        Row:  %4s' % row
        print ''

# subimage of channel for plot

        ymin = int(max([int(row)-npix/2,0]))
        ymax = int(min([int(row)+npix/2+1,img.shape[0]]))
        xmin = int(max([int(column)-npix/2,0]))
        xmax = int(min([int(column)+npix/2+1,img.shape[1]]))

# intensity scale

        nstat = 2; pixels = []
        for i in range(ymin,ymax+1):
            for j in range(xmin,xmax+1):
                pixels.append(img[i,j])
        pixels = array(sort(pixels),dtype=float32)
        if int(float(len(pixels)) / 10 + 0.5) > nstat:
            nstat = int(float(len(pixels)) / 10 + 0.5)
        if not zmin:
            zmin = median(pixels[:nstat])
        if not zmax:
            zmax = median(pixels[-nstat:])
        if 'log' in zscale:
            img = log10(img)
            zmin = log10(zmin)
            zmax = log10(zmax)
        if ('sq' in zscale):
            img = sqrt(img)
            zmin = sqrt(zmin)
            zmax = sqrt(zmax)
        pimg = img[ymin:ymax,xmin:xmax]

# plot limits

        ymin = float(ymin) - 0.5
        ymax = float(ymax) - 0.5
        xmin = float(xmin) - 0.5
        xmax = float(xmax) - 0.5

# plot style

        try:
            params = {'backend': 'png',
                      'axes.linewidth': 2.5,
                      'axes.labelsize': 24,
                      'axes.font': 'sans-serif',
                      'axes.fontweight' : 'bold',
                      'text.fontsize': 12,
                      'legend.fontsize': 12,
                      'xtick.labelsize': 16,
                      'ytick.labelsize': 16}
            pylab.rcParams.update(params)
        except:
            pass

    if status == 0:
        pylab.figure(figsize=[10,7])
        plotimage(cmdLine)

# render plot

        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
    return

# -----------------------------------------------------------
# plot channel image

def plotimage(cmdLine):

    global aid, cid, did, eid, fid

# print FFI and source location data on plot
    
    pylab.clf()
    pylab.axes([0.73,0.09,0.25,0.4])
    pylab.text(0.1,1.0,'      KepID: %s' % pkepid,fontsize=12)
    pylab.text(0.1,0.9,' RA (J2000): %s' % pra,fontsize=12)
    pylab.text(0.1,0.8,'Dec (J2000): %s' % pdec,fontsize=12)
    pylab.text(0.1,0.7,'     KepMag: %s' % pkepmag,fontsize=12)
    pylab.text(0.1,0.6,'   SkyGroup: %2s' % skygroup,fontsize=12)
    pylab.text(0.1,0.5,'     Season: %2s' % str(season),fontsize=12)
    pylab.text(0.1,0.4,'    Channel: %2s' % channel,fontsize=12)
    pylab.text(0.1,0.3,'     Module: %2s' % module,fontsize=12)
    pylab.text(0.1,0.2,'     Output: %1s' % output,fontsize=12)
    pylab.text(0.1,0.1,'     Column: %4s' % column,fontsize=12)
    pylab.text(0.1,0.0,'        Row: %4s' % row,fontsize=12)
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    xlim(0.0,1.0)
    ylim(-0.05,1.12)

# clear button

    pylab.axes([0.73,0.87,0.25,0.09])
    pylab.text(0.5,0.5,'CLEAR',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    aid = connect('button_press_event',clicker1)

# dump custom aperture to file button

    pylab.axes([0.73,0.77,0.25,0.09])
    pylab.text(0.5,0.5,'DUMP',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    cid = connect('button_press_event',clicker3)

# print window to png file button

    pylab.axes([0.73,0.67,0.25,0.09])
    pylab.text(0.5,0.5,'PRINT',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    did = connect('button_press_event',clicker4)

# print close plot file button

    pylab.axes([0.73,0.57,0.25,0.09])
    pylab.text(0.5,0.5,'CLOSE',fontsize=24,weight='heavy',
               horizontalalignment='center',verticalalignment='center')
    pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    pylab.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    xlim(0.0,1.0)
    ylim(0.0,1.0)
    eid = connect('button_press_event',clicker5)

# plot the image window

    ax = pylab.axes([0.08,0.09,0.63,0.88])
    pylab.subplots_adjust(0.06,0.1,0.93,0.88)
    pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
    labels = ax.get_yticklabels()
    setp(labels, 'rotation', 90)
    imshow(pimg,aspect='auto',interpolation='nearest',origin='lower',
           vmin=zmin,vmax=zmax,extent=(xmin,xmax,ymin,ymax),cmap=colmap)
    pylab.gca().set_autoscale_on(False)
    xlabel('Pixel Column Number', {'color' : 'k'})
    ylabel('Pixel Row Number', {'color' : 'k'})
    grid()

# plot the mask

    if colmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                'gray','pink','RdGy']:
        sqcol = 'g'
        alpha = 0.5
    else:
        sqcol = '#ffffee'
        alpha = 0.8
    for pixel in mask:
        m = int(pixel.split(',')[0])
        n = int(pixel.split(',')[1])
        y = [m-0.5,m+0.5,m+0.5,m-0.5,m-0.5]
        x = [n-0.5,n-0.5,n+0.5,n+0.5,n-0.5]
        pylab.fill(x,y,sqcol,alpha=alpha,ec=sqcol)
#        print x,y
    fid = connect('key_press_event',clicker6)

# render plot

    if cmdLine: 
        pylab.show()
    else: 
        pylab.ion()
        pylab.plot([])
        pylab.ioff()	

    return

# -----------------------------------------------------------
# target data retrieval from MAST based upon KepID

def MASTKepID(id,season):

    global skygroup, column, row

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
    if len(out) > 0:
        kepid = out[2]
        ra = out[6]
        dec = out[7]
        kepmag = out[42]
        skygroup = out[73]
        channel = out[94 + season * 5]
        module = out[95 + season * 5]
        output = out[96 + season * 5]
        row = out[97 + season * 5]
        column = out[98 + season * 5]
    else:
        txt = 'ERROR -- no target found with KepID %s' % id
        sys.exit(txt)

    return kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column

# -------------------------------------
# detector location retrieval based upon RA and Dec

def MASTRADec(ra,dec,darcsec,season):

    global skygroup, column, row

# WCS data

    cd1_1 = 0.000702794927969
    cd1_2 = -0.000853190160515
    cd2_1 = -0.000853190160515
    cd2_2 = -0.000702794927969
    cd = array([[cd1_1,cd1_2],[cd2_1,cd2_2]])
    cd = linalg.inv(cd)

# coordinate limits

    x1 = 1.0e30
    x2 = x1
    darcsec /= 3600.0
    ra1 = ra - darcsec / 15.0 / cos(dec * pi / 180)
    ra2 = ra + darcsec / 15.0 / cos(dec * pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

# build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&kic_degree_ra=' + str(ra1) + '..' + str(ra2)
    url += '&kic_dec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=CSV'

# retrieve results from MAST: nearest KIC source to supplied coordinates

    z = ''
    x = 1.0e30
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')
            r = (float(out[6].split(' ')[0]) + \
                float(out[6].split(' ')[1]) / 60.0 + \
                float(out[6].split(' ')[2]) / 3600.0) * 15.0
            d = float(out[7].split(' ')[0]) + \
                float(out[7].split(' ')[1]) / 60.0 + \
                float(out[7].split(' ')[2]) / 3600.0
            a = sqrt((abs(r - ra) / 15.0 / cos(d * pi / 180))**2 + abs(d - dec)**2)
            if a < x:
                x = a
                z = line.split(',')

    if len(z) > 0:
        kepid = None
        kepmag = None
        skygroup = out[73]
        channel = out[94 + season * 5]
        module = out[95 + season * 5]
        output = out[96 + season * 5]
    else:
        txt = 'ERROR -- row and column could not be calculated. Is location on silicon?'
        sys.exit(txt)
            
# convert coordinates to decimal for the two targets, determine distance from input

    zra,zdec = sex2dec(z[6],z[7])
    dra = zra - ra
    ddec = zdec - dec
    drow = cd[0,0] * dra + cd[0,1] * ddec
    dcol = cd[1,0] * dra + cd[1,1] * ddec
   
# pixel coordinates of the nearest KIC target

    row = z[97 + season * 5]
    column = z[98 + season * 5]

# pixel coordinate of target

    row = str(int(float(row) + drow + 0.5))
    column = str(int(float(column) + dcol + 0.5))

    return kepid,ra,dec,kepmag,skygroup,channel,module,output,row,column

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

# -----------------------------------
# convert decimal RA and Dec to sexagesimal

def dec2sex(ra,dec):

    if ra < 0.0 or ra > 360.0 or dec < -90.0 or dec > 90.0: 
        txt = 'ERROR -- badly defined RA and Dec provided'
        sys.exit()

    tmp = ra / 15
    ra_h = str(int(tmp))
    tmp = (tmp - float(ra_h)) * 60.0
    ra_m = str(int(tmp))
    tmp = (tmp - float(ra_m)) * 6000.0
    print tmp, float(int(tmp + 0.5))
    ra_s = '%.2f' % (float(int(tmp + 0.5)) / 100)

    if dec < 0.0:
        tmp = -dec
    else:
        tmp = dec
    dec_h = str(int(tmp))
    tmp = (tmp - float(dec_h)) * 60.0
    dec_m = str(int(tmp))
    tmp = (tmp - float(dec_m)) * 600.0
    dec_s = '%.1f' % (int(tmp + 0.5) / 10)
    if dec < 0.0:
        dec_h = '-' + dec_h
    outra = ra_h + ':' + ra_m + ':' + ra_s
    outdec = dec_h + ':' + dec_m + ':' + dec_s

    return outra, outdec

# -----------------------------------------------------------
# open HDU structure

def openfits(file,mode):

    status = 0
    try:
        struct = pyfits.open(file,mode=mode)
    except:
	txt = 'ERROR -- cannot open ' + file + ' as a FITS file'
	sys.exit(txt)
    return struct, status

# -----------------------------------------------------------
# close HDU structure

def closefits(struct):

    status = 0
    try:
	struct.close()
    except:
	txt = 'ERROR -- cannot close HDU structure'
	sys.exit(txt)
    return status

# -----------------------------------------------------------
# read image from HDU structure

def readimage(struct,hdu):

    status = 0
    try:
	imagedata = struct[hdu].data
    except:
	txt = 'ERROR -- cannot read image data from HDU ' + str(hdu)
	sys.exit(txt)
    return imagedata, status

# -----------------------------------------------------------
# clear all pixels from pixel mask

def clicker1(event):

    global mask, aid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 488 and event.y < 537):
                disconnect(aid)
                disconnect(cid)
                disconnect(did)
                disconnect(eid)
                disconnect(fid)
                mask = []
                pylab.clf()
                plotimage(cmdLine)

    return

# -----------------------------------------------------------
# dump custom aperture definition file

def clicker3(event):

    global aid, cid, did, eid, fid

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 432 and event.y < 480):
                masktxt  = 'NEW|'
                masktxt += skygroup + '|'
                masktxt += '{' + re.sub('\s+',':',str(ra))
                masktxt += ',' + re.sub('\s+',':',str(dec))
                masktxt += '},TAD_NO_HALO,TAD_NO_UNDERSHOOT_COLUMN|'
                masktxt += row + '|'
                masktxt += column + '|'
                for coord in sorted(set(mask)):
                    masktxt += str(int(coord.split(',')[0]) - int(row)) + ','
                    masktxt += str(int(coord.split(',')[1]) - int(column)) + ';'
                if (os.path.isfile(maskfile)):
                    os.remove(maskfile)
                out = open(maskfile,'a')
                out.write(masktxt[:-1]+'\n')
                out.close()
                print 'Wrote custom aperture definition file ' + maskfile
    return

# -----------------------------------------------------------
# print plot to png with left-mouse click

def clicker4(event):

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 377 and event.y < 425):
                pylab.savefig(plotfile)
                print 'Wrote plot hardcopy file ' + plotfile
    return

# -----------------------------------------------------------
# close plot and exit program

def clicker5(event):

    global mask, aid, cid, did, eid, fid, done

    if event.inaxes:
        if event.button == 1:
            if (event.x > 585 and event.x < 783 and
                event.y > 320 and event.y < 368):
                disconnect(aid)
                disconnect(cid)
                disconnect(did)
                disconnect(eid)
                disconnect(fid)
    return

# -----------------------------------------------------------
# this function will be called with every click of the mouse

def clicker6(event):

    global mask, aid, cid, did, eid, fid

    if event.inaxes:
        if event.key == 'x':
            if colmap in ['Greys','binary','bone','gist_gray','gist_yarg',
                        'gray','pink','RdGy']:
                sqcol = 'g'
                alpha = 0.5
            else:
                sqcol = '#ffffee'
                alpha = 0.8
            m = float(int(event.xdata + 0.5))
            n = float(int(event.ydata + 0.5))
            txt = str(int(n))+','+str(int(m))
            if txt in mask:
                tmpmask = []
                for pixel in mask:
                    if pixel != txt:
                        tmpmask.append(pixel)
                mask = tmpmask
            else:
                mask.append(txt)
            plotimage(cmdLine)

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot(cmdLine):

    pylab.figure(1,figsize=[5,10])
    ion()
    a=outer(ones(10),arange(0,1,0.01))
    subplots_adjust(top=0.99,bottom=0.00,left=0.01,right=0.8)
    maps=[m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        subplot(l,1,i+1)
        pylab.setp(pylab.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
        imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
        pylab.text(100.85,0.5,m,fontsize=10)

# render plot

    if cmdLine: 
        pylab.show()
    else: 
        pylab.ion()
        pylab.plot([])
        pylab.ioff()
	
    status = 1
    return status

# -----------------------------------------------------------
# main
if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot sub-areas of Kepler Full Frame Images and define custom target apertures')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('ffifile', help='name of input FFI FITS file', type=str)
    parser.add_argument('--kepid', default='', help='Kepler ID of target from Kepler Input Catalog', type=str)
    parser.add_argument('--ra', default='', help='Right Ascension of target J2000 [hours or deg]', type=str)
    parser.add_argument('--dec', default='', help='Declination of target J2000 [deg]', type=str)
    parser.add_argument('--aperfile', default='', help='name of ASCII custom aperture definition file', type=str)
    parser.add_argument('--imin', default=1.5E5, help='minimum of image intensity scale [e-]', type=float)
    parser.add_argument('--imax', default=5.0E6, help='minimum of image intensity scale [e-]', type=float)
    parser.add_argument('--iscale', default='logarithmic', help='type of image intensity scale', 
                        type=str, choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', default='PuBu', help='image colormap', type=str)
    parser.add_argument('--npix', default=30, help='pixel dimension of subimage', type=float)
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', 
                        dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True

    kepffi(args.ffifile,args.kepid,args.ra,args.dec,args.aperfile,args.imin,args.imax,args.iscale,args.cmap,
           args.npix,args.verbose,args.logfile,args.status,cmdLine)
    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepffi.par")
    t = iraf.IrafTaskFactory(taskname="kepffi", value=parfile, function=kepffi)


