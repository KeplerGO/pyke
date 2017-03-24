# -----------------------------------------------------------
# core code

def kepsimulate(infile,plotfile,kepid,ra,dec,imscale,colmap,clobber,verbose,logfile,status,cmdLine=False): 

# input arguments

    status = 0
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFIELD -- '
    call += 'kepid='+str(kepid)+' '
    call += 'ra='+str(ra)+' '
    call += 'dec='+str(dec)+' '
    call += 'outfile='+outfile+' '
    call += 'plotfile='+plotfile+' '
    call += 'imscale='+imscale+' '
    call += 'colmap='+colmap+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPSIMULATE started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# -----------------------------------------------------------
# main

if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Fitting PRF model to Target Pixel image')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
ype=str)
    parser.add_argument('--kepid', '-k', default='', help='Kepler ID of target from Kepler Input Catalog', type=str)
    parser.add_argument('--ra', '-r', default='', help='Right Ascension of target J2000 [hours or deg]', type=str)
    parser.add_argument('--dec', '-d', default='', help='Declination of target J2000 [deg]', type=str)
    parser.add_argument('--outfile', '-o', help='Name of output target pixel file', default='simulation.fits', dest='outfile', type=str)
    parser.add_argument('--plotfile', '-p', help='Name of output PNG plot file', default='None', dest='plotfile', type=str)
    parser.add_argument('--imscale', '-i', help='Type of image intensity scale', default='linear', dest='imscale', type=str,choices=['linear','logarithmic','squareroot'])
    parser.add_argument('--cmap', '-c', help='Image colormap', default='YlOrBr', dest='cmap', type=str,choices=['Accent','Blues','BrBG','BuGn','BuPu','Dark2','GnBu','Greens','Greys','OrRd','Oranges','PRGn','Paired','Pastel1','Pastel2','PiYG','PuBu','PuBuGn','PuOr','PuRd','Purples','RdBu','RdGy','RdPu','RdYlBu','RdYlGn','Reds','Set1','Set2','Set3','Spectral','YlGn','YlGnBu','YlOrBr','YlOrRd','afmhot','autumn','binary','bone','brg','bwr','cool','copper','flag','gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_yarg','gnuplot','gnuplot2','gray','hot','hsv','jet','ocean','pink','prism','rainbow','seismic','spectral','spring','summer','terrain','winter','browse'])
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', default='kepsimulatephot.log', help='Name of ascii log file', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)

    args = parser.parse_args()
    cmdLine=True
    kepsimulate(args.kepid,args.ra,args.dec,args.outfile,args.plotfile,args.imscale,args.cmap,
             args.verbose,args.logfile,args.status,cmdLine)
    
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepsimulate.par")
    t = iraf.IrafTaskFactory(taskname="kepsimulate", value=parfile, function=kepsimulate)
