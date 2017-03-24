import numpy, sys, time
from astropy.io import fits as pyfits
import kepio, kepmsg, kepkey, kepstat

def kepstitch(infiles,outfile,clobber,verbose,logfile,status): 

# startup parameters

    status = 0
    lct = []; bjd = []

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPSTITCH -- '
    call += 'infiles='+infiles+' '
    call += 'outfile='+outfile+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPSTITCH started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# parse input file list

    infiles, status = kepio.parselist(infiles,logfile,verbose)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPSTITCH: ' + outfile + ' exists. Use clobber=yes'
	    kepmsg.err(logfile,message,verbose)
	    status = 1

# open output file

    if status == 0:
	    outstr, status = kepio.openfits(infiles[0],'readonly',logfile,verbose)
	    nrows1 = outstr[1].data.shape[0]

# fudge non-compliant FITS keywords with no values

    if status == 0:
	    outstr = kepkey.emptykeys(outstr,file,logfile,verbose)
	    head0 = outstr[0].header
	    head1 = outstr[1].header

# open input files

    nfiles = 0
    if status == 0:
	    for infile in infiles: 
		    instr, status = kepio.openfits(infile,'readonly',logfile,verbose)

# append table data

		    if nfiles > 0:
			    nrows2 = instr[1].data.shape[0]
			    nrows = nrows1 + nrows2
			    outtab = pyfits.new_table(outstr[1].columns,nrows=nrows)
			    for name in outstr[1].columns.names:
                                try:
				    outtab.data.field(name)[nrows1:]=instr[1].data.field(name)
                                except:
                                    message = 'ERROR -- KEPSTITCH: column ' + name + ' missing from some files.'
                                    kepmsg.warn(logfile,message)
                                    pass
			    outstr[1] = outtab
			    outstr[0].header = head0
			    outstr[1].header = head1
			    nrows1 = nrows

# start and stop times of data

                    fitsvers = 1.0
		    lc_start, status = kepkey.get(infile,instr[1],'LC_START',logfile,verbose)
		    lc_end, status = kepkey.get(infile,instr[1],'LC_END',logfile,verbose)
                    try:
                        startbjd = instr[1].header['STARTBJD']
                    except:
                        startbjd, status = kepkey.get(infile,instr[1],'TSTART',logfile,verbose)
                        fitsvers = 2.0                        
                    try:
                        endbjd = instr[1].header['ENDBJD']
                    except:
                        endbjd, status = kepkey.get(infile,instr[1],'TSTOP',logfile,verbose)
                        fitsvers = 2.0
		    lct.append(lc_start); lct.append(lc_end)
		    bjd.append(startbjd); bjd.append(endbjd)

# close input files

		    status = kepio.closefits(instr,logfile,verbose)
		    nfiles += 1

# maxmimum and minimum times in file sample

    if status == 0:
	    lc_start = kepstat.min(lct)
	    lc_end = kepstat.max(lct)
	    startbjd = kepstat.min(bjd)
	    endbjd = kepstat.max(bjd)
	    status = kepkey.change('LC_START',lc_start,outstr[1],outfile,logfile,verbose)
	    status = kepkey.change('LC_END',lc_end,outstr[1],outfile,logfile,verbose)
            if fitsvers == 1.0:
                status = kepkey.change('STARTBJD',startbjd,outstr[1],outfile,logfile,verbose)
                status = kepkey.change('ENDBJD',endbjd,outstr[1],outfile,logfile,verbose)
            else:
                status = kepkey.change('TSTART',startbjd,outstr[1],outfile,logfile,verbose)
                status = kepkey.change('TSTOP',endbjd,outstr[1],outfile,logfile,verbose)                

# comment keyword in output file

    if status == 0:
	    status = kepkey.comment(call,outstr[0],outfile,logfile,verbose)

# close output file

    if status == 0:
	    outstr.writeto(outfile)
	    status = kepio.closefits(outstr,logfile,verbose)

## end time

    if (status == 0):
	    message = 'KEPSTITCH completed at'
    else:
	    message = '\nKEPSTITCH aborted at'
    kepmsg.clock(message,logfile,verbose)

# main
if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Append multiple month short cadence and/or multiple quarter long cadence data')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')
    parser.add_argument('infiles', help='List of input files', type=str)

    parser.add_argument('outfile', help='Name of FITS file to output', type=str)


    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    

    kepstitch(args.infiles,args.outfile,args.clobber,args.verbose,args.logfile,args.status)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepstitch.par")
    t = iraf.IrafTaskFactory(taskname="kepstitch", value=parfile, function=kepstitch)
