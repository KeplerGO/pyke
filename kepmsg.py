#!/usr/bin/env python

import sys, time, string

# -----------------------------------------------------------
# write message to log file and shell

def log(file,message,verbose):

    if (verbose):

# print to shell

	print (message)

# print to log file

	if (file):
	    output = open(file,'a')
	    output.write(message+'\n')
	    output.close()

# -----------------------------------------------------------
# write error message to log file and shell

def err(file,message,verbose):

    if verbose:
        log(file,message,verbose)
    else:
        pass
    return 1

# -----------------------------------------------------------
# write warning message to log file and shell

def warn(file,message):

    log(file,message,True)
    return 0

# -----------------------------------------------------------
# write error message and time to shell and exit

def abort(message,file,verbose):

    clock('Abort time is: ',file,verbose)
    if verbose:
        log(file,message,True)
    else:
        print message
        sys.exit(2)

# -----------------------------------------------------------
# write error message to shell and exit

def exit(message):

    sys.exit(message)

# -----------------------------------------------------------
# write time to log file and shell

def clock(text,file,verbose):

    if (verbose):
	message = text + ': ' + time.asctime(time.localtime())
	log(file,message,verbose)

# -----------------------------------------------------------
# write message to log file

def file(file,message,verbose):

    if (file and verbose):
        output = open(file,'a')
        output.write(message+'\n')
        output.close()

# -----------------------------------------------------------
# test the logfile name is good

def test(file):

    newlog = string.join(file.split(),"")
    if (len(newlog) == 0): newlog = 'kepler.log'
    if (newlog != file):
	file = newlog
	message = 'WARNING: logfile renamed to '+file+'\n\n'
	log(False,message,True)

    return file
