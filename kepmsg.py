import sys, time, string

def log(filename, message, verbose):
    """write message to log file and shell."""
    if verbose:
        print(message)
        if filename:
            output = open(filename, 'a')
            output.write(message + '\n')
            output.close()

def err(filename, message, verbose):
    if verbose:
        log(filename, message, verbose)

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
        print(message)
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
