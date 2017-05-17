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

def warn(filename, message):
    log(filename, message, True)
    return 0

def abort(message, filename, verbose):
    """write error message and time to shell and exit"""
    clock('Abort time is: ', filename, verbose)
    if verbose:
        log(filename, message, True)
    else:
        print(message)
        sys.exit(2)

def exit(message):
    """write error message to shell and exit"""
    sys.exit(message)

def clock(text, filename, verbose):
    """write time to log file and shell"""
    if (verbose):
        message = text + ': ' + time.asctime(time.localtime())
        log(filename, message,verbose)
