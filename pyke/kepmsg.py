import sys
import time
import logging

def log(filename, message, verbose):
    """write message to log file and shell."""
    if verbose:
        print(message)
    if filename:
        output = open(filename, 'a')
        output.write(message + '\n')
        output.close()

def err(filename, message, verbose):
    log(filename, message, verbose)
    raise Exception(message)

def warn(filename, message, verbose):
    log(filename, message, verbose)

def exit(message):
    """write error message to shell and exit"""
    sys.exit(message)

def clock(text, filename, verbose):
    """write time to log file and shell"""
    if verbose:
        message = text + ': ' + time.asctime(time.localtime())
        log(filename, message, verbose)
