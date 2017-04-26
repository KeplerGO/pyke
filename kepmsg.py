#!/usr/bin/env python

import string
import sys
import time


def log(filename, message, verbose):
    """write message to log file and shell"""
    # print to shell
    if (verbose):
        print(message)

    # print to log file
    if (filename):
        output = open(filename, 'a')
        output.write('{}\n'.format(message))
        output.close()


def err(filename, message, verbose):
    """write error message to log file and shell"""
    if verbose:
        log(filename, message, verbose)
    else:
        pass
    return 1


def warn(filename, message):
    """write warning message to log file and shell"""
    log(filename, message, True)
    return 0


def abort(message, filename, verbose):
    """write error message and time to shell and exit"""
    clock('Abort time is: ', filename, verbose)
    if verbose:
        log(filename, message, True)
    else:
        print message
        sys.exit(2)


def exit(message):
    """write error message to shell and exit"""
    sys.exit(message)


def clock(text, filename, verbose):
    """write time to log file and shell"""
    if (verbose):
        message = text + ': ' + time.asctime(time.localtime())
    log(filename, message, verbose)


def file(filename, message, verbose):
    """write message to log file"""
    if (filename and verbose):
        output = open(filename, 'a')
        output.write('{}\n'.format(message))
        output.close()


def test(filename):
    """test the logfile name is good"""
    newlog = string.join(filename.split(), "")
    if (len(newlog) == 0):
        newlog = 'kepler.log'
    if (newlog != filename):
        filename = newlog
    message = 'WARNING: logfile renamed to {}\n\n'.format(filename)
    log(False, message, True)

    return filename
