import sys, time
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import kepio, kepmsg, kepkey

# global variables

labelsize = 24; ticksize = 16; xsize = 17; ysize = 7
lcolor = '#0000ff'; lwidth = 1.0; fcolor = '#ffff00'; falpha = 0.2
instr = ''; cadence = 1800.0; barytime0 = 0
nrm = 1; barytime = []; flux = []; xmin = 0; xmax = 1
ymin = 0; ymax = 1; xr = 1; yr = 1; xlab = ''; ylab = ''
mask = []; aid = None; bid = None; cid = None; did = None; eid = None; fid = None
clobb = True; outf = ''; verb = True; logf = ''; rinf = ''
cmdLine = False

def keprange(infile, rinfile, outfile, column, clobber, verbose, logfile,
             status, cLine=False):

# startup parameters

    status = 0
    global instr, cadence, barytime0, nrm, barytime, flux
    global xmin, xmax, ymin, ymax, xr, yr, xlab, ylab
    global clobb, outf, verb, logf, rinf, col, bjdref, cade, cmdLine

# log the call

    if rinfile.lower() == 'none':
        rinfile = ''
    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPRANGE -- '
    call += 'infile='+infile+' '
    call += 'rinfile='+rinfile+' '
    call += 'outfile='+outfile+' '
    call += 'column='+column+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)
    clobb = clobber
    outf = outfile
    verb = verbose
    logf = logfile
    rinf = rinfile
    cmdLine = cLine

# start time

    kepmsg.clock('KEPRANGE started at: ',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile):
        message = 'ERROR -- KEPRANGE: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

## open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,
                                                                infile,
                                                                logfile,
                                                                verbose,
                                                                status)
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence
    cade = cadenom

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr, infile, logfile, verbose)

# input data

    if status == 0:
        table = instr[1].data

# filter out NaNs

    work1 = []; work2 = []
    col = column
    if status == 0:
        barytime, status = kepio.readtimecol(infile, table, logfile, verbose)
    if status == 0:
        try:
            flux = instr[1].data.field(col)
        except:
            message = 'ERROR -- KEPRANGE: no column named ' + col + ' in table ' +  infile + '[1]'
            status = kepmsg.err(file,message,verbose)
    if status == 0:
        for i in range(len(barytime)):
            if (np.isfinite(barytime[i]) and np.isfinite(flux[i]) and flux[i] != 0.0):
                work1.append(barytime[i] + bjdref)
                work2.append(flux[i])
        barytime = np.array(work1, dtype=np.float64)
        flux = np.array(work2, dtype=np.float32) / cadenom

# clean up x-axis unit

    if status == 0:
        barytime0 = float(int(tstart / 100) * 100.0)
        barytime = barytime - barytime0
        xlab = 'BJD $-$ %d' % barytime0

# clean up y-axis units

    if status == 0:
        nrm = len(str(int(flux.max())))-1
        flux = flux / 10**nrm
        ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

# data limits

        xmin = barytime.min()
        xmax = barytime.max()
        ymin = flux.min()
        ymax = flux.max()
        xr = xmax - xmin
        yr = ymax - ymin
        flux[0] = 0.0
        flux[-1] = 0.0

# plot new light curve

    if status == 0:
        plt.figure(figsize=[xsize,ysize])
        plt.clf()
        plotlc(cmdLine)

    #if status == 0:
    #    status = kepio.closefits(instr, logfile, verbose)


# -----------------------------------------------------------
# plot light curve and tool

def plotlc(cmdLine):

    global aid, bid, cid, did, eid, fid, mask

# load button

    plt.ion()
    plt.clf()
    plt.axes([0.06,0.02,0.22,0.1])
    plt.text(0.5,0.5,'LOAD',fontsize=24,weight='heavy',
             horizontalalignment='center', verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    aid = plt.connect('button_press_event',clicker1)

# save button

    plt.axes([0.2933,0.02,0.22,0.1])
    plt.text(0.5,0.5,'SAVE',fontsize=24,weight='heavy',
             horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    bid = plt.connect('button_press_event',clicker2)

# clear button

    plt.axes([0.5266,0.02,0.22,0.1])
    plt.text(0.5,0.5,'CLEAR',fontsize=24,weight='heavy',
             horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    cid = plt.connect('button_press_event',clicker3)

# print button

    plt.axes([0.76,0.02,0.22,0.1])
    plt.text(0.5,0.5,'PRINT',fontsize=24,weight='heavy',
             horizontalalignment='center',verticalalignment='center')
    plt.setp(plt.gca(),xticklabels=[],xticks=[],yticklabels=[],yticks=[])
    plt.fill([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0],'#ffffee')
    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    did = plt.connect('button_press_event',clicker4)

# light curve

    plt.axes([0.06,0.213,0.92,0.77])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ltime = []; ldata = []
    work1 = instr[1].data.field(0) + bjdref
    work2 = instr[1].data.field(col) / cade
    for i in range(len(work1)):
        if np.isfinite(work1[i]) or np.isfinite(work2[i]):
            ltime.append(work1[i])
            ldata.append(work2[i])
        else:
            ltime = np.array(ltime, dtype=np.float64) - barytime0
            ldata = np.array(ldata, dtype=np.float64) / 10**nrm
            plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)
            ltime = []; ldata = []
    ltime = np.array(ltime, dtype=np.float64) - barytime0
    ldata = np.array(ldata, dtype=np.float64) / 10**nrm
    plt.plot(ltime,ldata,color=lcolor,linestyle='-',linewidth=lwidth)
    plt.fill(barytime,flux,fc=fcolor,linewidth=0.0,alpha=falpha)
    plt.xlabel(xlab, {'color' : 'k'})
    plt.ylabel(ylab, {'color' : 'k'})
    plt.grid()

# plt.plot masks

    for i in range(len(mask)):
        t = float(mask[i])
        plt.plot([t,t],[ymin,ymax],color='g',linestyle='-',linewidth=0.5)
    nt = 0
    for i in range(int(len(mask)/2)):
        t1 = float(mask[nt])
        t2 = float(mask[nt+1])
        nt += 2
        plt.fill([t1,t1,t2,t2,t1],[ymin,ymax,ymax,ymin,ymin],
                 fc='g',linewidth=0.0,alpha=0.5)
# plot ranges

    plt.xlim(xmin-xr*0.01,xmax+xr*0.01)
    if ymin-yr*0.01 <= 0.0:
        plt.ylim(1.0e-10,ymax+yr*0.01)
    else:
        plt.ylim(ymin-yr*0.01,ymax+yr*0.01)

# ranges

    eid = plt.connect('key_press_event',clicker6)

# render plot

    plt.ion()
    plt.show()

# -----------------------------------------------------------
# load mask from ascii file

def clicker1(event):

    global mask, aid, bid, cid, did, eid, fid, cmdLine

    if event.inaxes:
        if event.button == 1:
            if (event.x > 83 and event.x < 383 and
                event.y > 12 and event.y < 68):
                if kepio.fileexists(rinf):
                    mask = []
                    lines, status = kepio.openascii(rinf,'r',logf,verb)
                    for line in lines:
                        line = line.strip().split(',')
                        try:
                            float(line[0])
                            float(line[1])
                            if barytime0 > 2.4e6:
                                mask.append(float(line[0]) - barytime0)
                                mask.append(float(line[1]) - barytime0)
                            else:
                                mask.append(float(line[0]) - barytime0 - 2.4e6)
                                mask.append(float(line[1]) - barytime0 - 2.4e6)
                        except:
                            message = 'ERROR -- KEPRANGE: ascii format of ranges '
                            message += 'file not recognized.'
                            status = kepmsg.err(logf,message,False)
                    plt.disconnect(aid)
                    plt.disconnect(bid)
                    plt.disconnect(cid)
                    plt.disconnect(did)
                    plt.disconnect(eid)
                    plt.disconnect(fid)
                    plotlc(cmdLine)
                else:
                    print('WARNING -- KEPRANGE: input ranges file does not exist or was not provided')
    return

# -----------------------------------------------------------
# save mask to ascii file

def clicker2(event):

    global mask, aid, bid, cid, did, eid, fid, clobb, cmdLine

    if clobb: status = kepio.clobber(outf,logf,verb)
    if kepio.fileexists(outf):
        message = 'ERROR -- KEPRANGE: ' + outf + ' exists. Use --clobber'
        status = kepmsg.err(logf,message,verb)
    else:
        if event.inaxes:
            if event.button == 1:
                if (event.x > 402 and event.x < 702 and
                    event.y > 12 and event.y < 68):
                    nt = 0; txt = ''
                    for i in range(int(len(mask)/2)):
                        t1 = float(mask[nt]) + barytime0
                        t2 = float(mask[nt+1]) + barytime0
                        if t1 < 2.4e6: t1 += 2.4e6
                        if t2 < 2.4e6: t2 += 2.4e6
                        txt += str(t1) + ',' + str(t2) + '\n'
                        nt += 2
                    txt = txt.strip()
                    print txt
                    kepmsg.file(outf,txt,True)
                    print '\nWrote ASCII file ' + outf
                    plotlc(cmdLine)
    return

# -----------------------------------------------------------
# clear time domain mask

def clicker3(event):

    global mask, aid, bid, cid, did, eid, fid, cmdLine

    if event.inaxes:
        if event.button == 1:
            if (event.x > 723 and event.x < 1022 and
                event.y > 12 and event.y < 68):
                mask = []
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                plotlc(cmdLine)
    return

# -----------------------------------------------------------
# print time domain mask

def clicker4(event):

    global mask, aid, bid, cid, did, eid, fid, cmdLine

    if event.inaxes:
        if event.button == 1:
            if (event.x > 1042 and event.x < 1342 and
                event.y > 12 and event.y < 68):
                nt = 0; txt = ''
                for i in range(int(len(mask)/2)):
                    t1 = float(mask[nt]) + barytime0
                    t2 = float(mask[nt+1]) + barytime0
                    if t1 < 2.4e6: t1 += 2.4e6
                    if t2 < 2.4e6: t2 += 2.4e6
                    txt += str(t1) + ',' + str(t2) + '\n'
                    nt += 2
                txt = txt.strip()
                print('\n' + txt)
                plt.disconnect(aid)
                plt.disconnect(bid)
                plt.disconnect(cid)
                plt.disconnect(did)
                plt.disconnect(eid)
                plt.disconnect(fid)
                plotlc(cmdLine)
    return

# -----------------------------------------------------------
# right-click create time ranges

def clicker6(event):

    global mask, aid, bid, cid, did, eid, fid, cmdLine

    if event.inaxes:
#        event.key
        if event.key == '' or \
        event.key == '1' or \
        event.key == '2' or \
        event.key == '3' or \
        event.key == '4' or \
        event.key == '5' or \
        event.key == '6' or \
        event.key == '7' or \
        event.key == '8' or \
        event.key == '9' or \
        event.key == '0' or \
        event.key == '-' or \
        event.key == '=' or \
        event.key == 'q' or \
        event.key == 'w' or \
        event.key == 'e' or \
        event.key == 'r' or \
        event.key == 't' or \
        event.key == 'y' or \
        event.key == 'u' or \
        event.key == 'i' or \
        event.key == 'o' or \
        event.key == 'p' or \
        event.key == '[' or \
        event.key == ']' or \
        event.key == 'a' or \
        event.key == 'd' or \
        event.key == 'f' or \
        event.key == 'g' or \
        event.key == 'h' or \
        event.key == 'j' or \
        event.key == 'k' or \
        event.key == 'l' or \
        event.key == ';' or \
        event.key == 'x' or \
        event.key == 'c' or \
        event.key == 'v' or \
        event.key == 'b' or \
        event.key == 'n' or \
        event.key == 'm' or \
        event.key == ',' or \
        event.key == '.' or \
        event.key == '/':
            if (event.x > 83 and event.x < 1337 and
                event.y > 122 and event.y < 558):
                if len(mask) % 2 == 0:
                    mask.append(event.xdata)
                else:
                    if event.xdata > mask[-1]:
                        mask.append(event.xdata)
                    else:
                        mask.append(mask[-1])
                        mask[-2] = event.xdata
                plotlc(cmdLine)
        if event.key == 'z':
            kepio.closefits(instr, "keprangelog.txt", True)
            plt.close()
            return

# main
if '--shell' in sys.argv:
    import argparse
    parser = argparse.ArgumentParser(description='Interactively define and store time ranges via a GUI')
    parser.add_argument('--shell', action='store_true',
                        help='Are we running from the shell?')
    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('--rinfile', default='',
                        help='Name of input ASCII time ranges file', type=str)
    parser.add_argument('--outfile', default='',
                        help='Name of output ASCII time ranges file',
                        type=str)
    parser.add_argument('--column', default='SAP_FLUX',
                        help='Name of diagnostic FITS column', type=str)
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true',
                        help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file',
                        default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)',
                        default=0, dest='status', type=int)
    args = parser.parse_args()
    cmdLine=True
    keprange(args.infile, args.rinfile, args.outfile, args.column,
             args.clobber, args.verbose, args.logfile, args.status, cmdLine)
else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$keprange.par")
    t = iraf.IrafTaskFactory(taskname="keprange", value=parfile, function=keprange)
