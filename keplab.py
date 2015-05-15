#!/usr/bin/env python

import pylab, numpy
from pylab import *
from matplotlib import *
from numpy import *

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot():

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
    ioff()
    status = 1
    return status

# -----------------------------------------------------------
# image intensity min and max

def ImageMinMax(img,frac):

    if status == 0:
        nstat = 2; pixels = []
        for i in range(img.shape[0]):
            for j in range(img.shape[1]):
                pixels.append(img[i,j])
        pixels = array(sort(pixels),dtype=float32)
        if int(float(len(pixels)) / 10 + 0.5) > nstat:
            nstat = int(float(len(pixels)) / 10 + 0.5)
        zmin = median(pixels[:nstat])
        zmax = median(pixels[-nstat:])

    return zmin, zmax
