#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
import matplotlib
matplotlib.use('TkAgg')

from .version import __version__
from .keparray import *
from .kepbls import *
from .kepclean import *
from .kepclip import *
from .kepconvert import *
from .kepcotrend import *
from .kepdetrend import *
from .kepdiffim import *
from .kepdraw import *
from .kepdynamic import *
from .kepextract import *
from .kepfilter import *
from .kepfit import *
from .kepflatten import *
from .kepfold import *
from .kepfourier import *
from .kepperiodogram import *
from .kepfunc import *
from .kephead import *
from .kepimages import *
from .kepio  import *
from .kepkey import *
from .kepmask import *
from .kepmsg import *
from .kepoutlier import *
from .keppca import *
from .keppixseries import *
from .kepplot import *
from .kepprf import *
from .kepprfphot import *
from .keprange import *
from .kepsff import *
from .kepsmooth import *
from .kepstat import *
from .kepstddev import *
from .kepstitch import *
from .keptimefix import *
from .keptrial import *
from .keptrim import *
from .kepwindow import *
from .prf import *
from .lightcurve import *
from .targetpixelfile import *
from .utils import *
