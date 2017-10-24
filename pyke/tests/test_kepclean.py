from..kepclean import kepclean
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits as pyfits
import numpy as np
import os
from ..kepio import delete

LC_filename  =  get_pkg_data_filename("data/golden-lc-with-nans.fits")
TPF_filename  =  get_pkg_data_filename("data/test-tpf-with-nans.fits")

def test_kepclean():
    kepclean(LC_filename, "kepclean.fits", overwrite=True)
    h = pyfits.open('kepclean.fits')
    j = pyfits.open(LC_filename)
    #Check it creates a NEW extention
    assert len(j) == len(h)
    #Check it creates a WINDOW FUNCTION extention
    assert h[1].header['NANCLEAN'] == True
    #Check that there are the correct number of points in the POWER field
    assert len(h[1].data) < len(j[1].data)
    #Check length of cleaned file
    assert len(h[1].data) == 1087
    assert len(j[1].data) == 1287
    h.close()
    j.close()
    delete("kepclean.fits", "kepclean.log", False)

def test_kepclean_tpf():
    #This file doesn't have any frames which are ALL nans so it shouldn't change
    kepclean(TPF_filename, "kepclean.fits", overwrite=True)
    h = pyfits.open('kepclean.fits')
    j = pyfits.open(TPF_filename)
    assert len(j) == len(h)
    assert h[1].header['NANCLEAN'] == True
    assert len(h[1].data) == len(j[1].data)
    h.close()
    j.close()
    delete("kepclean.fits", "kepclean.log", False)
