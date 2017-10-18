import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepwindow import kepwindow
from ..kepdynamic import kepdynamic
from ..kepperiodogram import kepperiodogram
from ..kepfold import kepfold
from ..kepio import delete
import numpy as np

# fake light curve with transits of period = 2.02 days
# and duration 0.18 days
fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_kepwindow():
    kepwindow(fake_lc, "kepwindow.fits", overwrite=True, noninteractive=True)
    h = pyfits.open('kepwindow.fits')
    #Check it creates a NEW extention
    assert len(pyfits.open(fake_lc)) == len(h) - 1
    #Check it creates a WINDOW FUNCTION extention
    assert h[-1].header['EXTNAME'] == 'WINDOW FUNCTION'
    #Check that there are the correct number of points in the POWER field
    assert len(h[-1].data['POWER']) == 66
    h.close()
    delete("kepwindow.fits", "kepwindow.log", False)

def test_kepdynamic():
    kepdynamic(fake_lc, "kepdynamic.fits", overwrite=True, noninteractive=True)
    h = pyfits.open('kepdynamic.fits')
    #Check it creates a NEW extention
    assert len(pyfits.open(fake_lc)) == len(h) - 1
    #Check it creates a WINDOW FUNCTION extention
    assert h[-1].header['EXTNAME'] == 'DYNAMIC FT'
    #Check that the data is the correct shape
    assert np.shape(h[-1].data) == (2000,10)
    h.close()
    delete("kepdynamic.fits", "kepdynamic.log", False)

def test_kepperiodogram():
    kepperiodogram(fake_lc, "kepperiodogram.fits", overwrite=True, noninteractive=True)
    h = pyfits.open('kepperiodogram.fits')
    #Check it creates a NEW extention
    assert len(pyfits.open(fake_lc)) == len(h) - 1
    #Check it creates a WINDOW FUNCTION extention
    assert h[-1].header['EXTNAME'] == 'POWER SPECTRUM'
    #Check it finds the right period
    assert np.isclose(h[-1].header['PERIOD'],2.056372801152145,rtol=1E-5)
    h.close()

def test_kepfold():
    kepfold("kepperiodogram.fits", "kepfold.fits", bindata=True,
    overwrite=True, noninteractive=True)
    h = pyfits.open('kepfold.fits')
    #Check it creates a NEW extention
    assert len(pyfits.open("kepperiodogram.fits")) == len(h) - 2
    #Check it creates two FOLDED extentions
    assert h[-2].header['EXTNAME'] == 'FOLDED LC'
    assert h[-1].header['EXTNAME'] == 'BINNED FOLDED LC'
    #Check there is the right amount of binned data
    assert len(h[-1].data['FLUX_BINNED']) == 1000
    delete("kepperiodogram.fits", "kepperiodogram.log", False)
    delete("kepfold.fits", "kepfold.log", False)
