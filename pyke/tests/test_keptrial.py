import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..keptrial import keptrial
from ..kepio import delete

# fake light curve with transits of period = 2.02 days
# and duration 0.18 days
fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_keptrial():
    keptrial(fake_lc, "keptrial.fits", datacol='SAP_FLUX', errcol='SAP_FLUX_ERR',
             fmin=0.4, fmax=0.67, ntrials=50, nfreq=10, overwrite=True)
    f = pyfits.open("keptrial.fits")
    assert abs(f[3].header['PERIOD'] - 2.02) < 1e-3
    delete("keptrial.fits", "log_keptrial.txt", False)
