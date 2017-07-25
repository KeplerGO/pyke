import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepbls import kepbls
from ..kepio import delete

# fake light curve with transits of period = 2.02 days
# and duration 0.18 days
fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_kepbls():
    kepbls(fake_lc, "kepbls.fits", datacol='SAP_FLUX', errcol='SAP_FLUX_ERR',
           minper=1, maxper=3, mindur=0, maxdur=.5, nsearch=100, nbins=100,
           overwrite=True)
    f = pyfits.open("kepbls.fits")
    best_guess = f[3].data['SIG_RES'] == max(f[3].data['SIG_RES'])
    assert abs(f[3].data['PERIOD'][best_guess] - 2.02) < 0.1
    assert abs(f[3].data['DURATION'][best_guess] - 0.18) < 0.01
    delete("kepbls.fits", "log_kepextract.txt", False)
