import pytest
import numpy as np
from numpy.testing import assert_array_equal
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepsmooth import kepsmooth
from ..kepio import delete

# fake light curve with transits of period = 2.02 days
# and duration 0.18 days
fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_kepsmooth():
    kepsmooth(fake_lc, outfile="kepsmooth.fits", fscale=0.05, overwrite=True)
    f = pyfits.open("kepsmooth.fits")
    g = pyfits.open(fake_lc)
    assert_array_equal(f[1].data['SAP_FLUX'], g[1].data['SAP_FLUX'])
    f.close()
    g.close()
    delete("kepsmooth.fits", "log_kepsmooth.txt", False)
