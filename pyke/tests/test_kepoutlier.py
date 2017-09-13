import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_array_almost_equal
from ..kepoutlier import kepoutlier
from ..kepio import delete

fake_lc = get_pkg_data_filename("data/golden-lc.fits")
sigma_clipped_lc = get_pkg_data_filename("data/golden-lc-kepoutlier.fits")

def test_kepoutlier():
    kepoutlier(fake_lc, outfile="kepoutlier.fits", datacol="PDCSAP_FLUX",
            nsig=2.0, stepsize=5, overwrite=True)
    f = pyfits.open("kepoutlier.fits")
    g = pyfits.open(sigma_clipped_lc)
    assert_array_almost_equal(f[1].data['PDCSAP_FLUX'],
                              g[1].data['PDCSAP_FLUX'])
    assert_array_almost_equal(f[1].data['TIME'],
                              g[1].data['TIME'])
    f.close()
    g.close()
    delete("kepoutlier.fits", "log_kepoutlier.txt", False)
