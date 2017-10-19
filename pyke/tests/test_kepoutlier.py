import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_array_almost_equal
from ..kepoutlier import kepoutlier
from ..kepio import delete

clean_fake_lc = get_pkg_data_filename("data/golden-lc.fits")
nan_fake_lc = get_pkg_data_filename("data/golden-lc-with-nans.fits")
sigma_clipped_lc = get_pkg_data_filename("data/golden-lc-kepoutlier.fits")
nan_sigma_clipped_lc = get_pkg_data_filename("data/golden-lc-with-nans-kepoutlier.fits")

@pytest.mark.parametrize("fake_lc, true_lc",
                         [(clean_fake_lc, sigma_clipped_lc),
                          (nan_fake_lc, nan_sigma_clipped_lc)])
def test_kepoutlier(fake_lc, true_lc):
    kepoutlier(fake_lc, outfile="kepoutlier.fits", datacol="SAP_FLUX",
            nsig=2.0, stepsize=5, overwrite=True)
    f = pyfits.open("kepoutlier.fits")
    g = pyfits.open(true_lc)
    assert_array_almost_equal(f[1].data['SAP_FLUX'],
                              g[1].data['SAP_FLUX'])
    assert_array_almost_equal(f[1].data['TIME'],
                              g[1].data['TIME'])
    f.close()
    g.close()
    delete("kepoutlier.fits", "log_kepoutlier.txt", False)
