import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_array_almost_equal
from ..kepclip import kepclip
from ..kepio import delete

fake_lc = get_pkg_data_filename("data/golden-lc.fits")
clipped_lc = get_pkg_data_filename("data/golden-lc-clipped.fits")

def test_kepclip():
    kepclip(fake_lc, "2457508.314460751,2457519.5668576737",
            outfile="kepclip.fits", overwrite=True)
    f = pyfits.open("kepclip.fits")
    g = pyfits.open(clipped_lc)
    assert_array_almost_equal(f[1].data['PDCSAP_FLUX'],
                              g[1].data['PDCSAP_FLUX'])
    assert_array_almost_equal(f[1].data['TIME'],
                              g[1].data['TIME'])
    delete("kepclip.fits", "log_kepclip.txt", False)
