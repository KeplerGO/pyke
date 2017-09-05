import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepstitch import kepstitch
from ..kepio import delete

# fake light curve with transits of period = 2.02 days
# and duration 0.18 days
fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_kepbls():
    lcs = [fake_lc, fake_lc]
    kepstitch(lcs, outfile="kepstitch.fits", overwrite=True)
    f = pyfits.open("kepstitch.fits")
    g = pyfits.open(fake_lc)
    g_len = len(g[1].data['PDCSAP_FLUX'])
    assert (f[1].data['PDCSAP_FLUX'][:g_len] == g[1].data['PDCSAP_FLUX']).all() == True
    assert (f[1].data['PDCSAP_FLUX'][g_len:] == g[1].data['PDCSAP_FLUX']).all() == True
    f.close()
    g.close()
    delete("kepstitch.fits", "log_kepstitch.txt", False)
