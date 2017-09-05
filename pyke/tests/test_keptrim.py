import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..keptrim import keptrim
from ..kepio import delete

# 3 x 3 tpf with ones in the center pixel and zeros everywhere else
tpf_one_center = get_pkg_data_filename("data/test-tpf-non-zero-center.fits")

def test_keptrim():
    keptrim(tpf_one_center, 1013, 918, 1, outfile='tpf.fits', overwrite=True)
    f = pyfits.open('tpf.fits')
    assert (f[1].data['FLUX'] == 1).all()
    assert f[1].data['FLUX'].shape == (len(f[1].data['FLUX']), 1, 1)
    f.close()
    delete("tpf.fits", "log_keptrim.txt", False)
