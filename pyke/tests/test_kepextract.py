from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepextract import kepextract
from ..kepio import delete

# 3 x 3 target pixel file filled with zeros everywhere
tpf_all_zeros = get_pkg_data_filename("data/test-tpf-all-zeros.fits")
# 3 x 3 tpf with ones in the center pixel and zeros everywhere else
tpf_one_center = get_pkg_data_filename("data/test-tpf-non-zero-center.fits")
# mask file selecting the [1, 1] pixel
maskfile = get_pkg_data_filename("data/center-mask.txt")

def test_kepextract_all_zeros():
    kepextract(tpf_all_zeros, "all-zeros-lc.fits")
    f = pyfits.open("all-zeros-lc.fits")
    assert f[1].data['SAP_FLUX'].all() == 0
    delete("all-zeros-lc.fits", "test_kepextract_all_zeros.txt", False)

def test_kepextract_all_zeros_with_mask():
    kepextract(tpf_all_zeros, "all-zeros-lc.fits", maskfile=maskfile)
    f = pyfits.open("all-zeros-lc.fits")
    assert f[1].data['SAP_FLUX'].all() == 0
    delete("all-zeros-lc.fits", "test_kepextract_all_zeros_with_mask.txt", False)

def test_kepextract_one_center():
    kepextract(tpf_one_center, "center-ones-lc.fits")
    f = pyfits.open("center-ones-lc.fits")
    assert f[1].data['SAP_FLUX'].all() == 1
    delete("center-ones-lc.fits", "test_kepextract_one_center.txt", False)

def test_kepextract_one_center_with_mask():
    kepextract(tpf_one_center, "center-ones-lc.fits", maskfile=maskfile)
    f = pyfits.open("center-ones-lc.fits")
    assert f[1].data['SAP_FLUX'].all() == 1
    delete("center-ones-lc.fits", "test_kepextract_one_center_with_mask.txt", False)
