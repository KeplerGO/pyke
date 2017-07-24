import pytest
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

@pytest.mark.parametrize("tpf, maskfile, bkg, answer",
                         [(tpf_all_zeros, 'ALL', False, 0),
                          (tpf_all_zeros, 'ALL', True, 0),
                          (tpf_all_zeros, maskfile, False, 0),
                          (tpf_all_zeros, maskfile, True, 0),
                          (tpf_one_center, 'ALL', False, 1),
                          (tpf_one_center, 'ALL', True, 1),
                          (tpf_one_center, maskfile, False, 1),
                          (tpf_one_center, maskfile, True, 1)])
def test_kepextract(tpf, maskfile, bkg, answer):
    kepextract(tpf, "lc.fits", maskfile=maskfile, bkg=bkg, overwrite=True)
    f = pyfits.open("lc.fits")
    assert f[1].data['SAP_FLUX'].all() == answer
    delete("lc.fits", "log_kepextract.txt", False)
