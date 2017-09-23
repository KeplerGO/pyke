import pytest
from astropy.io import fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from ..kepextract import kepextract
from ..kepio import delete

# 3 x 3 target pixel file filled with zeros everywhere
tpf_all_zeros = get_pkg_data_filename("data/test-tpf-all-zeros.fits")
# 3 x 3 tpf with ones in the center pixel and zeros everywhere else
tpf_one_center = get_pkg_data_filename("data/test-tpf-non-zero-center.fits")
# 3 x 3 tpf as follows:
# 0 1 0
# 1 1 1
# 0 1 0
tpf_star = get_pkg_data_filename("data/test-tpf-star.fits")
# 3 x 3 tpf as follows:
# 0   1   0
# 1  nan  1
# 0   1   0
tpf_nans = get_pkg_data_filename("data/test-tpf-with-nans.fits")
# mask file selecting pixel [1, 1]
maskfile = get_pkg_data_filename("data/center-mask.txt")

@pytest.mark.parametrize("tpf, maskfile, bkg, answer",
                         [(tpf_all_zeros, 'ALL', False, 0),
                          (tpf_all_zeros, 'ALL', True, 0),
                          (tpf_all_zeros, maskfile, False, 0),
                          (tpf_all_zeros, maskfile, True, 0),
                          (tpf_one_center, 'ALL', False, 1),
                          (tpf_one_center, 'ALL', True, 1),
                          (tpf_one_center, maskfile, False, 1),
                          (tpf_one_center, maskfile, True, 1),
                          (tpf_star, 'ALL', False, 5),
                          (tpf_star, 'ALL', True, -4),
                          (tpf_star, maskfile, False, 1),
                          (tpf_star, maskfile, True, 0),
                          (tpf_nans, 'ALL', False, 4),
                          (tpf_nans, 'ALL', True, 0)])
def test_kepextract(tpf, maskfile, bkg, answer):
    kepextract(tpf, outfile="lc.fits", maskfile=maskfile, bkg=bkg, overwrite=True)
    f = pyfits.open("lc.fits")
    assert (f[1].data['SAP_FLUX'] == answer).all()
    f.close()
    delete("lc.fits", "log_kepextract.txt", False)
