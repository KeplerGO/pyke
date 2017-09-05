import pytest
import numpy as np
from numpy.testing import assert_array_equal
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

@pytest.mark.parametrize("tpf, maskfile, bkg, qual, answer",
                         [(tpf_all_zeros, 'ALL', False, False, 0),
                          (tpf_all_zeros, 'ALL', False, True,  0),
                          (tpf_all_zeros, 'ALL', True,  False, 0),
                          (tpf_all_zeros, 'ALL', True,  True,  0),
                          (tpf_all_zeros, maskfile, False, False, 0),
                          (tpf_all_zeros, maskfile, False, True, 0),
                          (tpf_all_zeros, maskfile, True, False, 0),
                          (tpf_all_zeros, maskfile, True, True, 0),
                          (tpf_one_center, 'ALL', False, False, 1),
                          (tpf_one_center, 'ALL', False, True, 1),
                          (tpf_one_center, 'ALL', True, False, 1),
                          (tpf_one_center, 'ALL', True, True, 1),
                          (tpf_one_center, maskfile, False, False, 1),
                          (tpf_one_center, maskfile, False, True, 1),
                          (tpf_one_center, maskfile, True, False, 1),
                          (tpf_one_center, maskfile, True, True, 1),
                          (tpf_star, 'ALL', False, False, 5),
                          (tpf_star, 'ALL', False, True, 5),
                          (tpf_star, 'ALL', True, False, -4),
                          (tpf_star, 'ALL', True, True, -4),
                          (tpf_star, maskfile, False, False, 1),
                          (tpf_star, maskfile, False, True, 1),
                          (tpf_star, maskfile, True, False, 0),
                          (tpf_star, maskfile, True, True, 0),
                          (tpf_nans, 'ALL', False, False, 4),
                          (tpf_nans, 'ALL', False, True, 4),
                          (tpf_nans, 'ALL', True, False, 0),
                          (tpf_nans, 'ALL', True, True, 0)])
def test_kepextract(tpf, maskfile, bkg, qual, answer):
    xtpf = pyfits.open(tpf)
    time_ans = xtpf[1].data['TIME']
    mask_time_ans = np.isnan(time_ans)
    quality = xtpf[1].data['QUALITY']
    full_mask = (quality == 0) & (~mask_time_ans)

    kepextract(tpf, outfile="lc.fits", qual=qual, maskfile=maskfile, bkg=bkg,
               overwrite=True)
    f = pyfits.open("lc.fits")
    assert (f[1].data['SAP_FLUX'] == answer).all()

    if qual:
        len_ans = full_mask.sum()
        assert len(f[1].data['TIME']) == len_ans
        assert_array_equal(np.isnan(f[1].data['TIME']),
                           np.zeros(len_ans, dtype=bool))
        assert_array_equal(f[1].data['TIME'], time_ans[full_mask])
    else:
        mask_time_des = np.isnan(f[1].data['TIME'])
        assert_array_equal(mask_time_ans, mask_time_des)
        assert_array_equal(time_ans[~mask_time_ans], f[1].data['TIME'][~mask_time_des])

    f.close()
    delete("lc.fits", "log_kepextract.txt", False)
