from numpy.testing import assert_almost_equal
from ..lightcurve import KeplerCBVCorrector, KeplerLightCurveFile


def test_kepler_cbv_fit():
    tabby_quarter_8 = ("https://archive.stsci.edu/missions/kepler/lightcurves"
                       "/0084/008462852/kplr008462852-2011073133259_llc.fits")
    cbv = KeplerCBVCorrector(tabby_quarter_8)
    cbv_lc = cbv.correct()
    assert_almost_equal(cbv.coeffs, [0.08550738, 0.10808083])

    lcf = KeplerLightCurveFile(tabby_quarter_8)
    cbv_lcf = lcf.compute_cotrended_lightcurve()
    assert_almost_equal(cbv_lc.flux, cbv_lcf.flux)
