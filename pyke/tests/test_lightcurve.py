from numpy.testing import assert_almost_equal
from ..lightcurve import KeplerCBVCorrector, KeplerLightCurveFile

# 8th Quarter of Tabby's star
TABBY_Q8 = ("https://archive.stsci.edu/missions/kepler/lightcurves"
            "/0084/008462852/kplr008462852-2011073133259_llc.fits")

def test_kepler_cbv_fit():
    # comparing that the two methods to do cbv fit are the nearly the same
    cbv = KeplerCBVCorrector(TABBY_Q8)
    cbv_lc = cbv.correct()
    assert_almost_equal(cbv.coeffs, [0.08534423, 0.10814261], decimal=4)

    lcf = KeplerLightCurveFile(TABBY_Q8)
    cbv_lcf = lcf.compute_cotrended_lightcurve()
    assert_almost_equal(cbv_lc.flux, cbv_lcf.flux)


def test_KeplerLightCurve():
    lcf = KeplerLightCurveFile(TABBY_Q8)
    kplc = lcf.get_lightcurve('SAP_FLUX')

    assert kplc.channel == lcf.channel
    assert kplc.campaign is None
    assert kplc.quarter == lcf.quarter
    assert kplc.mission == 'Kepler'
