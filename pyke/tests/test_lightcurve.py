import pytest
import numpy as np
from numpy.testing import assert_almost_equal
from ..lightcurve import LightCurve, KeplerCBVCorrector, KeplerLightCurveFile

# 8th Quarter of Tabby's star
TABBY_Q8 = ("https://archive.stsci.edu/missions/kepler/lightcurves"
            "/0084/008462852/kplr008462852-2011073133259_llc.fits")


def test_kepler_cbv_fit():
    # comparing that the two methods to do cbv fit are the nearly the same
    cbv = KeplerCBVCorrector(TABBY_Q8)
    cbv_lc = cbv.correct()
    assert_almost_equal(cbv.coeffs, [0.08534423, 0.10814261], decimal=3)
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


@pytest.mark.parametrize("quality_bitmask,answer", [('hardest', 2661),
    ('hard', 2706), ('default', 2917), (None, 3279),
    (1, 3279), (100, 3252), (2096639, 2661)])
def test_bitmasking(quality_bitmask, answer):
    '''Test whether the bitmasking behaves like it should'''
    lcf = KeplerLightCurveFile(TABBY_Q8, quality_bitmask=quality_bitmask)
    flux = lcf.get_lightcurve('SAP_FLUX').flux
    assert len(flux) == answer


def test_lightcurve_fold():
    """Test the ``LightCurve.fold()`` method."""
    lc = LightCurve(time=[1, 2, 3], flux=[1, 1, 1])
    assert_almost_equal(lc.fold(period=1).time[0], 0)
    assert_almost_equal(lc.fold(period=1, phase=-0.1).time[0], 0.1)


def test_cdpp():
    """Test the basics of the CDPP noise metric."""
    # A flat lightcurve should have a CDPP close to zero
    assert_almost_equal(LightCurve(np.arange(200), np.ones(200)).cdpp(norm_factor=1), 0)
    # An artificial lightcurve with sigma=100ppm should have cdpp=100ppm
    lc = LightCurve(np.arange(10000), np.random.normal(loc=1, scale=100e-6, size=10000))
    assert_almost_equal(lc.cdpp(transit_duration=1, norm_factor=1), 100, decimal=-0.5)


def test_cdpp_tabby():
    """Compare the cdpp noise metric against the pipeline value."""
    lcf = KeplerLightCurveFile(TABBY_Q8)
    # Tabby's star shows dips after cadence 1000 which may confuse the cdpp
    lc = LightCurve(lcf.PDCSAP_FLUX.time, lcf.PDCSAP_FLUX.flux)
    assert_almost_equal(lc.cdpp(), lcf.header(ext=1)['CDPP6_0'])


def test_lightcurve_plot():
    """Sanity check to verify that lightcurve plotting works"""
    lcf = KeplerLightCurveFile(TABBY_Q8)
    lcf.plot()
    lcf.SAP_FLUX.plot()
