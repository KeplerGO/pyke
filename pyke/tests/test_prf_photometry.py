import pytest
import math
import numpy as np
from scipy.stats import mode
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from oktopus import PoissonPosterior, UniformPrior, GaussianPrior, JointPrior
from ..prf_photometry import KeplerPRF, SceneModel, PRFPhotometry, get_initial_guesses


def test_prf_normalization():
    """Does the PRF model integrate to the requested flux across the focal plane?"""
    for channel in [1, 20, 40, 60, 84]:
        for col in [123, 678]:
            for row in [234, 789]:
                shape = (18, 14)
                flux = 100
                prf = KeplerPRF(channel=channel, column=col, row=row, shape=shape)
                prf_sum = prf.evaluate(flux, col + shape[0]/2, row + shape[1]/2, 1, 1, 0).sum()
                assert np.isclose(prf_sum, flux, rtol=0.1)


def test_prf_vs_aperture_photometry():
    """Is the PRF photometry result consistent with simple aperture photometry?"""
    tpf_fn = get_pkg_data_filename("data/ktwo201907706-c01-first-cadence.fits.gz")
    tpf = fits.open(tpf_fn)
    col, row = 173, 526
    prf = KeplerPRF(channel=tpf[0].header['CHANNEL'],
                    column=col, row=row,
                    shape=tpf[1].data.shape)
    scene = SceneModel(prfs=prf)
    fluxo, colo, rowo, _ = get_initial_guesses(data=tpf[1].data,
                                                 ref_col=prf.col_coord[0],
                                                 ref_row=prf.row_coord[0])
    bkg = mode(tpf[1].data, None)[0]
    prior = JointPrior(UniformPrior(lb=0.1*fluxo, ub=fluxo),
                       UniformPrior(lb=prf.col_coord[0], ub=prf.col_coord[-1]),
                       UniformPrior(lb=prf.row_coord[0], ub=prf.row_coord[-1]),
                       GaussianPrior(mean=1, var=1e-2),
                       GaussianPrior(mean=1, var=1e-2),
                       GaussianPrior(mean=0, var=1e-2),
                       UniformPrior(lb=.2*bkg, ub=5*bkg))
    logL = PoissonPosterior(tpf[1].data, mean=scene, prior=prior)
    result = logL.fit(x0=prior.mean, method='powell')
    prf_flux, prf_col, prf_row, prf_scale_col, prf_scale_row, prf_rotation, prf_bkg = logL.opt_result.x
    assert result.success is True
    assert np.isclose(prf_col, colo, rtol=1e-1)
    assert np.isclose(prf_row, rowo, rtol=1e-1)
    assert np.isclose(prf_bkg, np.percentile(tpf[1].data, 10), rtol=0.1)

    # Test KeplerPRFPhotometry class
    kepler_phot = PRFPhotometry(scene_model=scene, prior=prior)
    tpf_flux = tpf[1].data.reshape((1, tpf[1].data.shape[0], tpf[1].data.shape[1]))
    kepler_phot.fit(tpf_flux=tpf_flux)
    opt_params = kepler_phot.opt_params.reshape(-1)
    assert np.isclose(opt_params[0], prf_flux, rtol=0.1)
    assert np.isclose(opt_params[1], prf_col, rtol=1e-1)
    assert np.isclose(opt_params[2], prf_row, rtol=1e-1)
    assert np.isclose(opt_params[-1], prf_bkg, rtol=0.1)
