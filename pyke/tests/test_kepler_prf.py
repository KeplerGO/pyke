import os
import math
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from .. import models
from ..likelihood import PoissonLikelihood
from ..posterior import PoissonPosterior
from ..prior import UniformPrior, GaussianPrior, JointPrior
from .. import core
from ..kepler_prf import KeplerPRF


def test_prf_normalization():
    """Does the PRF model integrate to the requested flux across the focal plane?"""
    for channel in [1, 20, 40, 60, 84]:
        for col in [123, 678]:
            for row in [234, 789]:
                shape = (18, 14)
                flux = 100
                prf = KeplerPRF(channel=channel, column=col, row=row, shape=shape)
                prf_sum = prf.evaluate(flux, col + shape[0]/2, row + shape[1]/2, 0).sum()
                assert np.isclose(prf_sum, flux, rtol=0.1)


def test_prf_vs_aperture_photometry():
    """Is the PRF photometry result consistent with simple aperture photometry?"""
    tpf_fn = get_pkg_data_filename("data/ktwo201907706-c01-first-cadence.fits.gz")
    tpf = fits.open(tpf_fn)
    col, row = 173, 526
    prf = KeplerPRF(channel=tpf[0].header['CHANNEL'],
                    column=col, row=row,
                    shape=tpf[1].data.shape)
    fluxo, colo, rowo, sigmao = models.get_initial_guesses(data=tpf[1].data,
                                                           X=prf.x,
                                                           Y=prf.y)
    bkgo = np.mean(tpf[1].data)
    aperture_flux = tpf[1].data.sum() - bkgo
    prior = JointPrior(GaussianPrior(mean=fluxo, var=math.sqrt(fluxo)),
                       UniformPrior(lb=prf.x[0], ub=prf.x[-1]),
                       UniformPrior(lb=prf.y[0], ub=prf.y[-1]),
                       GaussianPrior(mean=bkgo, var=bkgo))
    logL = PoissonPosterior(tpf[1].data, prf.evaluate, prior=prior)
    fitresult = logL.fit((fluxo, colo, rowo, bkgo))
    print(fitresult)
    prf_flux, prf_col, prf_row, prf_bkg = fitresult.x
    assert np.isclose(prf_col, col+9, rtol=1e-3)
    assert np.isclose(prf_row, row+9, rtol=1e-3)
    assert np.isclose(prf_bkg, np.percentile(tpf[1].data, 10), rtol=0.1)
    assert np.isclose(aperture_flux, prf_flux, rtol=0.1)
