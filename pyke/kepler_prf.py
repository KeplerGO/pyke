from . import DEFAULT_PRFDIR
from abc import ABC, abstractmethod
import os
import glob
import math
import scipy
import numpy as np
from astropy.io import fits as pyfits
from oktopus.models import get_initial_guesses
from oktopus.likelihood import PoissonLikelihood
from pyke.utils import channel_to_module_output


__all__ = ['KeplerPRFPhotometry', 'KeplerPRF']


class PRFPhotometry(ABC):
    """An abstract base class for a general PRF/PSF photometry algorithm
    for target pixel files."""

    @abstractmethod
    def do_photometry(self, tpf, initial_guesses=None):
        """Perform photometry on a given target pixel file.

        Parameters
        ----------
        tpf : pyke.TargetPixelFile instance
            A target pixel file instance
        initial_guesses : None or array-like
            A vector of initial estimates for the PRF/PSF model
        """
        pass

    @abstractmethod
    def generate_residuals_movie(self):
        """Creates a movie showing the residuals (image - fitted stars)
        for every cadence.
        """
        pass


class PRFModel(ABC):
    """An abstract base class for a general PRF/PSF parametric model."""

    @abstractmethod
    def evaluate(self, params):
        """Builds the PRF model parametrized by params.

        Parameters
        ----------
        *params : list-like
            Parameter values used to build a PRF model.

        Returns
        -------
        prf_model : 2D array
            PRF/PSF model.
        """
        pass


class KeplerPRFPhotometry(PRFPhotometry):
    """
    This class performs PRF Photometry on a target pixel file from
    NASA's Kepler/K2 missions.

    Attributes
    ----------
    prf_model : instance of PRFModel
    """
    # Let's borrow as much as possible from photutils here. Ideally,
    # this could be a child class from BasicPSFPhotometry.

    def __init__(self, prf_model, loss_function=PoissonLikelihood):
        self.prf_model = prf_model
        self.loss_function = loss_function
        self.opt_params = []
        self.residuals = []
        self.uncertainties = []

    def do_photometry(self, tpf, initial_guesses=None):
        if initial_guesses is None:
            # this must be clever enough to find the number of stars
            # great way to go is to use photutils.detection.DAOStarFinder
            initial_guesses, _ = get_inital_guesses(tpf.flux)

        for t in range(len(tpf.time)):
            logL = self.loss_function(tpf.flux[t], self.prf_model)
            opt_result = logL.fit(*initial_guesses).x
            residuals_opt_result = tpf.flux - self.prf_model(*opt_result)
            initial_guesses = opt_result
            self.opt_params.append(opt_result)
            self.residuals.append(residuals_opt_result)
            self.uncertainties.append(logL.uncertainties())

        self.opt_params = self.opt_params.reshape((tpf.shape[0], len(initial_guesses)))
        self.uncertainties = self.uncertainties.reshape((tpf.shape[0], len(initial_guesses)))

    def generate_residuals_movie(self):
        pass


class KeplerPRF(object):
    """
    Kepler's Pixel Response Function

    This class provides the necessary interface to load Kepler PSF
    calibration files and to create a model that can be fit as a function
    of flux and centroid position.

    Attributes
    ----------
    prf_files_dir : str
        Relative or aboslute path to a directory containing the Pixel Response
        Function calibration files produced during Kepler data comissioning.
    channel : int
        KeplerTargetPixelFile.channel
    shape : (int, int)
        KeplerTargetPixelFile.shape
    column : int
        KeplerTargetPixelFile.column
    row : int
        KeplerTargetPixelFile.row
    """

    def __init__(self, channel, shape, column, row, prf_files_dir=DEFAULT_PRFDIR):
        self.prf_files_dir = prf_files_dir
        self.channel = channel
        self.shape = shape
        self.column = column
        self.row = row
        self._prepare_prf()

    def prf_to_detector(self, params):
        """
        Interpolates the PRF model onto detector coordinates.

        Parameters
        ----------
        flux : float or array-like
            Total integrated flux of the PRF
        centroid_col : float or array-like
            Column coordinate of the centroid
        centroid_row : float or array-like
            Row coordinate of the centroid

        Returns
        -------
        prf_model : 2D array
            Two dimensional array representing the PRF values parametrized
            by `params`.
        """
        nsrcs = len(params) // 3
        if nsrcs > 1:
            F, xo, yo, self.prf_model = np.zeros(nsrcs), np.zeros(nsrcs), np.zeros(nsrcs), []
            for i in range(nsrcs):
                F[i] = params[i]
                xo[i] = params[i + nsrcs]
                yo[i] = params[i + 2 * nsrcs]
                self.prf_model.append(F[i] * self.interpolate(self.y - yo[i], self.x - xo[i]))
            return np.sum(self.prf_model, axis=0)
        else:
            F, xo, yo = params
            self.prf_model = F * self.interpolate(self.y - yo, self.x - xo)
            return self.prf_model

    def evaluate(self, *params):
        return self.prf_to_detector(params[:-1]) + params[-1]

    def _read_prf_calibration_file(self, path, ext):
        prf_cal_file = pyfits.open(path)
        data = prf_cal_file[ext].data
        # looks like these data below are the same for all prf calibration files
        crval1p = prf_cal_file[ext].header['CRVAL1P']
        crval2p = prf_cal_file[ext].header['CRVAL2P']
        cdelt1p = prf_cal_file[ext].header['CDELT1P']
        cdelt2p = prf_cal_file[ext].header['CDELT2P']
        prf_cal_file.close()
        return data, crval1p, crval2p, cdelt1p, cdelt2p

    def _prepare_prf(self):
        n_hdu = 5
        min_prf_weight = 1e-6
        module, output = channel_to_module_output(self.channel)
        # determine suitable PRF calibration file
        if module < 10:
            prefix = 'kplr0'
        else:
            prefix = 'kplr'
        prf_file_path = os.path.join(self.prf_files_dir,
                                     prefix + str(module) + '.' + str(output) + '*_prf.fits')
        prffile = glob.glob(prf_file_path)[0]

        # read PRF images
        prfn = [0] * n_hdu
        crval1p = np.zeros(n_hdu, dtype='float32')
        crval2p = np.zeros(n_hdu, dtype='float32')
        cdelt1p = np.zeros(n_hdu, dtype='float32')
        cdelt2p = np.zeros(n_hdu, dtype='float32')
        for i in range(n_hdu):
            prfn[i], crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i] = self._read_prf_calibration_file(prffile, i+1)
        prfn = np.array(prfn)
        PRFx = np.arange(0.5, np.shape(prfn[0])[1] + 0.5)
        PRFy = np.arange(0.5, np.shape(prfn[0])[0] + 0.5)
        PRFx = (PRFx - np.size(PRFx) / 2) * cdelt1p[0]
        PRFy = (PRFy - np.size(PRFy) / 2) * cdelt2p[0]

        # interpolate the calibrated PRF shape to the target position
        ydim, xdim = self.shape[0], self.shape[1]
        prf = np.zeros(np.shape(prfn[0]), dtype='float32')
        prfWeight = np.zeros(n_hdu, dtype='float32')
        ref_column = self.column + (xdim - 1.) / 2.
        ref_row = self.row + (ydim - 1.) / 2.
        for i in range(n_hdu):
            prfWeight[i] = math.sqrt((ref_column - crval1p[i]) ** 2
                                     + (ref_row - crval2p[i]) ** 2)
            if prfWeight[i] < min_prf_weight:
                prfWeight[i] = min_prf_weight
            prf += prfn[i] / prfWeight[i]
        prf /= (np.nansum(prf) * cdelt1p[0] * cdelt2p[0])

        # location of the data image centered on the PRF image (in PRF pixel units)
        self.x = np.arange(self.column + .5, self.column + xdim + .5)
        self.y = np.arange(self.row + .5, self.row + ydim + .5)
        self.interpolate = scipy.interpolate.RectBivariateSpline(PRFx, PRFy, prf)
