from . import DEFAULT_PRFDIR
from .utils import channel_to_module_output
from abc import abstractmethod
import os
import glob
import math
import scipy
import numpy as np
import tqdm
from astropy.io import fits as pyfits
from oktopus.posterior import PoissonPosterior


__all__ = ['KeplerPRFPhotometry', 'KeplerSceneModel', 'KeplerPRF', 'get_initial_guesses']


class KeplerPRFPhotometry(object):
    """
    This class performs PRF Photometry on a target pixel file from
    NASA's Kepler/K2 missions.

    Attributes
    ----------
    scene_model : instance of KeplerSceneModel
    priors : instance of oktopus.JointPrior
        Priors on the parameters that will be estimated
    loss_function : subclass of oktopus.LossFunction
        Noise distribution associated with each random measurement
    """

    def __init__(self, scene_model, prior, loss_function=PoissonPosterior):
        self.scene_model = scene_model
        self.prior = prior
        self.loss_function = loss_function
        self.opt_params = np.array([])
        self.residuals = np.array([])
        self.uncertainties = np.array([])

    def fit(self, tpf_flux, x0=None, cadences='all', **kwargs):
        """
        Fits the scene model to the given data in ``tpf_flux``.

        Parameters
        ----------
        tpf_flux : array-like
            A pixel flux time-series, e.g, KeplerTargetPixelFile.flux.
            Such that (time, row, column) represents the shape of ``tpf_flux``.
        x0 : array-like or None
            Initial guesses on the parameters.
        cadences : array-like of ints or str
            A list or array that contains the cadences which will be fitted.
            Default is to fit all cadences.
        kwargs : dict
            Dictionary of additional parameters to be passed to
            `scipy.optimize.minimize`.
        """
        if x0 is None:
            x0 = self.prior.mean

        if cadences == 'all':
            cadences = range(tpf_flux.shape[0])

        for t in tqdm.tqdm(cadences):
            loss = self.loss_function(tpf_flux[t], self.scene_model, prior=self.prior)
            opt_params = loss.fit(x0=x0, **kwargs).x
            residuals = tpf_flux[t] - self.scene_model(*opt_params)
            self.opt_params = np.append(self.opt_params, opt_params)
            self.residuals = np.append(self.residuals, residuals)
        self.opt_params = self.opt_params.reshape((tpf_flux.shape[0], len(x0)))
        self.residuals = self.residuals.reshape(tpf_flux.shape)

    def generate_residuals_movie(self):
        pass


class KeplerSceneModel(object):
    """
    This class builds a generic model for a Kepler scene.

    Attributes
    ----------
    prf_model : instance of KeplerPRF
        An instance of the KeplerPRF class
    n_sources : int
        Number of sources to be modeled with ``prf_model``
    bkg_model : callable
        A function that models the background variation.
        Default is a constant background
    """

    def __init__(self, prf_model, n_sources, bkg_model=lambda bkg: np.array([bkg])):
        self.prf_model = prf_model
        self.n_sources = n_sources
        self.bkg_model = bkg_model

    def __call__(self, *params):
        return self.evaluate(*params)

    def evaluate(self, *params):
        """
        Parameters
        ----------
        flux : scalar or array-like
            Total integrated flux of the PRF model
        centroid_col, centroid_row : scalar or array-like
            Column and row coordinates of the centroid
        scale_col, scale_row : scalar or array-like
            Pixel scale in the column and row directions
        bkg_params : scalar or array-like
            Parameters for the background model
        """
        self.mixture_model = []
        for i in range(self.n_sources):
            self.mixture_model.append(self.prf_model(params[i],
                                                     params[i + self.n_sources],
                                                     params[i + 2 * self.n_sources]))
        self.scene_model = np.sum(self.mixture_model, axis=0) + self.bkg_model(params[-1])

        return self.scene_model


class KeplerPRF(object):
    """
    Kepler's Pixel Response Function

    This class provides the necessary interface to load Kepler PRF
    calibration files and to create a model that can be fit as a function
    of flux, centroid positions, width, and rotation angle.

    Attributes
    ----------
    channel : int
        KeplerTargetPixelFile.channel
    shape : (int, int)
        KeplerTargetPixelFile.shape[1:]
    column : int
        KeplerTargetPixelFile.column
    row : int
        KeplerTargetPixelFile.row
    prf_files_dir : str
        Relative or aboslute path to a directory containing the Pixel Response
        Function calibration files produced during Kepler data comissioning.

    Notes
    -----
    The PRF calibration files can be downloaded here: http://archive.stsci.edu/missions/kepler/fpc/prf/
    """

    def __init__(self, channel, shape, column, row, prf_files_dir=DEFAULT_PRFDIR):
        self.prf_files_dir = prf_files_dir
        self.channel = channel
        self.shape = shape
        self.column = column
        self.row = row
        self.col_coord, self.row_coord, self.interpolate = self._prepare_prf()

    def __call__(self, *params):
        return self.evaluate(*params)

    def evaluate(self, *params):
        return self.prf_to_detector(*params)

    def prf_to_detector(self, flux, centroid_col, centroid_row, scale_col,
                        scale_row, rotation_angle):
        """
        Interpolates the PRF model onto detector coordinates.

        Parameters
        ----------
        flux : float or array-like
            Total integrated flux of the PRF
        centroid_col, centroid_row : float or array-like
            Column and row coordinates of the centroid
        scale_col, scale_row : float or array-like
            Pixel scale in the column and row directions
        rotation_angle : float
            Rotation angle in radians

        Returns
        -------
        prf_model : 2D array
            Two dimensional array representing the PRF values parametrized
            by `params`.
        """
        cosa = math.cos(rotation_angle)
        sina = math.sin(rotation_angle)

        delta_col = self.col_coord - centroid_col
        delta_row = self.row_coord - centroid_row

        delta_row, delta_col = np.meshgrid(delta_row, delta_col)
        rot_col = delta_col * cosa - delta_row * sina
        rot_row = delta_col * sina + delta_row * cosa

        self.prf_model = self.interpolate(rot_row.flatten() * scale_row,
                                          rot_col.flatten() * scale_col, grid=False).reshape(self.shape)
        return self.prf_model

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
        PRFcol = np.arange(0.5, np.shape(prfn[0])[1] + 0.5)
        PRFrow = np.arange(0.5, np.shape(prfn[0])[0] + 0.5)
        PRFcol = (PRFcol - np.size(PRFcol) / 2) * cdelt1p[0]
        PRFrow = (PRFrow - np.size(PRFrow) / 2) * cdelt2p[0]

        # interpolate the calibrated PRF shape to the target position
        rowdim, coldim = self.shape[0], self.shape[1]
        prf = np.zeros(np.shape(prfn[0]), dtype='float32')
        prfWeight = np.zeros(n_hdu, dtype='float32')
        ref_column = self.column + (coldim - 1.) / 2.
        ref_row = self.row + (rowdim - 1.) / 2.
        for i in range(n_hdu):
            prfWeight[i] = math.sqrt((ref_column - crval1p[i]) ** 2
                                     + (ref_row - crval2p[i]) ** 2)
            if prfWeight[i] < min_prf_weight:
                prfWeight[i] = min_prf_weight
            prf += prfn[i] / prfWeight[i]
        prf /= (np.nansum(prf) * cdelt1p[0] * cdelt2p[0])

        # location of the data image centered on the PRF image (in PRF pixel units)
        col_coord = np.arange(self.column + .5, self.column + coldim + .5)
        row_coord = np.arange(self.row + .5, self.row + rowdim + .5)
        interpolate = scipy.interpolate.RectBivariateSpline(PRFcol, PRFrow, prf)

        return col_coord, row_coord, interpolate


def get_initial_guesses(data, ref_col, ref_row):
    """
    Compute the initial guesses for total flux, centroids position, and PSF
    width using the sample moments of the data.

    Parameters
    ----------
    data : 2D array-like
        Image data
    ref_col, ref_row : scalars
        Reference column and row (coordinates of the bottom left corner)

    Return
    ------
    flux0, col0, row0, sigma0: floats
        Inital guesses for flux, centroid position, and width
    """

    flux0 = np.nansum(data)
    yy, xx = np.indices(data.shape)
    yy = ref_row + yy
    xx = ref_col + xx
    col0 = np.nansum(xx * data) / flux0
    row0 = np.nansum(yy * data) / flux0
    marg_col = data[:, int(np.round(col0 - ref_col))]
    marg_row = data[int(np.round(row0 - ref_row)), :]
    sigma_y = math.sqrt(np.abs((np.arange(marg_row.size) - row0) ** 2 * marg_row).sum() / marg_row.sum())
    sigma_x = math.sqrt(np.abs((np.arange(marg_col.size) - col0) ** 2 * marg_col).sum() / marg_col.sum())
    sigma0 = math.sqrt((sigma_x**2 + sigma_y**2)/2.0)

    return flux0, col0, row0, sigma0
