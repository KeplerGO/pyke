from . import DEFAULT_PRFDIR
from .utils import channel_to_module_output
from abc import abstractmethod
import math
import scipy
import numpy as np
import tqdm
from inspect import signature
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
            result = loss.fit(x0=x0, **kwargs)
            opt_params = result.x
            residuals = tpf_flux[t] - self.scene_model(*opt_params)
            self.opt_params = np.append(self.opt_params, opt_params)
            self.residuals = np.append(self.residuals, residuals)
        self.opt_params = self.opt_params.reshape((tpf_flux.shape[0], len(x0)))
        self.residuals = self.residuals.reshape(tpf_flux.shape)

    def get_residuals(self):
        return self.residuals

    def get_fitted_parameters_matrix(self):
        return self.opt_params


class KeplerSceneModel(object):
    """
    This class builds a generic model for a Kepler scene.

    Attributes
    ----------
    prfs : list of callables
        A list of prfs
    bkg_model : callable
        A function that models the background variation.
        Default is a constant background
    """

    def __init__(self, prfs, bkg_model=lambda bkg: np.array([bkg])):
        self.prfs = np.asarray([prfs]).reshape(-1)
        self.bkg_model = bkg_model
        self._prepare_scene_model()

    def __call__(self, *params):
        return self.evaluate(*params)

    def _prepare_scene_model(self):
        self.n_models = len(self.prfs)
        self.bkg_order = len(signature(self.bkg_model).parameters)

        model_orders = [0]
        for i in range(self.n_models):
            model_orders.append(len(signature(self.prfs[i]).parameters))
        self.n_params = np.cumsum(model_orders)

    def evaluate(self, *params):
        """
        Parameters
        ----------
        flux : scalar or array-like
            Total integrated flux of the PRF model
        center_col, center_row : scalar or array-like
            Column and row coordinates of the center
        scale_col, scale_row : scalar or array-like
            Pixel scale in the column and row directions
        rotation_angle : float
            Rotation angle in radians
        bkg_params : scalar or array-like
            Parameters for the background model
        """
        self.mm = []
        for i in range(self.n_models):
            self.mm.append(self.prfs[i](*params[self.n_params[i]:self.n_params[i+1]]))
        self.scene_model = np.sum(self.mm, axis=0) + self.bkg_model(*params[-self.bkg_order:])

        return self.scene_model


class KeplerPRF(object):
    """
    Kepler's Pixel Response Function

    This class provides the necessary interface to load Kepler PRF
    calibration files and to create a model that can be fit as a function
    of flux, center positions, width, and rotation angle.

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

    Examples
    --------
    Objects from the KeplerPRF class are defined by a channel number, a pair of
    dimensions (the size of the image), and a reference coordinate (bottom left
    corner). In this example, we create a KeplerPRF object located at channel
    #44 with dimension equals 10 x 10, reference row and column coordinate
    equals (5, 5). After the object has been created, we may translate it to a
    given center coordinate. Additionally, we can specify total flux, pixel
    scales, and rotation around the object's center.

    >>> import math
    >>> import matplotlib.pyplot as plt
    >>> from pyke import KeplerPRF
    >>> kepprf = KeplerPRF(channel=44, shape=(10, 10), column=5, row=5)
    >>> prf = kepprf(flux=1000, center_col=10, center_row=10,
    ...              scale_row=0.7, scale_col=0.7, rotation_angle=math.pi/2)
    >>> plt.imshow(prf, origin='lower')
    """

    def __init__(self, channel, shape, column, row):
        self.channel = channel
        self.shape = shape
        self.column = column
        self.row = row
        self.col_coord, self.row_coord, self.interpolate = self._prepare_prf()

    def __call__(self, flux, center_col, center_row, scale_col, scale_row,
                 rotation_angle):
        return self.evaluate(flux, center_col, center_row,
                             scale_col, scale_row, rotation_angle)

    def evaluate(self, flux, center_col, center_row, scale_col, scale_row,
                 rotation_angle):
        """
        Interpolates the PRF model onto detector coordinates.

        Parameters
        ----------
        flux : float
            Total integrated flux of the PRF
        center_col, center_row : float
            Column and row coordinates of the center
        scale_col, scale_row : float
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

        delta_col = self.col_coord - center_col
        delta_row = self.row_coord - center_row
        delta_row, delta_col = np.meshgrid(delta_row, delta_col)

        rot_col = delta_col * cosa - delta_row * sina
        rot_row = delta_col * sina + delta_row * cosa

        self.prf_model = flux * self.interpolate(rot_row.flatten() * scale_row,
                                                 rot_col.flatten() * scale_col, grid=False).reshape(self.shape)
        return self.prf_model

    def gradient(self, flux, center_col, center_row):
        """
        This function returns the gradient of the KeplerPRF model with respect
        to flux, center_col, and center_row, in the particular case where the
        scales in the row and column dimensions are both equal to one and the
        rotation angle is zero.

        Parameters
        ----------
        flux : float
            Total integrated flux of the PRF
        center_col, center_row : float
            Column and row coordinates of the center

        Returns
        -------
        grad_prf : list
            Returns a list of arrays where the elements are the derivative
            of the KeplerPRF model with respect to flux, center_col, and
            center_row, respectively.
        """
        delta_col = self.col_coord - center_col
        delta_row = self.row_coord - center_row
        delta_row, delta_col = np.meshgrid(delta_row, delta_col)

        deriv_flux = self.interpolate(delta_col, delta_row, grid=False).reshape(self.shape)

        deriv_center_col = - flux * self.interpolate(delta_col, delta_row, dx=1,
                                                     grid=False).reshape(self.shape)

        deriv_center_row = - flux * self.interpolate(delta_col, delta_row, dy=1,
                                                     grid=False).reshape(self.shape)

        return [deriv_flux, deriv_center_col, deriv_center_row]


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
        prfs_url_path = "http://archive.stsci.edu/missions/kepler/fpc/prf/extracted/"
        prffile = prfs_url_path + prefix + str(module) + '.' + str(output) + '_2011265_prf.fits'

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
    Compute the initial guesses for total flux, centers position, and PSF
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
        Inital guesses for flux, center position, and width
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