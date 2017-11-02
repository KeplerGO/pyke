import numpy as np
import scipy
import matplotlib.pyplot as plt
from astropy.io import fits
from .lightcurve import LightCurve
from .utils import KeplerQualityFlags, plot_image


__all__ = ['KeplerTargetPixelFile']


class TargetPixelFile(object):
    """
    TargetPixelFile class
    """
    def to_lightcurve(self, method=None, subtract_bkg=False, **kwargs):
        """Returns a raw light curve of the TPF.

        Attributes
        ----------
        method : str or None
            Method to detrend the light curve.
        kwargs : dict
            Keyword arguments passed to the detrending method.

        Returns
        -------
        lc : LightCurve object
            Array containing the summed or detrended flux within the aperture
            for each cadence.
        """
        pass


class KeplerTargetPixelFile(TargetPixelFile):
    """
    Defines a TargetPixelFile class for the Kepler/K2 Mission.
    Enables extraction of raw lightcurves and centroid positions.

    Attributes
    ----------
    path : str
        Path to fits file.

    quality_bitmask : int
        Bitmask specifying quality flags of cadences that should be ignored.

    References
    ----------
    .. [1] Kepler: A Search for Terrestrial Planets. Kepler Archive Manual.
        http://archive.stsci.edu/kepler/manuals/archive_manual.pdf
    """

    def __init__(self, path, aperture_mask=None,
                 quality_bitmask=KeplerQualityFlags.DEFAULT_BITMASK,
                 **kwargs):
        self.path = path
        self.hdu = fits.open(self.path, **kwargs)
        self.quality_bitmask = quality_bitmask
        self.quality_mask = self._quality_mask(quality_bitmask)
        self.aperture_mask = None

    def _quality_mask(self, quality_bitmask):
        """Returns a boolean mask which flags all good-quality cadences.

        Parameters
        ----------
        quality_bitmask : int
            Bitmask. See ref. [1], table 2-3.
        """
        return (self.hdu[1].data['QUALITY'] & quality_bitmask) == 0

    @property
    def keplerid(self):
        return self.hdu[0].header['KEPLERID']

    @property
    def module(self):
        return self.hdu[0].header['MODULE']

    @property
    def channel(self):
        return self.hdu[0].header['CHANNEL']

    @property
    def output(self):
        return self.hdu[0].header['OUTPUT']

    @property
    def column(self):
        return self.hdu['TARGETTABLES'].header['1CRV5P']

    @property
    def row(self):
        return self.hdu['TARGETTABLES'].header['2CRV5P']

    def plot(self, nframe=100, **kwargs):
        pflux = self.flux[nframe]
        plot_image(pflux, title='Kepler ID: {}'.format(self.keplerid),
                   extent=(self.column, self.column + self.shape[2],
                           self.row, self.row + self.shape[1]), **kwargs)

    @property
    def aperture_mask(self):
        return self._aperture_mask

    @aperture_mask.setter
    def aperture_mask(self, mask):
        if mask is not None:
            self._aperture_mask = mask
        else:
            self._aperture_mask = np.ones((self.shape[1], self.shape[2]),
                                          dtype=bool)

    @property
    def aperture_npix(self):
        """Number of pixels in the aperture"""
        return self.aperture_mask.sum()

    @property
    def n_good_cadences(self):
        """Returns the number of good-quality cadences."""
        return self.quality_mask.sum()

    @property
    def shape(self):
        """Return the cube dimension shape."""
        return self.flux.shape

    @property
    def time(self):
        """Returns the time for all good-quality cadences."""
        return self.hdu[1].data['TIME'][self.quality_mask]

    @property
    def nan_time_mask(self):
        """Returns a boolean mask flagging cadences whose time is `nan`."""
        return ~np.isfinite(self.time)

    @property
    def flux(self):
        """Returns the flux for all good-quality cadences."""
        return self.hdu[1].data['FLUX'][self.quality_mask] + self.flux_bkg

    @property
    def flux_err(self):
        """Returns the flux uncertainty for all good-quality cadences."""
        return self.hdu[1].data['FLUX_ERR'][self.quality_mask]

    @property
    def flux_bkg(self):
        """Returns the background flux for all good-quality cadences."""
        return self.hdu[1].data['FLUX_BKG'][self.quality_mask]

    @property
    def quality(self):
        """Returns the quality flag integer of every good cadence."""
        return self.hdu[1].data['QUALITY'][self.quality_mask]

    def estimate_background_per_pixel(self, method='mode'):
        """Returns the median value of the fluxes for every cadence as an
        estimate for the background density."""
        if method == 'median':
            return np.nanmedian(self.flux[:, self.aperture_mask], axis=1)
        elif method == 'mean':
            return np.nanmean(self.flux[:, self.aperture_mask], axis=1)
        elif method == 'mode':
            return scipy.stats.mode(self.flux[:, self.aperture_mask], axis=1,
                                    nan_policy='omit')[0]
        else:
            raise ValueError("method {} is not available".format(method))

    def to_fits(self):
        """Save the TPF to fits"""
        raise NotImplementedError

    def _get_aperture_flux(self):
        return np.nansum(self.flux[:, self.aperture_mask], axis=1)

    def get_bkg_lightcurve(self):
        return np.nansum(self.flux_bkg[:, self.aperture_mask], axis=1)

    def to_lightcurve(self, subtract_bkg=False):
        """Performs apperture photometry and optionally detrends the lightcurve.

        Attributes
        ----------
        subtract_bkg : bool
            Whether or not to subtract the background.

        Returns
        -------
        lc : LightCurve object
            Array containing the summed flux within the aperture for each
            cadence.
        """

        aperture_flux = self._get_aperture_flux()
        if subtract_bkg:
            aperture_flux = aperture_flux - self.get_bkg_lightcurve()

        return LightCurve(flux=aperture_flux, time=self.time)
