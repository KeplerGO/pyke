import numpy as np
from astropy.io import fits
from astropy.stats.funcs import median_absolute_deviation
import scipy.ndimage
from .lightcurve import LightCurve


__all__ = ['KeplerTargetPixelFile']


class TargetPixelFile(object):
    """
    TargetPixelFile class
    """

    def open(self, **kwargs):
        pass

    def get_data(self):
        pass

    def to_lightcurve(self, aperture_mask=None, method=None, **kwargs):
        """Returns a raw light curve of the TPF.

        Attributes
        ----------
        aperture_mask: boolean ndarray or None
            Aperture under which the flux will be summed up. If ``None``,
            then an aperture is computed using ``aperture_mask``.
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

    max_quality : int
        Maximum tolerated quality for the cadences.

    References
    ----------
    .. [1] Kepler: A Search for Terrestrial Planets. Kepler Archive Manual.
        http://archive.stsci.edu/kepler/manuals/archive_manual.pdf
    """

    def __init__(self, path, max_quality=1, aperture_mask=None, **kwargs):
        self.path = path
        self.hdu = fits.open(self.path, **kwargs)
        self.max_quality = max_quality
        self.aperture_mask = None
        self._aperture_flux = None
        self._good_quality_cadences = self.good_quality_cadences()

    def good_quality_cadences(self, max_quality=None):
        """Returns a boolean mask flagging cadences whose quality is at most
        `max_quality`.

        Parameters
        ----------
        max_quality : int or None
            Maximum tolerated quality. See ref. [1], table 2-3.
        """

        if max_quality is None:
            max_quality = self.max_quality
        return self.quality < max_quality

    @property
    def module(self):
        return self.hdu[0].header['MODULE']

    @property
    def channel(self):
        return self.hdu[0].header['CHANNEL']

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
    def n_cadences(self):
        """Returns the number of good-quality cadences."""
        return self.good_quality_cadences.sum()

    @property
    def shape(self):
        """Return the cube dimension shape."""
        return self.flux.shape

    @property
    def time(self):
        """Returns the time for all good-quality cadences."""
        return self.hdu[1].data['TIME'][self.good_quality_cadences]

    @property
    def nan_time_mask(self):
        """Returns a boolean mask flagging cadences whose time is `nan`."""
        return ~np.isfinite(self.time)

    @property
    def flux(self):
        """Returns the flux for all good-quality cadences."""
        return self.hdu[1].data['FLUX'][self.good_quality_cadences]

    @property
    def bkg(self):
        """Returns the median value of the fluxes for every cadence as an
        estimate for the background."""
        return np.nanmedian(self.flux[:, self.aperture_mask], axis=1)

    @property
    def quality(self):
        """Returns the quality flag integer of every cadence."""
        return self.hdu[1]['QUALITY']

    def to_fits(self):
        """Save the TPF to fits"""
        raise NotImplementedError

    def make_aperture_mask(self, snr_threshold=5, margin=4):
        """Returns an aperture photometry mask.

        Parameters
        ----------
        snr_threshold : float
            Background detection threshold.
        """

        # Find the pixels that are above the threshold in the median flux image
        median = np.nanmedian(self.flux, axis=0)
        mad = median_absolute_deviation(median[np.isfinite(median)])
        # 1.4826 turns MAD into STDEV for a Gaussian
        mad_cut = 1.4826 * mad * snr_threshold

        region = np.where(median > mad_cut, 1, 0)
        # Label all contiguous regions above the threshold
        labels = scipy.ndimage.label(region)[0]
        # Central pixel coordinate
        centralpix = [1 + median.shape[0] // 2, 1 + median.shape[1] // 2]

        # find brightest pix within margin of central pix
        central_img = median[centralpix[0] - margin: centralpix[0] + margin,
                             centralpix[1] - margin: centralpix[1] + margin]
        # unravel_index converts indices into a tuple of coordinate arrays
        brightestpix = np.unravel_index(central_img.argmax(), central_img.shape)
        bpixy, bpixx = brightestpix

        # Which label corresponds to the brightest pixel?
        regnum = labels[centralpix[0] - margin + bpixy, centralpix[1] - margin + bpixx]

        return labels == regnum

    def centroids(self):
        """Returns the centroids for every cadence under a given aperture
        mask.

        Attributes
        ----------

        Returns
        -------
        xc, yc: ndarrays
            centroid positions for every cadence
        """

        xc = np.zeros(self.n_cadences)
        yc = np.zeros(self.n_cadences)
        y, x = np.mgrid[:img.shape[0], :img.shape[1]]

        x = x[self.aperture_mask]
        y = y[self.aperture_mask]

        if self._aperture_flux is None:
            self._aperture_flux = self._get_aperture_flux()

        for i in range(self.n_cadences):
            flux_i = self.flux[i][self.aperture_mask]
            xc[i] = np.nansum(flux_i * x) / self._aperture_flux[i]
            yc[i] = np.nansum(flux_i * y) / self._aperture_flux[i]

        return xc, yc

    def _get_aperture_flux(self):
        return np.nansum(self.flux[:, self.aperture_mask], axis=1)

    def to_lightcurve(self, method=None, subtract_bkg=False, **kwargs):
        """Performs aperture photometry and optionally detrends the lightcurve.
        """

        if self._aperture_flux is None:
            self._aperture_flux = self._get_aperture_flux()

        if subtract_bkg:
            # number of pixels in the aperture
            npix = self.aperture_mask.sum()
            self._aperture_flux = self._aperture_flux - npix * self.bkg

        if method is None:
            return LightCurve(flux=self._aperture_flux, time=self.time)
        else:
            return LightCurve(flux=self._aperture_flux,
                              time=self.time).detrend(method=method, **kwargs)
