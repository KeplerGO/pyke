import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import (PercentileInterval, ImageNormalize,
                                   SqrtStretch, LogStretch, LinearStretch)

from .lightcurve import LightCurve


__all__ = ['KeplerTargetPixelFile', 'KeplerQualityFlags']


class KeplerQualityFlags(object):
    """
    This class encodes the meaning of the various Kepler QUALITY bitmask flags,
    as documented in the Kepler Archive Manual (Table 2.3).
    """
    AttitudeTweak = 1
    SafeMode = 2
    CoarsePoint = 4
    EarthPoint = 8
    ZeroCrossing = 16
    Desat = 32
    Argabrightening = 64
    ApertureCosmic = 128
    ManualExclude = 256
    SensitivityDropout = 1024
    ImpulsiveOutlier = 2048
    ArgabrighteningOnCCD = 4096
    CollateralCosmic = 8192
    DetectorAnomaly = 16384
    NoFinePoint = 32768
    NoData = 65536
    RollingBandInAperture = 131072
    RollingBandInMask = 262144
    PossibleThrusterFiring = 524288
    ThrusterFiring = 1048576

    # Which is the recommended QUALITY mask to identify bad data?
    DEFAULT_BITMASK = (AttitudeTweak | SafeMode | CoarsePoint | EarthPoint |
                       Desat | ApertureCosmic | ManualExclude | NoData | ThrusterFiring)

    # Pretty string descriptions for each flag
    STRINGS = {
        1: "Attitude tweak",
        2: "Safe mode",
        4: "Coarse point",
        8: "Earth point",
        16: "Zero crossing",
        32: "Desaturation event",
        64: "Argabrightening",
        128: "Cosmic ray in optimal aperture",
        256: "Manual exclude",
        1024: "Sudden sensitivity dropout",
        2048: "Impulsive outlier",
        4096: "Argabrightening on CCD",
        8192: "Cosmic ray in collateral data",
        16384: "Detector anomaly",
        32768: "No fine point",
        65536: "No data",
        131072: "Rolling band in optimal aperture",
        262144: "Rolling band in full mask",
        524288: "Possible thruster firing",
        1048576: "Thruster firing"
    }

    @classmethod
    def decode(cls, quality):
        """Converts a Kepler QUALITY value into a list of human-readable strings.

        This function takes the QUALITY bitstring that can be found for each
        cadence in Kepler/K2's pixel and light curve files and converts into
        a list of human-readable strings explaining the flags raised (if any).

        Parameters
        ----------
        quality : int
            Value from the 'QUALITY' column of a Kepler/K2 pixel or lightcurve file.

        Returns
        -------
        flags : list of str
            List of human-readable strings giving a short description of the
            quality flags raised.  Returns an empty list if no flags raised.
        """
        result = []
        for flag in cls.STRINGS.keys():
            if quality & flag > 0:
                result.append(cls.STRINGS[flag])
        return result


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
        self._aperture_flux = None

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

    def plot(self, nframe=100, scale='linear', **kwargs):
        pflux = self.flux[nframe]
        vmin, vmax = PercentileInterval(95.).get_limits(pflux)
        if scale == 'linear':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
        elif scale == 'sqrt':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
        elif scale == 'log':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
        else:
            raise ValueError("scale {} is not available.".format(scale))

        plt.imshow(pflux, origin='lower', norm=norm,
                   extent=(self.column, self.column + self.shape[2],
                           self.row, self.row + self.shape[1]), **kwargs)
        plt.xlabel('Pixel Column Number')
        plt.ylabel('Pixel Row Number')
        plt.title('Kepler ID: {}'.format(self.keplerid))
        plt.colorbar(norm=norm)

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
        return self.hdu[1].data['FLUX'][self.quality_mask]

    @property
    def flux_err(self):
        """Returns the flux for all good-quality cadences."""
        return self.hdu[1].data['FLUX_ERR'][self.quality_mask]

    @property
    def quality(self):
        """Returns the quality flag integer of every good cadence."""
        return self.hdu[1].data['QUALITY'][self.quality_mask]

    def estimate_background(self):
        """Returns the median value of the fluxes for every cadence as an
        estimate for the background."""
        return np.nanmedian(self.flux[:, self.aperture_mask], axis=1)

    def to_fits(self):
        """Save the TPF to fits"""
        raise NotImplementedError

    def _get_aperture_flux(self):
        return np.nansum(self.flux[:, self.aperture_mask], axis=1)

    def to_lightcurve(self, method=None, subtract_bkg=False, **kwargs):
        """Performs apperture photometry and optionally detrends the lightcurve.

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

        if self._aperture_flux is None:
            self._aperture_flux = self._get_aperture_flux()

        if subtract_bkg:
            # number of pixels in the aperture
            self._aperture_flux = self._aperture_flux - self.aperture_npix * self.estimate_background

        if method is None:
            return LightCurve(flux=self._aperture_flux, time=self.time)
        else:
            return LightCurve(flux=self._aperture_flux,
                              time=self.time).detrend(method=method, **kwargs)
