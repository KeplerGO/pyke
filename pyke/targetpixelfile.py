import numpy as np
import scipy
import matplotlib.pyplot as plt
from astropy.io import fits
from .lightcurve import KeplerLightCurve, LightCurve
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
    quality_bitmask : str or int
        Bitmask specifying quality flags of cadences that should be ignored.
        If a string is passed, it has the following meaning:

            * "default": recommended quality mask
            * "hard": removes more flags, known to remove good data
            * "hardest": removes all data that has been flagged

    References
    ----------
    .. [1] Kepler: A Search for Terrestrial Planets. Kepler Archive Manual.
        http://archive.stsci.edu/kepler/manuals/archive_manual.pdf
    """

    def __init__(self, path, quality_bitmask=KeplerQualityFlags.DEFAULT_BITMASK,
                 **kwargs):
        self.path = path
        self.hdu = fits.open(self.path, **kwargs)
        self.quality_bitmask = quality_bitmask
        self.quality_mask = self._quality_mask(quality_bitmask)

    def _quality_mask(self, bitmask):
        """Returns a boolean mask which flags all good-quality cadences.

        Parameters
        ----------
        bitmask : str or int
            Bitmask. See ref. [1], table 2-3.
        """
        if bitmask is None:
            return np.ones(len(self.hdu[1].data['TIME']), dtype=bool)
        elif isinstance(bitmask, str):
            bitmask = KeplerQualityFlags.OPTIONS[bitmask]
        return (self.hdu[1].data['QUALITY'] & bitmask) == 0

    def header(self, ext=0):
        """Returns the header for a given extension."""
        return self.hdu[ext].header

    @property
    def keplerid(self):
        return self.header()['KEPLERID']

    @property
    def module(self):
        return self.header()['MODULE']

    @property
    def channel(self):
        return self.header()['CHANNEL']

    @property
    def output(self):
        return self.header()['OUTPUT']

    @property
    def column(self):
        return self.hdu['TARGETTABLES'].header['1CRV5P']

    @property
    def row(self):
        return self.hdu['TARGETTABLES'].header['2CRV5P']

    @property
    def pipeline_mask(self):
        """Returns the aperture mask used by the Kepler pipeline"""
        return self.hdu[-1].data > 2

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
    def cadenceno(self):
        """Return the cadence number for all good-quality cadences."""
        return self.hdu[1].data['CADENCENO'][self.quality_mask]

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
        """Returns the flux uncertainty for all good-quality cadences."""
        return self.hdu[1].data['FLUX_ERR'][self.quality_mask]

    @property
    def flux_bkg(self):
        """Returns the background flux for all good-quality cadences."""
        return self.hdu[1].data['FLUX_BKG'][self.quality_mask]

    @property
    def flux_bkg_err(self):
        return self.hdu[1].data['FLUX_BKG_ERR'][self.quality_mask]

    @property
    def quality(self):
        """Returns the quality flag integer of every good cadence."""
        return self.hdu[1].data['QUALITY'][self.quality_mask]

    @property
    def quarter(self):
        """Quarter number"""
        try:
            return self.header(ext=0)['QUARTER']
        except KeyError:
            return None

    @property
    def campaign(self):
        """Campaign number"""
        try:
            return self.header(ext=0)['CAMPAIGN']
        except KeyError:
            return None

    @property
    def mission(self):
        """Mission name"""
        return self.header(ext=0)['MISSION']

    def to_fits(self):
        """Save the TPF to fits"""
        raise NotImplementedError

    def to_lightcurve(self, aperture_mask=None):
        """Performs aperture photometry.

        Attributes
        ----------
        aperture_mask : array-like
            A boolean array describing the aperture such that `False` means
            that the pixel will be masked out.
            The default behaviour is to use all pixels.

        Returns
        -------
        lc : KeplerLightCurve object
            Array containing the summed flux within the aperture for each
            cadence.
        """
        if aperture_mask is None:
            mask = ~np.isnan(self.hdu[1].data['FLUX'][100])
            aperture_mask = np.ones((self.shape[1], self.shape[2]),
                                    dtype=bool) * mask

        centroid_col, centroid_row = self.centroids(aperture_mask)

        return KeplerLightCurve(flux=np.nansum(self.flux[:, aperture_mask], axis=1),
                                time=self.time, flux_err=self.flux_err,
                                centroid_col=centroid_col,
                                centroid_row=centroid_row,
                                quality=self.quality,
                                channel=self.channel,
                                campaign=self.campaign,
                                quarter=self.quarter,
                                mission=self.mission,
                                cadenceno=self.cadenceno)

    def centroids(self, aperture_mask=None):
        """Returns centroids based on sample moments.

        Parameters
        ----------
        aperture_mask : array-like or None
            A boolean array describing the aperture such that `False` means
            that the pixel will be masked out. The default behaviour is to
            use all pixels.

        Returns
        -------
        col_centr, row_centr : tuple
            Arrays containing centroids for column and row at each cadence
        """
        if aperture_mask is None:
            mask = ~np.isnan(self.hdu[1].data['FLUX'][100])
            aperture_mask = np.ones((self.shape[1], self.shape[2]),
                                    dtype=bool) * mask

        yy, xx = np.indices(self.shape[1:]) + 0.5
        yy = self.row + yy
        xx = self.column + xx
        total_flux = np.nansum(self.flux[:, aperture_mask], axis=1)
        col_centr = np.nansum(xx * aperture_mask * self.flux) / total_flux
        row_centr = np.nansum(yy * aperture_mask * self.flux) / total_flux

        return col_centr, row_centr

    def plot(self, frame=None, cadenceno=None, **kwargs):
        """
        Plot a target pixel file at a given frame (index) or cadence number.

        Parameters
        ----------
        frame : int
            Frame number.
        cadenceno : int
            Alternatively, a cadence number can be provided.
            This argument has priority over frame number.
        """
        if cadenceno is not None:
            frame = np.argwhere(cadenceno == self.cadenceno)[0][0]
        elif frame is None:
            raise ValueError("Either frame or cadenceno must be provided.")

        pflux = self.flux[frame]
        plot_image(pflux, title='Kepler ID: {}'.format(self.keplerid),
                   extent=(self.column, self.column + self.shape[2],
                           self.row, self.row + self.shape[1]), **kwargs)

    def get_bkg_lightcurve(self, aperture_mask=None):
        if aperture_mask is None:
            mask = self.hdu[1].data['FLUX'][100] == self.hdu[1].data['FLUX'][100]
            aperture_mask = np.ones((self.shape[1], self.shape[2]), dtype=bool) * mask

        return LightCurve(flux=np.nansum(self.flux_bkg[:, aperture_mask], axis=1),
                          time=self.time, flux_err=self.flux_bkg_err)
