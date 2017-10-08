import numpy as np
from astropy.io import fits
from scipy import linalg
from oktopus import L1Norm

__all__ = ['LightCurve', 'KeplerLightCurveFile']


class LightCurve(object):
    """
    Implements a basic time-series class for a generic lightcurve.

    Parameters
    ----------
    time : numpy array-like
        Time-line.
    flux : numpy array-like
        Data flux for every time point.
    """

    def __init__(self, time, flux, flux_err=None, quality=None, centroid_col=None,
                 centroid_row=None):
        self.time = time
        self.flux = flux
        self.flux_err = flux_err
        self.quality = quality
        self.centroid_col = centroid_col
        self.centroid_row = centroid_row

    def detrend(self, method='arclength', **kwargs):
        """
        """
        if method == 'arclength':
            return ArcLengthDetrender().detrend(time=self.time,
                                                flux=self.flux,
                                                flux_err=self.flux_err,
                                                centroid_col=self.centroid_col,
                                                centroid_row=self.centroid_row)
        else:
            return FirstDifferenceDetrender().detrend(time=self.time, flux=self.flux,
                                                      flux_err=self.flux_err, **kwargs)

    def draw(self):
        raise NotImplementedError("Should we implement a LightCurveDrawer class?")

    def to_csv(self):
        raise NotImplementedError()

    def to_fits(self):
        raise NotImplementedError()


class KeplerLightCurveFile(object):

    def __init__(self, path, **kwargs):
        self.hdu = fits.open(path, **kwargs)

    def get_lightcurve(self, flux_type, centroid_type='MOM_CENTR'):
        if flux_type in self._flux_types():
            return LightCurve(self.hdu[1].data['TIME'], self.hdu[1].data[flux_type],
                              flux_err=self.hdu[1].data[flux_type + "_ERR"],
                              quality=self.hdu[1].data['QUALITY'],
                              centroid_col=self.hdu[1].data[centroid_type + "1"],
                              centroid_row=self.hdu[1].data[centroid_type + "2"])
        else:
            raise KeyError("{} is not a valid flux type. Available types are: {}".
                           format(flux_type, self._flux_types))
    @property
    def SAP_FLUX(self):
        return self.get_lightcurve('SAP_FLUX')

    @property
    def PDCSAP_FLUX(self):
        return self.get_lightcurve('PDCSAP_FLUX')

    def _flux_types(self):
        """Returns a list of available flux types for this light curve file"""
        return [n for n in  self.hdu[1].data.columns.names if 'FLUX' in n]


class Detrender(object):
    """
    """
    def detrend(self):
        """
        Returns a LightCurve object
        """
        pass

class FirstDifferenceDetrender(Detrender):
    """
    First difference detrending
    """
    def detrend(time, flux):
        return LightCurve(time, np.append(flux[1:], 0) - flux)

class LinearDetrender(Detrender):
    """
    """
    @staticmethod
    def detrend(time, flux):
        pass

class ArcLengthDetrender(Detrender):
    def detrend(time, flux, centroid_col, centroid_row):
        centroids = self._rotate_centroids(centroid_col, centroid_row)

        # fit a polynomial to rotated centroids

        # compute the length of the curve of this polynomial

        #-----

        # fit B-splines to the light curve

        # divide raw light curve by B-splines

        # fit polynomial to flux as a function of arclength

        # divide the raw light curve by the polynomial fit

        # iterate between fitting a B-spline to the light curve and
        # a polynomial to flux as a function of arclength
        pass

    def _rotate_centroids(centroid_col, centroid_row):
        centroids = np.array([centroid_col, centroid_row])
        _, eig_vec = linalg.eigh(np.cov(centroids))
        return np.dot(eig_vec, centroids)

    def _remove_centroid_outliers(centroid_col, centroid_row):


class EMDDetrender(Detrender):
    """
    Empirical Mode Decomposition Detrender
    """
    def detrend(time, flux):
        pass

class PolynomialDetrender(Detrender):
    """
    """
    def detrend(time, flux):
        pass
