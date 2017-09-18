import numpy as np
from astropy.io import fits

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

    def __init__(self, time, flux, flux_err):
        self.time = time
        self.flux = flux
        self.flux_err = flux_err

    def detrend(self, method='arclength', **kwargs):
        """
        """
        if method == 'arclength':
            return ArcLengthDetrender().detrend(self.time, self.flux, **kwargs)
        else:
            return FirstDifferenceDetrender().detrend(self.time, self.flux, **kwargs)

    def draw(self):
        raise NotImplementedError("Should we implement a LightCurveDrawer class?")

    def to_csv(self):
        raise NotImplementedError()

    def to_fits(self):
        raise NotImplementedError()


class KeplerLightCurveFile(LightCurve):

    def __init__(self, path, **kwargs):
        self.hdu = fits.open(path, **kwargs)
        self.time = self.hdu[1].data['TIME']
        self.flux = self.get_flux('PDCSAP_FLUX')

    def flux_types(self):
        """Returns a list of available flux types for this light curve file"""
        return [n for n in  self.hdu[1].data.columns.names if 'FLUX' in n]

    def get_flux(self, flux_type):
        if flux_type is in self._flux_types():
            return self.hdu[1].data[flux_type]
        else:
            raise KeyError("{} is not a valid flux type. Available types are: {}".
                           format(flux_type, self._flux_types))

    def set_flux_type(self, flux_type):
        self.flux = self.get_flux(flux_type)



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
        return LightCurve(time, flux - np.append(0, flux[1:]))

class LinearDetrender(Detrender):
    """
    """
    @staticmethod
    def detrend(time, flux):
        pass

class ArcLengthDetrender(Detrender):
    def detrend(time, flux):
        pass

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
