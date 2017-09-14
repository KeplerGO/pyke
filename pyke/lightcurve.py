import numpy as np

__all__ = ['LightCurve']

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

    def __init__(self, time, flux):
        self.time = time
        self.flux = flux

    def detrend(self, method='arclength', **kwargs):
        """
        """
        if method == 'arclength':
            return ArcLengthDetrender().detrend(self.time, self.flux, **kwargs)
        else:
            return FirstDifferenceDetrender().detrend(self.time, self.flux, **kwargs)

class Detrender(object):
    """
    """
    def detrend(self):
        """
        Returns a LightCurve object
        """
        pass

class SystematicsCorrector(object):
    """
    """
    def correct(**kwargs):
        """
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
