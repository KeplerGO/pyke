import numpy as np
from astropy.stats import LombScargle
import warnings
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib as mpl


__all__ = ['Periodogram']

class Periodogram(object):
    """Defines a periodogram class for Kepler/K2 data. Searches for periods using
    astropy.LombScargle and Box Least Squares.

    Attributes
    ----------
    time : array-like
        Time measurements
    flux : array-like
        Data flux for every time point
    flux_err : array-like
        Uncertainty on each flux data point
    minper : float
        Minimum period to search
    maxper : float
        Maximum period to search
    nterms : int
        Number of terms to use for Lomb-Scargle periodogram. (Default 1)
    """
    def __init__(self, time, flux, flux_err=None, minper=None, maxper=None, nterms=1):
        self.time = time
        self.flux = flux
        if flux_err is None:
            flux_err = np.copy(flux)*0.
        self.flux_err = flux_err
        self.minper = minper
        self.maxper = maxper
        self.nterms = nterms
        self.dur = time.max() - time.min()
        self.npoints = len(self.time)
        if self.minper is None:
            self.minper = np.median(self.time[1:] - self.time[0:-1])*4
        if self.maxper is None:
            self.maxper = np.nanmax(self.time - self.time.min())/4.
        self.clean()
        self.LombScargle()

    def clean(self):
        """
        Remove infinite values
        """
        ok = np.isfinite(self.flux)
        self.time = self.time[ok]
        self.flux = self.flux[ok]
        self.flux_err = self.flux_err[ok]

    def LombScargle(self,samples=40):
        """
        Creates a Lomb Scargle Periodogram

        Parameters
        ----------
        samples : int
            Number of samples to take in the
        """
        LS = LombScargle(self.time, self.flux, self.flux.max()-self.flux.min(), nterms=self.nterms)
        frequency, power = LS.autopower(maximum_frequency=1./self.minper, minimum_frequency=1./self.maxper, samples_per_peak=samples)
        self.lomb_per = 1. / frequency[np.argmax(power)]
        self.lomb_periods = 1. / frequency
        self.lomb_power = power

        s = np.argsort(self.time / self.lomb_per % 1)
        dp = self.npoints * (self.lomb_per * 5 / self.dur)
        if dp % 2 == 0:
            dp += 1
        smooth = signal.savgol_filter(self.flux[s], dp, 3)
        m = np.argmin(smooth)
        self.lomb_phase = (self.time[s] / period % 1)[m]

    def BLS(self):
        '''Not implemented'''
        pass

    def plot(self, ax=None, line=True, **kwargs):
        """
        Plot the periodogram object

        Parameters
        ----------
        ax : matplotlib frame
            Frame to plot the figure into. If unspecified, creates a new figure and frame.
        line : bool
            Plot a line at the bestfit period
        **kwargs : dict
            Dictionary of keyword values to pass to 'matplotlib.pyplot.plot'
        """
        if ax is None:
            fig, ax = plt.subplots()
        with mpl.style.use('ggplot'):
            ax.plot(self.lomb_periods, self.lomb_power,**kwargs)
            ax.set(xlabel='Period (days)', ylabel='Lomb Scargle Power')
            if line:
                plt.axvline(self.lomb_per,ls='--',lw=2,color='black')
                ax.text(self.lomb_per*1.1,self.lomb_power.max(),'Best Fit Period', ha='left',va='center')

    def per(self):
        '''Returns the best fit period.'''
        return self.lomb_per

    def phase(self):
        '''Returns the best fit phase.'''
        return self.lomb_phase
