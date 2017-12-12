import numpy as np
from astropy.stats import LombScargle
import warnings
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('ggplot')
mpl.rc('axes', color_cycle=["#4C72B0", "#55A868", "#C44E52",
                            "#8172B2", "#CCB974"])


__all__ = ['Periodogram']

class Periodogram(object):
    '''
    Implements different periodograms when passed a lightcurve object'''

    def __init__(self, time, flux, flux_err=None, minper = None, maxper = None, nterms=1):
        self.time = time
        self.flux = flux
        if flux_err is None:
            flux_err = np.copy(flux)*0.
        self.flux_err = flux_err
        self.minper = minper
        self.maxper = maxper
        self.nterms = nterms
        if self.minper is None:
            self.minper=np.median(self.time[1:]-self.time[0:-1])*4
        if self.maxper is None:
            self.maxper=np.nanmax(self.time-self.time.min())/4.
        self.clean()
        self.LombScargle()


    def clean(self):
        ok = np.isfinite(self.flux)
        self.time=self.time[ok]
        self.flux=self.flux[ok]
        self.flux_err=self.flux_err[ok]

    def LombScargle(self,samples = 40):
        LS = LombScargle(self.time, self.flux, self.flux.max()-self.flux.min(),nterms=self.nterms)
        frequency, power = LS.autopower(maximum_frequency=1./self.minper,minimum_frequency=1./self.maxper,samples_per_peak=samples)
        self.lomb_per = 1./frequency[np.argmax(power)]
        self.lomb_periods = 1./frequency
        self.lomb_power = power

    def plot_LombScargle(self, ax = None, line=True,**kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.lomb_periods, self.lomb_power,**kwargs)
        ax.set(xlabel='Period (days)', ylabel='Lomb Scargle Power')
        if line:
            plt.axvline(self.lomb_per,ls='--',lw=2,color='black')
            ax.text(self.lomb_per*1.1,self.lomb_power.max(),'Best Fit Period', ha='left',va='center')

    def per(self):
        '''Returns the best fit period'''
        return self.lomb_per
