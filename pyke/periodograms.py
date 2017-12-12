import numpy as np
from gatspy import periodic

class Periodogram(object):
    '''
    Implements different periodograms when passed a lightcurve object'''

    def __init__(self, time, flux, flux_err=None):
        self.time = time
        self.flux = flux
        self.flux_err = flux_err

    def LombScargle(self,minper=0.05,maxper=100):
        model = periodic.LombScargleFast(fit_period=True)
        model.optimizer.quiet=True
        model.optimizer.period_range = (minper, maxper)
        model.fit(self.time, self.flux, self.flux.max()-self.flux.min())
        self.per=model.best_period

    def per(self):
        '''Returns the best fit period'''
        return self.per
