import copy
import numpy as np
from scipy import linalg
from oktopus import L1Norm
from scipy import signal, interpolate
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from tqdm import tqdm
import oktopus
import requests
from bs4 import BeautifulSoup
from .utils import running_mean, channel_to_module_output, KeplerQualityFlags
from matplotlib import pyplot as plt

__all__ = ['LightCurve', 'KeplerLightCurve', 'KeplerLightCurveFile',
           'KeplerCBVCorrector', 'SPLDCorrector', 'SFFCorrector']


class LightCurve(object):
    """
    Implements a simple class for a generic light curve.

    Attributes
    ----------
    time : array-like
        Time measurements
    flux : array-like
        Data flux for every time point
    flux_err : array-like
        Uncertainty on each flux data point
    """

    def __init__(self, time, flux, flux_err=None):
        self.time = np.asarray(time)
        self.flux = np.asarray(flux)
        if flux_err is not None:
            self.flux_err = np.asarray(flux_err)
        else:
            self.flux_err = None

    def stitch(self, *others):
        """
        Stitches LightCurve objects.

        Parameters
        ----------
        *others : LightCurve objects
            Light curves to be stitched.

        Returns
        -------
        stitched_lc : LightCurve object
            Stitched light curve.
        """
        time = self.time
        flux = self.flux
        flux_err = self.flux_err

        for i in range(len(others)):
            time = np.append(time, others[i].time)
            flux = np.append(flux, others[i].flux)
            flux_err = np.append(flux_err, others[i].flux_err)

        return LightCurve(time=time, flux=flux, flux_err=flux_err)

    def flatten(self, window_length=101, polyorder=3, return_trend=False, **kwargs):
        """
        Removes low frequency trend using scipy's Savitzky-Golay filter.

        Parameters
        ----------
        window_length : int
            The length of the filter window (i.e. the number of coefficients).
            ``window_length`` must be a positive odd integer.
        polyorder : int
            The order of the polynomial used to fit the samples. ``polyorder``
            must be less than window_length.
        return_trend : bool
            If `True`, the method will return a tuple of two elements
            (flattened_lc, trend_lc) where trend_lc is the removed trend.
        **kwargs : dict
            Dictionary of arguments to be passed to `scipy.signal.savgol_filter`.

        Returns
        -------
        flatten_lc : LightCurve object
            Flattened lightcurve.
        If `return_trend` is `True`, the method will also return:
        trend_lc : LightCurve object
            Trend in the lightcurve data
        """
        lc_clean = self.remove_nans()  # The SG filter does not allow NaNs
        trend_signal = signal.savgol_filter(x=lc_clean.flux,
                                            window_length=window_length,
                                            polyorder=polyorder, **kwargs)
        flatten_lc = copy.copy(lc_clean)
        flatten_lc.flux = lc_clean.flux / trend_signal
        if flatten_lc.flux_err is not None:
            flatten_lc.flux_err = lc_clean.flux_err / trend_signal

        if return_trend:
            trend_lc = copy.copy(self)
            trend_lc.flux = trend_signal
            return flatten_lc, trend_lc
        return flatten_lc

    def fold(self, period, phase=0.):
        """Folds the lightcurve at a specified ``period`` and ``phase``.

        This method returns a new ``LightCurve`` object in which the time
        values range between -0.5 to +0.5.  Data points which occur exactly
        at ``phase`` or an integer multiple of `phase + n*period` have time
        value 0.0.

        Parameters
        ----------
        period : float
            The period upon which to fold.
        phase : float, optional
            Time reference point.

        Returns
        -------
        folded_lightcurve : LightCurve object
            A new ``LightCurve`` in which the data are folded and sorted by
            phase.
        """
        fold_time = ((self.time - phase + 0.5 * period) / period) % 1 - 0.5
        sorted_args = np.argsort(fold_time)
        if self.flux_err is None:
            return LightCurve(fold_time[sorted_args], self.flux[sorted_args])
        return LightCurve(fold_time[sorted_args], self.flux[sorted_args], flux_err=self.flux_err[sorted_args])

    def remove_nans(self):
        """Removes cadences where the flux is NaN.

        Returns
        -------
        clean_lightcurve : LightCurve object
            A new ``LightCurve`` from which NaNs fluxes have been removed.
        """
        lc = copy.copy(self)
        nan_mask = np.isnan(lc.flux)
        lc.time = self.time[~nan_mask]
        lc.flux = self.flux[~nan_mask]
        if lc.flux_err is not None:
            lc.flux_err = self.flux_err[~nan_mask]
        return lc

    def remove_outliers(self, sigma=5., **kwargs):
        """Removes outlier flux values using sigma-clipping.

        This method returns a new LightCurve object from which flux values
        are removed if they are separated from the mean flux by `sigma` times
        the standard deviation.

        Parameters
        ----------
        sigma : float, optional
            The number of standard deviations to use for clipping outliers.
            Defaults to 5.
        **kwargs : dict
            Dictionary of arguments to be passed to `astropy.stats.sigma_clip`.

        Returns
        -------
        clean_lightcurve : LightCurve object
            A new ``LightCurve`` in which outliers have been removed.
        """
        new_lc = copy.copy(self)
        outlier_mask = sigma_clip(data=new_lc.flux, sigma=sigma, **kwargs).mask
        new_lc.time = self.time[~outlier_mask]
        new_lc.flux = self.flux[~outlier_mask]
        if new_lc.flux_err is not None:
            new_lc.flux_err = self.flux_err[~outlier_mask]
        return new_lc

    def cdpp(self, transit_duration=13, savgol_window=101, savgol_polyorder=2,
             sigma_clip=5.):
        """Estimate the CDPP noise metric using the Savitzky-Golay (SG) method.

        A common estimate of the noise in a lightcurve is the scatter that
        remains after all long term trends have been removed. This is the idea
        behind the Combined Differential Photometric Precision (CDPP) metric.
        The official Kepler Pipeline computes this metric using a wavelet-based
        algorithm to calculate the signal-to-noise of the specific waveform of
        transits of various durations. In this implementation, we use the
        simpler "sgCDPP proxy algorithm" discussed by Gilliland et al
        (2011ApJS..197....6G) and Van Cleve et al (2016PASP..128g5002V).

        The steps of this algorithm are:
            1. Remove low frequency signals using a Savitzky-Golay filter with
               window length `savgol_window` and polynomial order `savgol_polyorder`.
            2. Remove outliers by rejecting data points which are separated from
               the mean by `sigma_clip` times the standard deviation.
            3. Compute the standard deviation of a running mean with
               a configurable window length equal to `transit_duration`.

        We use a running mean (as opposed to block averaging) to strongly
        attenuate the signal above 1/transit_duration whilst retaining
        the original frequency sampling.  Block averaging would set the Nyquist
        limit to 1/transit_duration.

        Parameters
        ----------
        transit_duration : int, optional
            The transit duration in cadences. This is the length of the window
            used to compute the running mean. The default is 13, which
            corresponds to a 6.5 hour transit in data sampled at 30-min cadence.
        savgol_window : int, optional
            Width of Savitsky-Golay filter in cadences (odd number).
            Default value 101 (2.0 days in Kepler Long Cadence mode).
        savgol_polyorder : int, optional
            Polynomial order of the Savitsky-Golay filter.
            The recommended value is 2.
        sigma_clip : float, optional
            The number of standard deviations to use for clipping outliers.
            The default is 5.

        Returns
        -------
        cdpp : float
            Savitzky-Golay CDPP noise metric in units parts-per-million (ppm).

        Notes
        -----
        This implementation is adapted from the Matlab version used by
        Jeff van Cleve but lacks the normalization factor used there:
        svn+ssh://murzim/repo/so/trunk/Develop/jvc/common/compute_SG_noise.m
        """
        if not isinstance(transit_duration, int):
            raise TypeError("transit_duration must be an integer")
        detrended_lc = self.flatten(window_length=savgol_window,
                                    polyorder=savgol_polyorder)
        cleaned_lc = detrended_lc.remove_outliers(sigma=sigma_clip)
        mean = running_mean(data=cleaned_lc.flux, window_size=transit_duration)
        cdpp_ppm = np.std(mean) * 1e6
        return cdpp_ppm

    def to_csv(self):
        raise NotImplementedError()

    def plot(self, ax=None, normalize=True, xlabel='Time - 2454833 (days)',
             ylabel='Normalized Flux', title=None, color='#363636', linestyle="",
             fill=False, grid=True, **kwargs):
        """Plots the light curve.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            A matplotlib axes object to plot into. If no axes is provided,
            a new one be generated.
        normalize : bool
            Normalized the lightcurve
        xlabel : str
            Plot x axis label
        ylabel : str
            Plot y axis label
        title : str
            Plot set_title
        color: str
            Color to plot flux points
        fill: bool
            Shade the region between 0 and flux
        grid: bool
            Add a grid to the plot
        **kwargs : dict
            Dictionary of arguments to be passed to `matplotlib.pyplot.plot`.

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            The matplotlib axes object.
        """
        if ax is None:
            fig, ax = plt.subplots(1)
        flux = self.flux
        flux_err = self.flux_err
        if normalize:
            if flux_err is not None:
                flux_err = flux_err / np.nanmedian(flux)
            flux = flux / np.nanmedian(flux)
        if flux_err is None:
            ax.plot(self.time, flux, marker='o', color=color, linestyle=linestyle,
                       **kwargs)
        else:
            ax.errorbar(self.time, flux, flux_err, color=color, linestyle=linestyle,
                        **kwargs)
        if fill:
            ax.fill(self.time, flux, fc='#a8a7a7', linewidth=0.0, alpha=0.3)
        if grid:
            ax.grid(alpha=0.3)
        if 'label' in kwargs:
            ax.legend()
        if title is not None:
            ax.set_title(title)
        ax.set_xlabel(xlabel, {'color': 'k'})
        ax.set_ylabel(ylabel, {'color': 'k'})
        return ax


class KeplerLightCurve(LightCurve):
    """Defines a light curve class for NASA's Kepler and K2 missions.

    Attributes
    ----------
    time : array-like
        Time measurements
    flux : array-like
        Data flux for every time point
    flux_err : array-like
        Uncertainty on each flux data point
    centroid_col, centroid_row : array-like, array-like
        Centroid column and row coordinates as a function of time
    quality : array-like
        Array indicating the quality of each data point
    quality_bitmask : int
        Bitmask specifying quality flags of cadences that should be ignored
    channel : int
        Channel number
    campaign : int
        Campaign number
    quarter : int
        Quarter number
    mission : str
        Mission name
    cadenceno : array-like
        Cadence numbers corresponding to every time measurement
    keplerid : int
        Kepler ID number
    """

    def __init__(self, time, flux, flux_err=None, centroid_col=None,
                 centroid_row=None, quality=None, quality_bitmask=None,
                 channel=None, campaign=None, quarter=None, mission=None,
                 cadenceno=None, keplerid=None):
        super(KeplerLightCurve, self).__init__(time, flux, flux_err)
        self.centroid_col = centroid_col
        self.centroid_row = centroid_row
        self.quality = quality
        self.quality_bitmask = quality_bitmask
        self.channel = channel
        self.campaign = campaign
        self.quarter = quarter
        self.mission = mission
        self.cadenceno = cadenceno
        self.keplerid = keplerid

    def correct(self, method='vanderburg', **kwargs):
        """Corrects a lightcurve for motion-dependent systematical errors.

        Parameters
        ----------
        method : str
        kwargs : dict

        Returns
        -------
        new_lc : KeplerLightCurve object
            Corrected lightcurve
        """
        if method == 'vanderburg':
            self.corrector = SFFCorrector(time=self.time, flux=self.flux,
                                          centroid_col=self.centroid_col,
                                          centroid_row=self.centroid_row)
            corrected_lc = self.corrector.correct(**kwargs)
        else:
            raise ValueError("method {} is not available.".format(method))
        new_lc = copy.copy(self)
        new_lc.flux = corrected_lc.flux
        return new_lc

    def to_fits(self):
        raise NotImplementedError()


class KeplerLightCurveFile(object):
    """Defines a class for a given light curve FITS file from NASA's Kepler and
    K2 missions.

    Attributes
    ----------
    path : str
        Directory path or url to a lightcurve FITS file.
    quality_bitmask : str or int
        Bitmask specifying quality flags of cadences that should be ignored.
        If a string is passed, it has the following meaning:

            * default: recommended quality mask
            * hard: removes more flags, known to remove good data
            * hardest: removes all data that has been flagged
    kwargs : dict
        Keyword arguments to be passed to astropy.io.fits.open.
    """
    def __init__(self, path, quality_bitmask=KeplerQualityFlags.DEFAULT_BITMASK,
                 **kwargs):
        self.path = path
        self.hdu = pyfits.open(self.path, **kwargs)
        self.quality_bitmask = quality_bitmask
        self.quality_mask = self._quality_mask(quality_bitmask)

    def get_lightcurve(self, flux_type, centroid_type='MOM_CENTR'):
        if flux_type in self._flux_types():
            return KeplerLightCurve(self.hdu[1].data['TIME'][self.quality_mask],
                                    self.hdu[1].data[flux_type][self.quality_mask],
                                    flux_err=self.hdu[1].data[flux_type + "_ERR"][self.quality_mask],
                                    centroid_col=self.hdu[1].data[centroid_type + "1"][self.quality_mask],
                                    centroid_row=self.hdu[1].data[centroid_type + "2"][self.quality_mask],
                                    quality=self.hdu[1].data['SAP_QUALITY'][self.quality_mask],
                                    quality_bitmask=self.quality_bitmask,
                                    channel=self.channel,
                                    campaign=self.campaign,
                                    quarter=self.quarter,
                                    mission=self.mission,
                                    cadenceno=self.cadenceno,
                                    keplerid=self.hdu[0].header['KEPLERID'])
        else:
            raise KeyError("{} is not a valid flux type. Available types are: {}".
                           format(flux_type, self._flux_types))

    def _quality_mask(self, bitmask):
        """Returns a boolean mask which flags all good-quality cadences.

        Parameters
        ----------
        bitmask : str or int
            Bitmask. See ref. [1], table 2-3.

        Returns
        -------
        boolean_mask : array of bool
            Boolean array in which `True` means the data is of good quality.
        """
        if bitmask is None:
            return np.ones(len(self.hdu[1].data['TIME']), dtype=bool)
        elif isinstance(bitmask, str):
            bitmask = KeplerQualityFlags.OPTIONS[bitmask]
        return (self.hdu[1].data['SAP_QUALITY'] & bitmask) == 0

    @property
    def SAP_FLUX(self):
        """Returns a KeplerLightCurve object for SAP_FLUX"""
        return self.get_lightcurve('SAP_FLUX')

    @property
    def PDCSAP_FLUX(self):
        """Returns a KeplerLightCurve object for PDCSAP_FLUX"""
        return self.get_lightcurve('PDCSAP_FLUX')

    @property
    def time(self):
        """Time measurements"""
        return self.hdu[1].data['TIME'][self.quality_mask]

    @property
    def cadenceno(self):
        """Cadence number"""
        return self.hdu[1].data['CADENCENO'][self.quality_mask]

    @property
    def channel(self):
        """Channel number"""
        return self.header(ext=0)['CHANNEL']

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

    def compute_cotrended_lightcurve(self, cbvs=[1, 2]):
        """Returns a LightCurve object after cotrending the SAP_FLUX
        against the cotrending basis vectors.

        Parameters
        ----------
        cbvs : list of ints
            The list of cotrending basis vectors to fit to the data. For example,
            [1, 2] will fit the first two basis vectors.
        """
        return KeplerCBVCorrector(self).correct(cbvs=cbvs)

    def header(self, ext=0):
        """Header of the object at extension `ext`"""
        return self.hdu[ext].header

    def _flux_types(self):
        """Returns a list of available flux types for this light curve file"""
        types = [n for n in self.hdu[1].data.columns.names if 'FLUX' in n]
        types = [n for n in types if not ('ERR' in n)]
        return types

    def plot(self, plottype=None, **kwargs):
        """Plot all the flux types in a light curve.

        Parameters
        ----------
        plottype : str or list of str
            List of FLUX types to plot. Default is to plot all available.
        """
        if not ('ax' in kwargs):
            fig, ax = plt.subplots(1)
            kwargs['ax'] = ax
        if not ('title' in kwargs):
            kwargs['title'] = 'KeplerID: {}'.format(self.SAP_FLUX.keplerid)
        if plottype is None:
            plottype = self._flux_types()
        if isinstance(plottype, str):
            plottype = [plottype]
        for idx, pl in enumerate(plottype):
            lc = self.get_lightcurve(pl)
            kwargs['color'] = 'C{}'.format(idx)
            lc.plot(label=pl, **kwargs)


class SFFCorrector(object):
    """Implements a similar procedure as described by Vanderburg and Johnson
    (2014). Briefly, the algorithm implemented in this class can be described
    as follows

       (1) Rotate the centroid measurements onto the subspace spanned by the
           eigenvectors of the centroid covariance matrix
       (2) Fit a polynomial to the rotated centroids
       (3) Compute the arclength of such polynomial
       (4) Fit a BSpline of the raw flux as a function of time
       (5) Normalize the raw flux by the fitted BSpline computed in step (4)
       (6) Bin and interpolate the normalized flux as function of the arclength
       (7) Divide the raw flux by the piecewise linear interpolation done in step [(6)
       (8) Set raw flux as the flux computed in step (7) and repeat

    Attributes
    ----------
    time : array-like
        Time measurements
    flux : array-like
        Data flux for every time point
    centroid_col, centroid_row : array-like, array-like
        Centroid column and row coordinates as a function of time
    """

    def __init__(self, time, flux, centroid_col, centroid_row):
        self.time = time
        self.flux = flux
        self.centroid_col = centroid_col
        self.centroid_row = centroid_row

    def correct(self, polyorder=5, niters=3, bins=15, windows=1):
        """
        Parameters
        ----------
        polyorder : int
            Degree of the polynomial which will be used to fit one
            centroid as a function of the other.
        niters : int
            Number of iterations of the aforementioned algorithm.
        bins : int
            Number of bins to be used in step (6).
        windows : int
            Number of windows to subdivide the data.
        """
        time = np.array_split(self.time, windows)
        flux = np.array_split(self.flux, windows)
        centroid_col = np.array_split(self.centroid_col, windows)
        centroid_row = np.array_split(self.centroid_row, windows)

        flux_hat = np.array([])
        for i in range(windows):
            # Rotate and fit centroids
            self.rot_col, self.rot_row = self.rotate_centroids(centroid_col[i], centroid_row[i])
            coeffs = np.polyfit(self.rot_row, self.rot_col, polyorder)
            self.poly = np.poly1d(coeffs)
            self.polyprime = np.poly1d(coeffs).deriv()

            # Compute the arclength s
            x = np.linspace(np.min(self.rot_row), np.max(self.rot_row), 10000)
            self.s = np.array([self.arclength(x1=xp, x=x) for xp in self.rot_row])

            for n in range(niters):
                # fit BSpline
                self.bspline = self.fit_bspline(time[i], flux[i])
                # Normalize raw flux
                self.normflux = flux[i] / self.bspline(time[i] - time[i][0])
                # Bin and interpolate normalized flux
                self.interp = self.bin_and_interpolate(self.s, self.normflux, bins)
                # Correcte the raw flux
                corrected_flux = self.normflux / self.interp(self.s)
                flux[i] = corrected_flux

            flux_hat = np.append(flux_hat, flux[i])

        return LightCurve(time=self.time, flux=flux_hat)

    def rotate_centroids(self, centroid_col, centroid_row):
        centroids = np.array([centroid_col, centroid_row])
        _, eig_vecs = linalg.eigh(np.cov(centroids))
        return np.dot(eig_vecs, centroids)

    def _plot_rotated_centroids(self):
        plt.plot(self.rot_row, self.rot_col, 'ko', markersize=3)
        plt.plot(self.rot_row, self.rot_col, 'ro', markersize=2)
        x = np.linspace(min(self.rot_row), max(self.rot_row), 200)
        plt.plot(x, self.poly(x), '--')

    def _plot_normflux_arclength(self):
        plt.plot(self.s, self.normflux, 'ko', markersize=3)
        plt.plot(self.s, self.normflux, 'bo', markersize=2)
        ss = np.sort(self.s)
        plt.plot(ss, self.interp(ss), '--')

    def arclength(self, x1, x):
        """Compute the arclength of the polynomial used to fit the centroid
        measurements.

        Parameters
        ----------
        x1 : float
            Upper limit of the integration domain.
        x : ndarray
            Domain at which the arclength integrand is defined.

        Returns
        -------
        arclength : float
            Result of the integral from x[0] to x1.
        """
        mask = x < x1
        return np.trapz(y=np.sqrt(1 + self.polyprime(x[mask]) ** 2), x=x[mask])

    def fit_bspline(self, time, flux, s=0):
        time = time - time[0]
        knots = np.arange(0, time[-1], 1.5)
        t, c, k = interpolate.splrep(time, flux, t=knots[1:], s=s, task=-1)
        return interpolate.BSpline(t, c, k)

    def bin_and_interpolate(self, s, normflux, bins):
        idx = np.argsort(s)
        knots = np.array([np.min(s)]
                         + [np.median(split) for split in np.array_split(s[idx], bins)]
                         + [np.max(s)])
        bin_means = np.array([normflux[idx][0]]
                             + [np.mean(split) for split in np.array_split(normflux[idx], bins)]
                             + [normflux[idx][-1]])
        return interpolate.interp1d(knots, bin_means)

    def breakpoints(self, campaign):
        """Return a break point as a function of the campaign number.
        """
        pass


class KeplerCBVCorrector(object):
    r"""Remove systematic trends from Kepler light curves by fitting
    cotrending basis vectors.

    .. math::

         \arg \min_{\theta \in \Theta} |f(t) - <\theta, \left[\rm{cbv}_1(t), ..., \rm{cbv}_n(t)\right]^{T}|^p

    Attributes
    ----------
    lc_file : KeplerLightCurveFile object or str
        An instance from KeplerLightCurveFile or a path for the .fits
        file of a NASA's Kepler/K2 light curve.
    loss_function : oktopus.Likelihood subclass
        A class that describes a cost function.
        The default is :class:`oktopus.LaplacianLikelihood`, which is tantamount
        to the L1 norm.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from pyke import KeplerCBVCorrector, KeplerLightCurveFile
    >>> fn = ("https://archive.stsci.edu/missions/kepler/lightcurves/"
    ...       "0084/008462852/kplr008462852-2011073133259_llc.fits")
    >>> cbv = KeplerCBVCorrector(fn)
    Downloading https://archive.stsci.edu/missions/kepler/lightcurves/0084/008462852/kplr008462852-2011073133259_llc.fits [Done]
    >>> cbv_lc = cbv.correct()
    Downloading http://archive.stsci.edu/missions/kepler/cbv/kplr2011073133259-q08-d25_lcbv.fits [Done]
    >>> sap_lc = KeplerLightCurveFile(fn).SAP_FLUX
    >>> plt.plot(sap_lc.time, sap_lc.flux, 'x', markersize=1, label='SAP_FLUX') # doctest: +SKIP
    >>> plt.plot(cbv_lc.time, cbv_lc.flux, 'o', markersize=1, label='CBV_FLUX') # doctest: +SKIP
    >>> plt.legend() # doctest: +SKIP
    """

    def __init__(self, lc_file, loss_function=oktopus.LaplacianLikelihood):
        self.lc_file = lc_file
        self.loss_function = loss_function

        if self.lc_file.mission == 'Kepler':
            self.cbv_base_url = "http://archive.stsci.edu/missions/kepler/cbv/"
        elif self.lc_file.mission == 'K2':
            self.cbv_base_url = "http://archive.stsci.edu/missions/k2/cbv/"

    @property
    def lc_file(self):
        return self._lc_file

    @lc_file.setter
    def lc_file(self, value):
        # this enables `lc_file` to be either a string
        # or an object from KeplerLightCurveFile
        if isinstance(value, str):
            self._lc_file = KeplerLightCurveFile(value)
        elif isinstance(value, KeplerLightCurveFile):
            self._lc_file = value
        else:
            raise ValueError("lc_file must be either a string or a"
                             " KeplerLightCurveFile instance, got {}.".format(value))

    @property
    def coeffs(self):
        """
        Returns the fitted coefficients.
        """
        return self._coeffs

    @property
    def opt_result(self):
        """
        Returns the result of the optimization process.
        """
        return self._opt_result

    def correct(self, cbvs=[1, 2]):
        """
        Correct the SAP_FLUX by fitting a number of cotrending basis vectors
        `cbvs`.

        Parameters
        ----------
        cbvs : list of ints
            The list of cotrending basis vectors to fit to the data. For example,
            [1, 2] will fit the first two basis vectors.
        """
        module, output = channel_to_module_output(self.lc_file.channel)
        cbv_file = pyfits.open(self.get_cbv_url())
        cbv_data = cbv_file['MODOUT_{0}_{1}'.format(module, output)].data

        cbv_array = []
        for i in cbvs:
            cbv_array.append(cbv_data.field('VECTOR_{}'.format(i))[self.lc_file.quality_mask])
        cbv_array = np.asarray(cbv_array)

        sap_lc = self.lc_file.SAP_FLUX
        median_sap_flux = np.nanmedian(sap_lc.flux)
        norm_sap_flux = sap_lc.flux / median_sap_flux - 1
        norm_err_sap_flux = sap_lc.flux_err / median_sap_flux

        def mean_model(*theta):
            coeffs = np.asarray(theta)
            return np.dot(coeffs, cbv_array)

        loss = self.loss_function(data=norm_sap_flux, mean=mean_model,
                                  var=norm_err_sap_flux)
        self._opt_result = loss.fit(x0=np.zeros(len(cbvs)), method='L-BFGS-B')
        self._coeffs = self._opt_result.x
        flux_hat = sap_lc.flux - median_sap_flux * mean_model(self._coeffs)

        return LightCurve(time=sap_lc.time, flux=flux_hat.reshape(-1))

    def get_cbv_url(self):
        # gets the html page and finds all references to 'a' tag
        # keeps the ones for which 'href' ends with 'fits'
        # this might slow things down in case the user wants to fit 1e3 stars
        soup = BeautifulSoup(requests.get(self.cbv_base_url).text, 'html.parser')
        cbv_files = [fn['href'] for fn in soup.find_all('a') if fn['href'].endswith('fits')]

        if self.lc_file.mission == 'Kepler':
            if self.lc_file.quarter < 10:
                quarter = 'q0' + str(self.lc_file.quarter)
            else:
                quarter = 'q' + str(self.lc_file.quarter)
            for cbv_file in cbv_files:
                if quarter + '-d25' in cbv_file:
                    break
        elif self.lc_file.mission == 'K2':
            if self.lc_file.campaign <= 8:
                campaign = 'c0' + str(self.lc_file.campaign)
            else:
                campaign = 'c' + str(self.lc_file.campaign)
            for cbv_file in cbv_files:
                if campaign in cbv_file:
                    break

        return self.cbv_base_url + cbv_file


class SPLDCorrector(object):
    r"""
    Implements the simple first order Pixel Level Decorrelation (PLD) proposed by
    Deming et. al. [1]_ and Luger et. al. [2]_, [3]_.

    Attributes
    ----------
    time : array-like
        Time array
    tpf_flux : array-like
        Pixel values series

    Notes
    -----
    This code serves only as a quick look into the PLD technique.
    Users are encouraged to check out the GitHub repos
    `everest <http://www.github.com/rodluger/everest>`_
    and `everest3 <http://www.github.com/rodluger/everest3>`_.

    References
    ----------
    .. [1] Deming et. al. Spitzer Secondary Eclipses of the Dense, \
           Modestly-irradiated, Giant Exoplanet HAT-P-20b using Pixel-Level Decorrelation.
    .. [2] Luger et. al. EVEREST: Pixel Level Decorrelation of K2 Light Curves.
    .. [3] Luger et. al. An Update to the EVEREST K2 Pipeline: short cadence, \
           saturated stars, and Kepler-like photometry down to K_p = 15.
    """

    def __init__(self, time, tpf_flux):
        self.time = time
        self.tpf_flux = tpf_flux

    def correct(self, window_length=None, polyorder=2):
        k = window_length
        if not k:
            k = int(len(self.time) / 2) - 1
        n_windows = int(len(self.time) / k)
        flux_hat = np.array([])
        for n in range(1, n_windows + 1):
            flux_hat = np.append(flux_hat,
                                 self._pld(self.tpf_flux[(n - 1) * k:n * k], polyorder))
        flux_hat = np.append(flux_hat, self._pld(self.tpf_flux[n * k:], polyorder))
        return LightCurve(self.time, flux_hat + np.nanmedian(np.nansum(self.tpf_flux, axis=(1, 2))))

    def _pld(self, tpf_flux, polyorder=2):
        if len(tpf_flux) == 0:
            return np.array([])
        pixels_series = tpf_flux.reshape((tpf_flux.shape[0], -1))
        lightcurve = np.nansum(pixels_series, axis=1).reshape(-1, 1)
        # design matrix
        X = pixels_series / lightcurve
        X = np.hstack((X, np.array([np.linspace(0, 1, tpf_flux.shape[0]) ** n for n in range(polyorder+1)]).T))
        opt_weights = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, lightcurve))
        model = np.dot(X, opt_weights)
        flux_hat = lightcurve - model
        return flux_hat
