3.1 (unreleased)
================

- [#66]  kepclip now supports Target Pixel Files
         kepclip no longer removes NaN values from fits files, users must run kepclean
         kepclean has been added, which removes NaN values from TPFs and LCs
- [#130] kepfold API updated to include optional keywords for period and BJD0.
         kepfold will otherwise search for period and BJD0 keywords in the input
         fits file.
         kepfold now includes the 'non-interactive' keyword
- [#133] kepfourier has been removed. All instances have been replaced with astropy's
         Lomb-Scargle method
         kepft has been renamed kepperiodogram
         kepft (kepperiodogram) now uses astropy's Lomb-Scargle to find periods
         kepft (kepperiodogram) now outputs the best fit period to a new fits extension
         kepdynamic has been sped up by using astropy LombScargle
         keptrial has been greatly sped up by using astropy LombScargle
         kepwindow has been sped up by using astropy LombScargle
         kepwindow api has been changed. Now includes the 'non-interactive' keyword.
         The nfreq and fmin keywords have been removed in favor of 'nyqfactor'.
- [#136] kepcotrend and keprange have new plot colors.
         The range of masked data is also plotted in kepcotrend.
         kepcotrend has had a minor bug fixed (a single cotrending basis vector can
         now be specified.)
         kepcotrend has had the keyword 'non-interactive' added (so plot windows
         can be suppressed.)
         kepcotrend now saves the plot in the same style as other routines
         (appending '.png' to the output file name.)
- [#138] a `compute_cotrend_lightcurve` method has been added in the KeplerLightCurveFile.


3.0 (2017-09-18)
================

PyKE 3.0.0 is a major release which provides the following key improvements:

- PyKE is now a [pip-installable](http://pyke.keplerscience.org/en/latest/install.html#installation)
  package and supports both Python 2 and 3.

- [Tasks](http://pyke.keplerscience.org/en/latest/overview.html#tasks) are now
  available both as command-line tools and Python functions.

- PyKE no longer depends on PyRAF.

- The tutorials have been updated for K2 and the [documentation](http://pyke.keplerscience.org)
  is now generated using Sphinx.

- The development has been moved to GitHub to encourage
  [user contributions](http://pyke.keplerscience.org/en/latest/contributing.html#contributing).

- The core code has been refactored to become more developer-friendly.
