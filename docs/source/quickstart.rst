Quickstart
----------

If you have a working version of Python 2 or 3 on your system
(we recommend `Anaconda Python <https://www.continuum.io/downloads>`_),
you can simply install the latest stable release of PyKE using ``pip``::

    $ pip install pyketools

With PyKE installed, you can directly visualize frames from a target pixel file.
For example, let's visualize the pixels of Kepler target KIC008462852
(a.k.a. Tabby's Star)::

    $ kepmask kplr008462852-2013098041711_lpd-targ.fits.gz --maskfile mask.txt

.. image:: _static/images/readme/kepmask.png

``kepmask`` is an interactive tool used to create a custom
aperture mask (by clicking on the desired pixels and hitting DUMP)
which can subsequently be used in other PyKE tasks.

For example, we can now use the ``kepextract`` task to perform aperture photometry using the pixels defined using ``kepmask`` above::

    $ kepextract kplr008462852-2013098041711_lpd-targ.fits.gz --maskfile mask.txt

This creates a file called ``kplr008462852-2013098041711_lpd-targ-kepextract.fits`` which contains a lightcurve in a format similar to those found in the official archive.
To visualize the resulting light curve, we can use ``kepdraw``::

    $ kepdraw kplr008462852-2013098041711_lpd-targ-kepextract.fits

.. image:: _static/images/readme/kepdraw.png
