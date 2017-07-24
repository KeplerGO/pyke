Quickstart
----------

If you have a working version of Python 2 or 3 on your system
(we recommend `Anaconda Python <https://www.continuum.io/downloads>`_),
you can simply install the latest stable release of PyKE using ``pip``::

    $ pip install pyketools

With PyKE installed, you can directly visualize frames from a target pixel file.
For example, let's visualize the pixels of Kepler target KIC008462852
(a.k.a. Taby's Star)::

    $ kepmask kplr008462852-2013098041711_lpd-targ.fits.gz --maskfile tabystar.txt

.. image:: _static/images/readme/kepmask.png

``kepmask`` is an interactive tool that allows one to create an arbitrary
aperture mask which can subsequently be used in another ``pyke`` tool,
such as ``kepextract``.

``kepextract`` performs simple aperture photometry in the pixels given by the mask created by ``kepmask``::

    $ kepextract kplr008462852-2013098041711_lpd-targ.fits.gz tabys-lc.fits --maskfile tabystar.txt

To visualize the light curve, we can use ``kepdraw``::

    $ kepdraw taby-lc.fits

.. image:: _static/images/readme/kepdraw.png
