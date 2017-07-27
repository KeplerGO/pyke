PyKE: Kepler, K2 & TESS Data Analysis Tools
============================================
|pypi-badge| |ci-badge| |doc-badge| |cov-badge| |bib-badge|

.. |pypi-badge| image:: https://img.shields.io/pypi/v/pyketools.svg
                :target: https://pypi.python.org/pypi/pyketools
.. |ci-badge| image:: https://travis-ci.org/KeplerGO/PyKE.svg?branch=master
              :target: https://travis-ci.org/KeplerGO/PyKE
.. |doc-badge| image:: https://readthedocs.org/projects/pyke/badge/?version=latest
              :target: https://pyke.readthedocs.io
.. |bib-badge| image:: https://img.shields.io/badge/NASA%20ADS-2012ascl.soft08004S-brightgreen.svg
              :target: http://adsabs.harvard.edu/abs/2012ascl.soft08004S
.. |cov-badge| image:: https://codecov.io/gh/KeplerGO/PyKE/branch/master/graph/badge.svg
              :target: https://codecov.io/gh/KeplerGO/PyKE


**Developed since 2012, PyKE offers a user-friendly way to inspect and analyze
the pixels and lightcurves obtained by NASA's Kepler, K2, and TESS missions.**

Documentation
-------------

Documentation is hosted at `pyke.keplerscience.org <https://pyke.keplerscience.org>`_.

What's new in PyKE v3? (July 2017)
----------------------------------


PyKE3 is the latest generation of the Kepler/K2/TESS toolkit.
It provides the following key improvements:

* PyKE is now a `pip-installable <http://pyke.readthedocs.io/en/latest/install.html#installing-pyke>`_ package and supports both Python 2 and 3;
* `tasks <http://pyke.readthedocs.io/en/latest/overview.html>`_ are now available both as command-line tools and Python functions;
* PyKE no longer depends on PyRAF and is TESS-ready.
* documentation and tutorials are now generated using Sphinx;
* the development has been moved to GitHub to encourage `user contributions <http://pyke.readthedocs.io/en/latest/contributing.html>`_.

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

.. image:: http://pyke.readthedocs.io/en/latest/_images/kepmask1.png

``kepmask`` is an interactive tool used to create a custom
aperture mask which can subsequently be used in other PyKE tasks.

For example, we can now use the ``kepextract`` task to perform aperture photometry using the pixels defined using ``kepmask`` above::

    $ kepextract kplr008462852-2013098041711_lpd-targ.fits.gz lightcurve.fits --maskfile mask.txt

This creates a file called ``lightcurve.fits`` which contains a lightcurve in a format similar to those found in the official archive.
To visualize the resulting light curve, we can use ``kepdraw``::

    $ kepdraw lightcurve.fits

.. image:: http://pyke.readthedocs.io/en/latest/_images/kepdraw1.png


Acknowledgement
---------------

If you find this code useful in your research, please consider `citing <http://adsabs.harvard.edu/abs/2012ascl.soft08004S>`_::

    Title: PyKE: Reduction and analysis of Kepler Simple Aperture Photometry data
    Authors: Still, Martin; Barclay, Tom
    Publication: Astrophysics Source Code Library, record ascl:1208.004
    Publication Date: 08/2012

*This package is developed by Martin Still, Tom Barclay, Ze Vinicius, Geert Barentsen, Michael Gully-Santiago, Ann Marie Cody, and Christine Hedges for the Kepler/K2 GO Office.*

Contributing
------------

Users are welcome to open `issues <https://github.com/KeplerGO/PyKE/issues>`_ or `pull requests <https://github.com/KeplerGO/PyKE/pulls>`_.
You can also contact the development team via keplergo@mail.arc.nasa.gov
