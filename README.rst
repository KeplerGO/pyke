PyKE: Kepler, K2 & TESS Data Analysis Tools
============================================
|pypi-badge| |ci-badge| |appveyor-badge| |doc-badge| |cov-badge| |doi-badge|

.. |pypi-badge| image:: https://img.shields.io/pypi/v/pyketools.svg
                :target: https://pypi.python.org/pypi/pyketools
.. |ci-badge| image:: https://travis-ci.org/KeplerGO/PyKE.svg?branch=master
              :target: https://travis-ci.org/KeplerGO/PyKE
.. |appveyor-badge| image:: https://ci.appveyor.com/api/projects/status/6jvv5d7a142gwm8a/branch/master?svg=true
                    :target: https://ci.appveyor.com/project/mirca/pyke
.. |doc-badge| image:: https://readthedocs.org/projects/pyke/badge/?version=latest
              :target: https://pyke.keplerscience.org
.. |cov-badge| image:: https://codecov.io/gh/KeplerGO/PyKE/branch/master/graph/badge.svg
              :target: https://codecov.io/gh/KeplerGO/PyKE
.. |doi-badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.835584.svg
              :target: https://doi.org/10.5281/zenodo.835584


**Developed since 2012, PyKE offers a user-friendly way to inspect and analyze
the pixels and lightcurves obtained by NASA's Kepler, K2, and TESS missions.**

Documentation
-------------

Documentation is hosted at `pyke.keplerscience.org <http://pyke.keplerscience.org>`_.

What's new in PyKE v3? (July 2017)
----------------------------------


PyKE3 is the latest generation of the Kepler/K2/TESS toolkit.
It provides the following key improvements:

* PyKE3 is now a `pip-installable <http://pyke.keplerscience.org/en/latest/install.html#installing-pyke>`_ package and supports both Python 2 and 3
* `tasks <http://pyke.keplerscience.org/en/latest/overview.html>`_ are now available both as command-line tools and Python functions
* PyKE3 no longer depends on PyRAF and is TESS-ready
* PyKE3 address performance issues specially noticeable in short-cadence data
* documentation and tutorials are now generated using Sphinx
* the development has been moved to GitHub to encourage `user contributions <http://pyke.keplerscience.org/en/latest/contributing.html>`_

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

.. we should use full url addresses for images henceforth, so that they will be correctly captured by PYPI

.. image:: http://pyke.keplerscience.org/en/latest/_images/kepmask1.png

``kepmask`` is an interactive tool used to create a custom
aperture mask which can subsequently be used in other PyKE tasks.

For example, we can now use the ``kepextract`` task to perform aperture photometry using the pixels defined using ``kepmask`` above::

    $ kepextract kplr008462852-2013098041711_lpd-targ.fits.gz --outfile lightcurve.fits --maskfile mask.txt

This creates a file called ``lightcurve.fits`` which contains a lightcurve in a format similar to those found in the official archive.
To visualize the resulting light curve, we can use ``kepdraw``::

    $ kepdraw lightcurve.fits

.. image:: http://pyke.keplerscience.org/en/latest/_images/kepdraw1.png


Contributing
------------

Users are welcome to open `issues <https://github.com/KeplerGO/PyKE/issues>`_ or `pull requests <https://github.com/KeplerGO/PyKE/pulls>`_.
You can also contact the development team via keplergo@mail.arc.nasa.gov


Citing
------

If you find this code useful in your research, please cite both (Vinícius et al. 2017) and (Still & Barclay, 2012)
using the BibTeX provided below. Also, please give us a GitHub star!

::

    @misc{pyke3,
      author       = {Zé Vinícius and
                      Geert Barentsen and
                      Michael Gully-Santiago and
                      Ann Marie Cody and
                      Christina Hedges and
                      Martin Still and
                      Tom Barclay},
      title        = {KeplerGO/PyKE},
      month        = jul,
      year         = 2017,
      doi          = {10.5281/zenodo.835583},
      url          = {https://doi.org/10.5281/zenodo.835583}
    }

    @misc{2012ascl.soft08004S,
      author       = {{Still}, M. and {Barclay}, T.},
      title        = "{PyKE: Reduction and analysis of Kepler Simple Aperture Photometry data}",
      keywords     = {Software},
      howpublished = {Astrophysics Source Code Library},
      year         = 2012,
      month        = aug,
      archivePrefix= "ascl",
      eprint       = {1208.004},
      adsurl       = {http://adsabs.harvard.edu/abs/2012ascl.soft08004S}
    }
