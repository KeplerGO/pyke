..

Welcome to PyKE!
================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

PyKE is a set of data analysis tools which offer a user-friendly way to inspect and analyze the pixels and lightcurves obtained by NASA's Kepler, K2, and TESS missions.

Installation
============

The easiest way to install PyKE is through ``pip``::

    pip install pyketools

Or if you would like to experiment our development version::

    git clone https://github.com/KeplerGO/PyKE.git
    cd PyKE
    pip install -e .

Getting started
===============

.. toctree::
    :maxdepth: 2

    tutorials/index

Command-line tools
==================

After PyKE has been installed, a set of command-line tools is available on your
favorite terminal. Most of the tools can be run with the following pattern::

    $ name-of-the-tool input-file output-file --options

where ``input-file`` is either a light curve or a target pixel file and
``output-file`` is a name given by the user.

And help documentation can be retrived by::

    $ name-of-the-tool --help

See the list below for detailed documentation.

.. toctree::
    :maxdepth: 2

    api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Acknowledgment
--------------

PyKE is developed by the Kepler/K2 Guest Observer Office.
PyKE was firstly designed by Dr. Martin Still and Dr. Tom Barclay.
Check out the previous project webpage: https://keplerscience.arc.nasa.gov/PyKE.shtml
