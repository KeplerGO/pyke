.. PyKE documentation master file, created by
   sphinx-quickstart on Thu May 11 16:52:07 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyKE!
================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

PyKE is a set of data analysis tools, maintained by the K2 Guest Observer
Office, for Kepler, K2, TESS, and future missions.

Installation
============

The easiest way to install PyKE is through ``pip``::

    pip install pyketools

Or if you would like to experiment our development version::

    git clone https://github.com/KeplerGO/PyKE.git
    cd PyKE
    pip install -e .

Command-line Tools Documentation
================================

After PyKE has been installed, a set of command-line tools is available on your
favorite terminal. Most of the tools can be run with the following pattern::

    $ name-of-the-tool input-file output-file --options

where ``input-file`` is either a light curve or a target pixel file and
``output-file`` is a name given by the user. See the list below for detailed
documentation.

.. toctree::
    :maxdepth: 2

    api/index

Tutorials
=========

.. toctree::
    :maxdepth: 2

    tutorials/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Acknowledgment
--------------

PyKE is developed by the Kepler/K2 Guest Observer Office.
PyKE was firstly designed by Dr. Martin Still and Dr. Tom Barclay.
Check out the previous project webpage: https://keplerscience.arc.nasa.gov/PyKE.shtml
