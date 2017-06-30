# PyKE3: Kepler, K2 & TESS Data Analysis Tools
[![Build Status](https://travis-ci.org/KeplerGO/PyKE.svg?branch=dev)](https://travis-ci.org/KeplerGO/PyKE) [![Documentation Status](https://readthedocs.org/projects/pyke/badge/?version=latest)](http://pyke.readthedocs.io/en/latest/?badge=latest) <a href="http://ascl.net/1208.004"><img src="https://img.shields.io/badge/ascl-1208.004-blue.svg?colorB=262255" alt="ascl:1208.004" /></a>

<!-- ***Tools for working with data from NASA's Kepler & TESS Space Telescopes.*** -->

Developed since 2012, PyKE offers a user-friendly way to inspect and analyze
the pixels and lightcurves obtained by NASA's Kepler, K2, and TESS missions.

## What's new in PyKE v3? (July 2017)

PyKE3 is the latest generation of the toolkit.
It provides the following key improvements:
* PyKE3 is now a true Python package installable via `pip`
* exposes true command-line tools
* support for both Python 2 and 3
* no dependency on PyRAF
* TESS-ready

## Example use

TBC

```
$ kepmask
```

```
$ kepextract
```


## Installation

If you have a working version of Python 2 or 3 on your system
(we recommend [Anaconda Python](https://www.continuum.io/downloads)),
you can simply install the latest stable release of PyKE using ``pip``:

    $ pip install pyketools

Alternatively, if you want to experiment with the latest development version of
PyKE, you can install it straight from GitHub:

    $ git clone https://github.com/KeplerGO/PyKE.git
    $ cd PyKE
    $ pip install -e .

Note: PyKE version 2 and older, which was in use between 2012 and 2016 and
required PyRAF, is available in the branch ``py27-pyraf``.

## Documentation

Documentation is hosted at [readthedocs](http://pyke.rtfd.io).


## Acknowledgement
If you find this code useful in your research, please consider [citing](http://adsabs.harvard.edu/abs/2012ascl.soft08004S):

```
Title:	PyKE: Reduction and analysis of Kepler Simple Aperture Photometry data
Authors:          Still, Martin; Barclay, Tom
Publication:      Astrophysics Source Code Library, record ascl:1208.004
Publication Date: 08/2012
```

*This package is developed by Martin Still, Tom Barclay, Ze Vinicius, Geert Barentsen, Michael Gully-Santiago, Ann Marie Cody, and Christine Hedges for the Kepler/K2 GO Office.*

## Contributing

Users are welcome to open [issues](https://github.com/KeplerGO/PyKE/issues) or [pull requests](https://github.com/KeplerGO/PyKE/pulls).
You can also contact the development team via keplergo@mail.arc.nasa.gov
