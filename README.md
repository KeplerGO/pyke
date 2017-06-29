# PyKE [![Build Status](https://travis-ci.org/KeplerGO/PyKE.svg?branch=dev)](https://travis-ci.org/KeplerGO/PyKE) [![Documentation Status](https://readthedocs.org/projects/pyke/badge/?version=latest)](http://pyke.readthedocs.io/en/latest/?badge=latest) <a href="http://ascl.net/1208.004"><img src="https://img.shields.io/badge/ascl-1208.004-blue.svg?colorB=262255" alt="ascl:1208.004" /></a>

***A suite of Python/PyRAF tools to analyze Kepler data.***

## Installation

The easiest way to install PyKE is through ``pip``:

    pip install pyketools

Or if you would like to experiment our development version:

    git clone https://github.com/KeplerGO/PyKE.git
    cd PyKE
    pip install -e .

Documentation is hosted at [readthedocs](http://pyke.rtfd.io).
Alternatively, it can be built locally:

    $ pip install sphinx numpydoc sphinx_rtd_theme
    $ cd pyke
    $ python setup.py build_sphinx

## Installation with PyRAF

To install PyKE within PyRAF, one may follow the instructions in the astroconda-iraf channel:
http://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf

After that, run the following commands on your favorite terminal:

1. ``mkiraf``
2. ``pyraf``
3. ``kepler``

The developer version of PyKE that is compatible with PyRAF is under the branch ``py27-pyraf``.


## Acknowledgement
If you find this code useful in your research, please consider [citing](http://adsabs.harvard.edu/abs/2012ascl.soft08004S):

```
Title:	PyKE: Reduction and analysis of Kepler Simple Aperture Photometry data
Authors:          Still, Martin; Barclay, Tom
Publication:      Astrophysics Source Code Library, record ascl:1208.004
Publication Date: 08/2012
```

*This package was mostly developed by [Tom Barclay](http://www.github.com/mrtommyb) and Martin Still.
Currently, this package is being developed and maintained by [Zé](http://www.github.com/mirca), [Gully](http://www.github.com/gully), and [Geert](http://www.github.com/barentsen).*


## Dependencies
```
numpy
scipy
astropy
matplotlib
tqdm
mdp (optional, needed for keppca)
```

## Support
Users are welcome to open [issues](https://github.com/KeplerGO/PyKE/issues) involving any aspects of this software
or submit [pull requests](https://github.com/KeplerGO/PyKE/pulls).

Feel free to contact us also through: keplergo@mail.arc.nasa.gov
