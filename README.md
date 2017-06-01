# PyKE
***A suite of Python/PyRAF tools to analyze Kepler data.***

For more information and documentation,
visit http://keplerscience.arc.nasa.gov/software.html#pyke


## Installation

The easiest way to install PyKE is through ``pip``::

    pip install pyke

Or if you would like to experiment our development version::

    git clone https://github.com/KeplerGO/PyKE.git
    cd PyKE
    pip install -e .


## Installation with PyRAF

The easiest way to install PyKE within PyRAF is through the astroconda-iraf channel:
http://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf

After that, run the following commands on a terminal of your preference:

1. ``mkiraf``
2. ``pyraf``
3. ``kepler``

The developer version of PyKE that is compatible with PyRAF is under the branch ``pyraf``.


## Acknowledgement
If you find this code useful in your research, please consider [citing](http://adsabs.harvard.edu/abs/2012ascl.soft08004S):

```
Title:	PyKE: Reduction and analysis of Kepler Simple Aperture Photometry data
Authors:          Still, Martin; Barclay, Tom
Publication:      Astrophysics Source Code Library, record ascl:1208.004
Publication Date: 08/2012
```

*This package was mostly developed by Tom Barclay ([@mrtommyb](http://www.github.com/mrtommyb)) and Martin Still.*


## Dependencies
```
numpy
scipy
astropy
matplotlib
tqdm
mdp
```

## Support
Users are welcome to open [issues](https://github.com/KeplerGO/PyKE/issues) involving any aspects of this software.

Feel free to contact us also through: keplergo@mail.arc.nasa.gov
