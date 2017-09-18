#!/usr/bin/env python
import os
import sys
from setuptools import setup

# Prepare and send a new release to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/pyketools*")
    sys.exit()


# Load the __version__ variable without importing the package already
exec(open('pyke/version.py').read())

entry_points = {'console_scripts': [
        'kepbls = pyke.kepbls:kepbls_main',
        'kepclip = pyke.kepclip:kepclip_main',
        'kepconvert = pyke.kepconvert:kepconvert_main',
        'kepcotrend = pyke.kepcotrend:kepcotrend_main',
        'kepdetrend = pyke.kepdetrend:kepdetrend_main',
        'kepdiffim = pyke.kepdiffim:kepdiffim_main',
        'kepdraw = pyke.kepdraw:kepdraw_main',
        'kepdynamic = pyke.kepdynamic:kepdynamic_main',
        'kepextract = pyke.kepextract:kepextract_main',
        'kepffi = pyke.kepffi:kepffi_main',
        'kepfilter = pyke.kepfilter:kepfilter_main',
        'kepflatten = pyke.kepflatten:kepflatten_main',
        'kepfold = pyke.kepfold:kepfold_main',
        'kepft = pyke.kepft:kepft_main',
        'kephead = pyke.kephead:kephead_main',
        'kepimages = pyke.kepimages:kepimages_main',
        'kepmask = pyke.kepmask:kepmask_main',
        'kepoutlier = pyke.kepoutlier:kepoutlier_main',
        'keppca = pyke.keppca:keppca_main',
        'keppixseries = pyke.keppixseries:keppixseries_main',
        'kepprf = pyke.kepprf:kepprf_main',
        'kepprfphot = pyke.kepprfphot:kepprfphot_main',
        'keprange = pyke.keprange:keprange_main',
        'kepsff = pyke.kepsff:kepsff_main',
        'kepsmooth = pyke.kepsmooth:kepsmooth_main',
        'kepstddev = pyke.kepstddev:kepstddev_main',
        'kepstitch = pyke.kepstitch:kepstitch_main',
        'keptimefix = pyke.keptimefix:keptimefix_main',
        'keptrial = pyke.keptrial:keptrial_main',
        'keptrim = pyke.keptrim:keptrim_main',
        'kepwindow = pyke.kepwindow:kepwindow_main',
]}

setup(name='pyketools',
      version=__version__,
      description="A backwards-incompatible, python3 compatible, pyraf-free "
                  "version of PyKE: a suite of tools to analyze Kepler/K2 "
                  "data",
      long_description=open('README.rst').read(),
      author='KeplerGO',
      author_email='keplergo@mail.arc.nasa.gov',
      license='MIT',
      packages=['pyke'],
      install_requires=['numpy>=1.11', 'astropy>=1.3', 'scipy>=0.17.0',
                        'matplotlib>=1.5.3', 'tqdm'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      entry_points=entry_points,
      include_package_data=True,
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
    )
