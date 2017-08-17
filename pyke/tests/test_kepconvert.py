import pytest
import csv
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits as pyfits
from astropy.time import Time as astropyTime
from pyke import kepconvert
from pyke import kepio
from ..kepio import delete

fake_lc = get_pkg_data_filename("data/golden-lc.fits")

def test_kepconvert():

    with pytest.raises(Exception,
                       message=("ERROR -- KEPCONVERT: input file myfile.fits"
                                " does not exist")):
        kepconvert("myfile.fits", "fits2excel",
                   "TIME,SAP_FLUX,SAP_FLUX_ERR,SAP_QUALITY",
                   outfile="output.txt", baddata=True,
                   overwrite=True, verbose=True)
    with pytest.raises(Exception,
                       message=("ERROR -- KEPCONVERT: conversion not supported"
                                ": fits2excel")):
        kepconvert(fake_lc, "fits2excel",
                   "TIME,SAP_FLUX,SAP_FLUX_ERR,SAP_QUALITY",
                   outfile="output.txt",
                   baddata=True, overwrite=True, verbose=True)
    # convert to csv
    kepconvert(fake_lc, "fits2csv", "TIME,SAP_FLUX", outfile="fake_lc.csv",
               baddata=True, overwrite=True, verbose=True)
    with open('fake_lc.csv', 'r') as csvfile:
        # check header
        line = csvfile.readline()
        first_line = line.split(',')
        assert first_line == ['TIME', 'SAP_FLUX\n']
    delete("fake_lc.csv", "kepconvert.log", False)
    # convert to ascii
    kepconvert(fake_lc, "fits2asc", "TIME,SAP_FLUX", outfile="fake_lc.txt",
               baddata=True, overwrite=True, verbose=True)
    with open('fake_lc.txt', 'r') as asciifile:
        lines = asciifile.readlines()
        first_line = lines[0]
        first_line = first_line.split(',')
        assert first_line == ['TIME', 'SAP_FLUX\n']
    delete("fake_lc.txt", "kepconvert.log", False)
    # time conversion
    with pytest.raises(Exception,
                       message=("ERROR -- KEPCONVERT: error converting time to"
                                " nasa: format must be one of ['jd', 'mjd',"
                                " 'decimalyear', 'unix', 'cxcsec', 'gps',"
                                " 'plot_date', 'datetime', 'iso', 'isot',"
                                " 'yday', 'fits', 'byear', 'jyear',"
                                " 'byear_str', 'jyear_str']")):
        kepconvert(fake_lc, "fits2csv", "TIME,SAP_FLUX", timeformat='nasa',
                   outfile="fake_lc.txt", baddata=True, overwrite=True,
                   verbose=True)
    kepconvert(fake_lc, "fits2csv", "TIME,SAP_FLUX", timeformat='unix',
               outfile="fake_lc.txt", baddata=True, overwrite=True,
               verbose=True)
    with open('fake_lc.txt', 'r') as asciifile:
        lines = asciifile.readlines()
        first_line = lines[0]
        first_line = first_line.split(',')
        assert first_line == ['TIME', 'SAP_FLUX\n']
        second_line = lines[1]
        second_line = second_line.split(',')
        assert second_line == (['1.461334722059185982e+09',
                                '1.500000000000000000e+01\n'])
        delete("fake_lc.txt", "kepconvert.log", False)

SUPPORTED_CONVERSION = ['jd', 'mjd', 'decimalyear', 'unix', 'cxcsec', 'gps',
                        'plot_date', 'datetime', 'iso', 'isot', 'yday', 'fits',
                        'byear', 'jyear', 'byear_str', 'jyear_str']

def test_time_conversion():

    fits = pyfits.open(fake_lc, 'readonly')
    tstart, _, _, _ = kepio.timekeys(fits, fake_lc, "kepconvert.log", False)

    for conversion in SUPPORTED_CONVERSION:
        dateobs = astropyTime(fits[1].header['DATE-OBS'])
        dateobs.format = conversion
        start = astropyTime(tstart, format='jd')
        start.format = conversion

    assert start.value == dateobs.value
