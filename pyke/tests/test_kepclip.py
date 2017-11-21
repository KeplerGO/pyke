from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import numpy as np
import os
from ..kepio import delete
from ..kepclip import kepclip


TPF_filename  =  get_pkg_data_filename("data/testtpf.fits")
LC_filename  =  get_pkg_data_filename("data/dip_llc.fits")
OUT_filename  =  "out.fits"


def test_kepclip_lc():
    """
    Test the light curve clipping.
    Uses cycle 8 of Tabby's star.
    Clips out the large dip and saves it to out.fits
    """

    h = fits.open(LC_filename)
    t = h[1].data['TIME'][200:400]
    ranges = '{:12.12f}'.format(t[0] + 2454833.) + ', {:12.12f}'.format(t[-1] + 2454833.)
    kepclip(infile=LC_filename, outfile=OUT_filename, ranges=ranges, overwrite=True)
    j = fits.open(OUT_filename)

    #The new file should be smaller
    assert os.path.getsize(LC_filename) > os.path.getsize(OUT_filename)
    #The new file should also contain the dip
    assert np.min(h[1].data['SAP_FLUX']) == np.min(j[1].data['SAP_FLUX'])
    #The median of the new file should be lower, as we are clipping out the dip.
    assert np.median(h[1].data['SAP_FLUX']) > np.median(j[1].data['SAP_FLUX'])
    #New file should have 200 entries
    assert len(j[1].data) == 200
    #Fits files should be the same shape
    assert len(h) == len(j)
    j.close()
    h.close()
    delete(OUT_filename, "log_kepclip.txt", False)

def test_kepclip_tpf():
    """
    Test the TPF clipping
    Uses a clipped TPF of Tabby's star.
    """
    h = fits.open(TPF_filename)
    t = h[1].data['TIME'][0:2]
    ranges = '{:12.12f}'.format(t[0] + 2454833.) + ', {:12.12f}'.format(t[-1] + 2454833.)

    kepclip(infile=TPF_filename, outfile=OUT_filename, ranges=ranges, overwrite=True)
    j = fits.open(OUT_filename)
    #The new file should be smaller
    assert os.path.getsize(TPF_filename) > os.path.getsize(OUT_filename)

    #Fits files should be the same shape
    assert len(h) == len(j)

    #The first frame in both files should be identical
    assert np.allclose(h[1].data['FLUX'][0], j[1].data['FLUX'][0])

    #Test if we clip out the middle file it also works
    t = h[1].data['TIME'][4:15]
    ranges = '{:12.12f}'.format(t[0] + 2454833.) + ', {:12.12f}'.format(t[-1] + 2454833.)

    #close the old file and clip it again
    j.close()
    kepclip(infile=TPF_filename, outfile=OUT_filename, ranges=ranges, overwrite=True)
    j = fits.open(OUT_filename)

    #The middle frame should be identical
    assert np.allclose(h[1].data['FLUX'][4], j[1].data['FLUX'][0])
    j.close()
    h.close()
    delete(OUT_filename, "log_kepclip.txt", False)
