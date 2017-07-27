import kepclip as kc
from pyke.kepio import delete
from astropy.io import fits
import numpy as np

def test_kepclip_tpf():
	kc.kepclip('testtpf.fits',
		'out.fits',
		'2455568.36367266,2455635.3033717',
		overwrite=True)
	h=fits.open('dip.fits')
	j=fits.open('out.fits')

	assert np.allclose(h[1].data['FLUX'][0],j[1].data['FLUX'][0])




def test_kepclip_tpf_middle():
	h=fits.open('testtpf.fits')
	t=h[1].data['TIME'][4:15]
	ranges='{0:.15f}'.format(t[0]+2454833)+',{0:.15f}'.format(t[-1]+2454833)
	
	kc.kepclip('dip.fits',
		'out.fits',
		ranges,
		overwrite=True)
	j=fits.open('out.fits')

	assert np.allclose(h[1].data['FLUX'][4],j[1].data['FLUX'][0])



def test_kepclip_tpf_timespecifyproblem():

	'''
		This misses out the first/last file if you pass it a list of exact times.
	'''

	h=fits.open('testtpf.fits')
	t=h[1].data['TIME'][4:17]
	ranges='{0:.15f}'.format(t[0]+2454833)+',{0:.15f}'.format(t[-1]+2454833)
	print t[0]+2454833,t[-1]+2454833
	kc.kepclip('dip.fits',
		'out.fits',
		ranges,
		overwrite=True)
	j=fits.open('out.fits')

	assert len(t) - len(j[1].data['TIME']) == 0
#	assert np.allclose(h[1].data['FLUX'][150],j[1].data['FLUX'][0])
