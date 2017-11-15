import pytest

HAS_MDP = True
try:
    import mdp
except ImportError:
    HAS_MDP = False


def test_import():
    from .. import keparray
    from .. import kepbls
    from .. import kepclip
    from .. import kepconvert
    from .. import kepcotrend
    from .. import kepdiffim
    from .. import kepdraw
    from .. import kepdynamic
    from .. import kepextract
    from .. import kepffi
    from .. import kepfilter
    from .. import kepfit
    from .. import kepflatten
    from .. import kepfold
    from .. import kepfourier
    from .. import kepperiodogram
    from .. import kepfunc
    from .. import kephead
    from .. import kepimages
    from .. import kepio
    from .. import kepkey
    from .. import kepmask
    from .. import kepmsg
    from .. import kepoutlier
    from .. import keppixseries
    from .. import kepplot
    from .. import kepprf
    from .. import kepprfphot
    from .. import keprange
    from .. import kepsff
    from .. import kepsmooth
    from .. import kepstat
    from .. import kepstddev
    from .. import kepstitch
    from .. import keptimefix
    from .. import keptrial
    from .. import keptrim
    from .. import kepwindow
    from .. import lightcurve
    from .. import targetpixelfile
    from .. import prf


@pytest.mark.skipif('not HAS_MDP')
def test_import_keppca():
    from .. import keppca
