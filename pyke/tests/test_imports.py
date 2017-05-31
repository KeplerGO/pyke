HAS_MDP = True
try:
    import mdp
except:
    HAS_MDP = False

def test_import():
    from .. import kepio
    from .. import kepmask
    from .. import kepmsg
    from .. import kepkey
    from .. import kepplot
    from .. import kepstat
    from .. import kepfunc
    from .. import keparray
    from .. import kepprf
    from .. import kepfit
    from .. import kepfold
    from .. import kepclip
    from .. import kepsmooth
    from .. import kepflatten
    from .. import kepdynamic
    from .. import kepstddev
    from .. import kepwindow
    from .. import keptrim

@pytest.mark.skipif(not HAS_MDP)
def test_import_keppca():
    from .. import keppca
