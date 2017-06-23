# Inserts the directory containing links
# to all Python-based tasks in the STSDAS tree to
# the default Python path.

import iraf,sys

# define path to top level Python directory in STSDAS

_path = iraf.osfn('kepler$')

# if directory not already in PYTHONPATH,...

if _path not in sys.path:
    # Add the directory to path.
    sys.path.insert(1,_path)
