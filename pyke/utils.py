import numpy as np
from argparse import HelpFormatter, SUPPRESS, OPTIONAL, ZERO_OR_MORE


class PyKEArgumentHelpFormatter(HelpFormatter):
    """Help message formatter which adds default values to argument help,
    except for boolean arguments.
    """

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if (action.default is not SUPPRESS
                and not isinstance(action.default, bool)):
                defaulting_nargs = [OPTIONAL, ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help


def module_output_to_channel(module, output):
    """Returns the cCCD channel number for a given module and output pair."""
    if module < 1 or module > 26:
        raise ValueError("Module number must be in range 1-26.")
    if output < 1 or output > 4:
        raise ValueError("Output number must be 1, 2, 3, or 4.")
    return _get_channel_lookup_array()[module, output]


def channel_to_module_output(channel):
    """Returns a (module, output) pair given a CCD channel number."""
    if channel < 1 or channel > 88:
        raise ValueError("Channel number must be in the range 1-88.")
    lookup = _get_channel_lookup_array()
    lookup[:,0] = 0
    modout = np.where(lookup == channel)
    return (modout[0][0], modout[1][0])


def _get_channel_lookup_array():
    """Returns a lookup table which maps (module, output) onto channel."""
    # In the array below, channel == array[module][output]
    # Note: modules 1, 5, 21, 25 are the FGS guide star CCDs.
    return np.array([
       [0,     0,    0,    0,    0],
       [1,    85,    0,    0,    0],
       [2,     1,    2,    3,    4],
       [3,     5,    6,    7,    8],
       [4,     9,   10,   11,   12],
       [5,    86,    0,    0,    0],
       [6,    13,   14,   15,   16],
       [7,    17,   18,   19,   20],
       [8,    21,   22,   23,   24],
       [9,    25,   26,   27,   28],
       [10,   29,   30,   31,   32],
       [11,   33,   34,   35,   36],
       [12,   37,   38,   39,   40],
       [13,   41,   42,   43,   44],
       [14,   45,   46,   47,   48],
       [15,   49,   50,   51,   52],
       [16,   53,   54,   55,   56],
       [17,   57,   58,   59,   60],
       [18,   61,   62,   63,   64],
       [19,   65,   66,   67,   68],
       [20,   69,   70,   71,   72],
       [21,   87,    0,    0,    0],
       [22,   73,   74,   75,   76],
       [23,   77,   78,   79,   80],
       [24,   81,   82,   83,   84],
       [25,   88,    0,    0,    0],
       ])
