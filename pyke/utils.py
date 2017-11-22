import numpy as np
import matplotlib.pyplot as plt
from argparse import HelpFormatter, SUPPRESS, OPTIONAL, ZERO_OR_MORE
from astropy.visualization import (PercentileInterval, ImageNormalize,
                                   SqrtStretch, LogStretch, LinearStretch)


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


class KeplerQualityFlags(object):
    """
    This class encodes the meaning of the various Kepler QUALITY bitmask flags,
    as documented in the Kepler Archive Manual (Table 2.3).
    """
    AttitudeTweak = 1
    SafeMode = 2
    CoarsePoint = 4
    EarthPoint = 8
    ZeroCrossing = 16
    Desat = 32
    Argabrightening = 64
    ApertureCosmic = 128
    ManualExclude = 256
    SensitivityDropout = 1024
    ImpulsiveOutlier = 2048
    ArgabrighteningOnCCD = 4096
    CollateralCosmic = 8192
    DetectorAnomaly = 16384
    NoFinePoint = 32768
    NoData = 65536
    RollingBandInAperture = 131072
    RollingBandInMask = 262144
    PossibleThrusterFiring = 524288
    ThrusterFiring = 1048576

    # Which is the recommended QUALITY mask to identify bad data?
    DEFAULT_BITMASK = (AttitudeTweak | SafeMode | CoarsePoint | EarthPoint |
                       Desat | ApertureCosmic | ManualExclude | NoData | ThrusterFiring)

    # This bitmask includes flags that are known to identify both good and bad cadences.
    # Use it wisely.
    CONSERVATIVE_BITMASK = (DEFAULT_BITMASK | SensitivityDropout | CollateralCosmic
                            | PossibleThrusterFiring)

    # Using this bitmask only QUALITY == 0 cadences will remain
    QUALITY_ZERO_BITMASK = 2096639

    # Pretty string descriptions for each flag
    STRINGS = {
        1: "Attitude tweak",
        2: "Safe mode",
        4: "Coarse point",
        8: "Earth point",
        16: "Zero crossing",
        32: "Desaturation event",
        64: "Argabrightening",
        128: "Cosmic ray in optimal aperture",
        256: "Manual exclude",
        1024: "Sudden sensitivity dropout",
        2048: "Impulsive outlier",
        4096: "Argabrightening on CCD",
        8192: "Cosmic ray in collateral data",
        16384: "Detector anomaly",
        32768: "No fine point",
        65536: "No data",
        131072: "Rolling band in optimal aperture",
        262144: "Rolling band in full mask",
        524288: "Possible thruster firing",
        1048576: "Thruster firing"
    }

    @classmethod
    def decode(cls, quality):
        """Converts a Kepler QUALITY value into a list of human-readable strings.

        This function takes the QUALITY bitstring that can be found for each
        cadence in Kepler/K2's pixel and light curve files and converts into
        a list of human-readable strings explaining the flags raised (if any).

        Parameters
        ----------
        quality : int
            Value from the 'QUALITY' column of a Kepler/K2 pixel or lightcurve file.

        Returns
        -------
        flags : list of str
            List of human-readable strings giving a short description of the
            quality flags raised.  Returns an empty list if no flags raised.
        """
        result = []
        for flag in cls.STRINGS.keys():
            if quality & flag > 0:
                result.append(cls.STRINGS[flag])
        return result

def plot_image(image, scale='linear', origin='lower', xlabel='Pixel Column Number',
               ylabel='Pixel Row Number', title=None, **kwargs):
        vmin, vmax = PercentileInterval(95.).get_limits(image)
        if scale == 'linear':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
        elif scale == 'sqrt':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
        elif scale == 'log':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
        else:
            raise ValueError("scale {} is not available.".format(scale))

        plt.imshow(image, origin=origin, norm=norm, **kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.colorbar(norm=norm)
