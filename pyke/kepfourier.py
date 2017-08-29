from .utils import PyKEArgumentHelpFormatter
import numpy as np

def ft(x, y, f1, f2, df, verbose):
    """
    Compute the Fourier transform of a signal ``y`` with support ``x``
    in the bandwidth between the frequencies ``f1`` to ``f2``. ``df``
    is the frequency of resolution.

    Parameters
    ----------
    x : ndarray
        Support of the signal (usually time)
    y : ndarray
        Signal values
    f1, f2 : float
        Initial and last frequencies
    df : float
        Frequency of resolution

    Returns
    -------
    fr : ndarray
        Frequency support
    power : ndarray
        Power spectrum
    """
    notnans = (~np.isnan(x)) & (~np.isnan(y))
    x = x[notnans]
    y = y[notnans]

    ft_real, ft_imag, power, fr = [], [], [], []
    nstep = 0
    len_x = len(x)
    for freq in np.arange(f1, f2, df):
        ft_real.append(0.0)
        ft_imag.append(0.0)
        omega = 2.0 * np.pi * freq
        for i in range(len_x):
            expo = omega * x[i]
            c = np.cos(expo)
            s = np.sin(expo)
            ft_real[-1] += y[i] * c
            ft_imag[-1] += y[i] * s
        fr.append(freq)
        power.append((ft_real[-1]**2 + ft_imag[-1]**2) / len_x ** 2)
        nstep += 1
        if verbose:
            print('Step: {0}  Period: {1} (d)  Power: {2}'
                  .format(nstep, 1.0 / fr[-1], power[-1]))

    fr = np.array(fr, dtype='float32')
    power = np.array(power, dtype='float32')

    return fr, power
