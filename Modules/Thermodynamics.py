import numpy as np
from scipy.special import tanh, sinh, cosh

def w_to_a(w, T):
    """
    Convert the frequency w to the average displacement a.
    BEWARE: The frequency w should be in Hartree.

    This is a inner utility function that replace the original w_to_a
    in the thermodynamic fortran module, do not use it directly.

    Parameters
    ----------
    w : array_like
        Frequency in Hartree.

    """
    n = len(w)
    a = np.zeros(n)
    if T == 0.0:
        a[:] = np.sqrt(1.0 / (2.0 * w))
    else:
        a[:] = np.sqrt((1.0 / tanh(0.5 * w * 315774.65221921849 / T)) / (2.0 * w))
    return a
