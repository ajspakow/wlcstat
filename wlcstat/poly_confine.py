"""
poly_confine

Module containing functions for evaluating the
statistical behavior of confined polymers

"""
import numpy as np
from numba import jit
import matplotlib.pyplot as plt


def eval_r2_free(delta, length_kuhn, b, a=1, n_max=1000):
    r"""
    Evaluate the average :math:`\langle R^{2} \rangle` for a polymer confined within
    a sphere with free ends

    Parameters
    ----------
    delta : float
        The position of the segment on the polymer for the average
    length_kuhn : float
        The total length of the chain (in Kuhn segments)
    b : float
        The value of the Kuhn length (dimensions of length
    a : float
        The radius of the confining sphere
    n_max : int
        The total number of :math:`n` modes in the Bessel function expansion

    Returns
    -------
    r2 : float
        The average :math:`\langle R^{2} \rangle` for a polymer

    """
    r2 = np.zeros_like(delta)
    norm = 0

    for n in range(1, n_max):
        cn = np.exp(-(n * np.pi) ** 2 * (b / a) ** 2 * (length_kuhn - delta) / 6)
        norm += np.exp(-(n * np.pi) ** 2 * (b / a) ** 2 * length_kuhn / 6) / (n * np.pi) ** 2
        for n0 in range(1, n_max):
            cn0 = np.exp(-(n0 * np.pi) ** 2 * (b / a) ** 2 * delta / 6)
            if n == n0:
                inn0 = 1 / 6 - 1 / (4 * np.pi ** 2 * n ** 2)
            else:
                inn0 = 4 * n * n0 * (-1) ** n * (-1) ** n0 / (n ** 2 - n0 ** 2) ** 2 / np.pi ** 2
            r2 += (-1) ** n * (-1) ** n0 * inn0 * cn * cn0 / (n * np.pi) / (n0 * np.pi)

    norm *= 8 * np.pi
    r2 *= 16 * np.pi * a ** 2 / norm

    return r2


def eval_r2_surf(delta, length_kuhn, b, a=1, n_max=1000):
    r"""
    Evaluate the average :math:`\langle R^{2} \rangle` for a polymer confined within
    a sphere with surface-attached ends

    Parameters
    ----------
    delta : float
        The position of the segment on the polymer for the average
    length_kuhn : float
        The total length of the chain (in Kuhn segments)
    b : float
        The value of the Kuhn length (dimensions of length
    a : float
        The radius of the confining sphere
    n_max : int
        The total number of :math:`n` modes in the Bessel function expansion

    Returns
    -------
    r2 : float
        The average :math:`\langle R^{2} \rangle` for a polymer

    """
    r2 = np.zeros_like(delta)
    norm = 0

    for n in range(1, n_max):
        cn = np.exp(-(n * np.pi) ** 2 * (b / a) ** 2 * (length_kuhn - delta) / 6)
        norm += np.exp(-(n * np.pi) ** 2 * (b / a) ** 2 * length_kuhn / 6) * (n * np.pi) ** 2
        for n0 in range(1, n_max):
            cn0 = np.exp(-(n0 * np.pi) ** 2 * (b / a) ** 2 * delta / 6)
            if n == n0:
                inn0 = 1 / 6 - 1 / (4 * np.pi ** 2 * n ** 2)
            else:
                inn0 = 4 * n * n0 * (-1) ** n * (-1) ** n0 / (n ** 2 - n0 ** 2) ** 2 / np.pi ** 2
            r2 += (-1) ** n * (-1) ** n0 * inn0 * cn * cn0 * (n * np.pi) * (n0 * np.pi)

    norm *= 8 * np.pi
    r2 *= 16 * np.pi * (a ** 2) / norm

    return r2
