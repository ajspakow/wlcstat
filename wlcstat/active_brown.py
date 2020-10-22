r"""
Active-Brownian dynamics

Notes
-----
Detailed derivations for MSD and MSCD are found in "Interplay of active and thermal fluctuations in polymer dynamics"

"""
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
# import mpmath

from functools import lru_cache
from pathlib import Path
import os


def mscd_active(t, D, ka, fa, Ndel, N, b=1, num_modes=20000):
    r"""
    Compute mscd for two points on an active-Brownian polymer.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    ka : float
        Dimensionless active-force rate constant
    fa : float
        Dimensionless active-force magnitude
    Ndel : float
        Distance from the last linkage site to the measured site. This ends up
        being (1/2)*separation between the loci (in Kuhn lengths).
    N : float
        The full length of the linear polymer (in Kuhn lengths).
    b : float
        The Kuhn length (in desired length units).
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    mscd : (M,) np.array<float>
        result
    """

    mscd = np.zeros_like(t)

    k1 = 3 * np.pi ** 2 / (N * (b ** 2))
    sum_coeff = 48 / k1
    exp_coeff = k1 * D / N
    sin_coeff = np.pi * Ndel / N

    for p in range(1, num_modes+1, 2):
        active_coeff = 1 + 0.5 * fa ** 2 * ka / (ka + p ** 2)
        mscd += active_coeff * (1 / p ** 2) * (1 - np.exp(-exp_coeff * (p ** 2) * t)) * np.sin(sin_coeff * p) ** 2

    return sum_coeff * mscd


def msd_active(t, D, ka, fa, N, b=1, num_modes=20000):
    r"""
    Compute msd for two points on an active-Brownian polymer.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    ka : float
        Dimensionless active-force rate constant
    fa : float
        Dimensionless active-force magnitude
    N : float
        The full length of the linear polymer (in Kuhn lengths).
    b : float
        The Kuhn length (in desired length units).
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    msd : (M,) np.array<float>
        result

    """
    t_rouse = (N * b) ** 2 / (D * 3 * np.pi ** 2)
    t_dim = t / t_rouse
    msd = np.zeros_like(t)

    k1 = 3 * np.pi ** 2 / (N * (b ** 2))
    sum_coeff = 12 / k1
    exp_coeff = k1 * D / N

    msd_com = D * t_rouse / N * (6 * t_dim + 3 * fa ** 2 * (t_dim - 1 / ka + 1 / ka * np.exp(- ka * t_dim)))

    for p in range(2, num_modes+1, 2):
        active_coeff = 1 + 0.5 * fa ** 2 * ka / (ka + p ** 2)
        msd += active_coeff * (1 / p ** 2) * (1 - np.exp(-exp_coeff * (p ** 2) * t))

    return sum_coeff * msd + msd_com
