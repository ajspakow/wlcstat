"""
wlc_lcpoly

Module containing functions for evaluating the
self-consistent field theory for polymer solutions, including nematic
liquid crystallinity

"""
import numpy as np
from numba import jit
import matplotlib.pyplot as plt


def m_lcpoly(length_kuhn, lam, alpha_max=50):
    r"""
    Calculate the single-polymer partition function

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    lam : float
        The value of the quadrupole field :math:`\lambda`
    alpha_max : int
        Maximum number of poles evaluated (default 50)

    Returns
    -------
    m_val : float
        Evaluated order parameter

    """

    if lam < 50:
        d_lam = 1e-8

        q_val_p = q_lcpoly(length_kuhn, lam + d_lam, alpha_max)
        q_val_m = q_lcpoly(length_kuhn, lam - d_lam, alpha_max)
        m_val = 1.5 * (np.log(q_val_p) - np.log(q_val_m)) / (2 * d_lam * length_kuhn)

    else:
        m_val = (1 - 3 / 2 / np.sqrt(lam) - 3 / 16 / lam ** (3 / 2) - 3 / 8 / lam ** 2
                - 207 / 256 / lam ** (5 / 2) - 123 / 64 / lam ** 3 - 10215 / 2048 / lam ** (7 / 2))

    return m_val


def q_lcpoly(length_kuhn, lam, alpha_max = 50):
    r"""
    Calculate the single-polymer partition function

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    lam : float
        The value of the quadrupole field :math:`\lambda`
    alpha_max : int
        Maximum number of poles evaluated (default 50)

    Returns
    -------
    q_val : float
        Evaluated single-polymer partition function

    """

    # Evaluate the partition function by determining the poles of the h_matrix
    if lam == 0:
        q_val = 1
    elif lam < 50 and lam != 0:
        poles = eval_poles_lcpoly(lam, alpha_max)
        q_val = 0

        for i_pole in range(len(poles)):
            l = 2 * alpha_max
            a_l = l / np.sqrt(4 * l ** 2 - 1)
            a_lp1 = (l + 1) / np.sqrt(4 * (l + 1) ** 2 - 1)
            beta_l = a_lp1 ** 2 + a_l ** 2 - 1 / 3
            dj = 1
            j = -poles[i_pole] + l * (l + 1) - beta_l * lam
            for l in range(2 * (alpha_max - 1), -1, -2):
                if l > 0:
                    a_l = l / np.sqrt(4 * l ** 2 - 1)
                else:
                    a_l = 0
                a_lp1 = (l + 1) / np.sqrt(4 * (l + 1) ** 2 - 1)
                a_lp2 = (l + 2) / np.sqrt(4 * (l + 2) ** 2 - 1)
                alpha_lp2 = a_lp2 * a_lp1
                beta_l = a_lp1 ** 2 + a_l ** 2 - 1 / 3
                dj = 1 + (alpha_lp2 * lam) ** 2 * dj / j ** 2
                j = -poles[i_pole] + l * (l + 1) - beta_l * lam - (alpha_lp2 * lam) ** 2 / j

            q_val = q_val + np.exp(- poles[i_pole] * length_kuhn) / dj

    # Evaluate the partition function using the large-lam expansion
    else:
        q_val = np.exp(length_kuhn * (2 / 3 * lam - 2 * np.sqrt(lam) + 1 + 1 / 4 / np.sqrt(lam)
                                      + 1 / 4 / lam + 23 / 64 / lam ** (3 / 2)
                                      + 41 / 64 / lam ** 2 + 681 / 512 / lam ** (5 / 2)))
    return q_val


def eval_poles_lcpoly(lam, alpha_max=50):
    r"""
    Evaluate the poles for a wormlike nematic liquid crystalline polymer

    Parameters
    ----------
    lam : float
        The value of the quadrupole field :math:`\lambda`
    alpha_max : int
        Maximum number of poles evaluated (default 50)

    Returns
    -------
    poles : float
        Evaluated poles for the given :math:`\lambda`

    """

    l_vec = 2 * np.arange(alpha_max)

    a_l = np.append(0, l_vec[1:alpha_max] / np.sqrt(4 * l_vec[1:alpha_max] ** 2 - 1))
    a_lp1 = (l_vec + 1) / np.sqrt(4 * (l_vec + 1) ** 2 - 1)
    a_lp2 = (l_vec + 2) / np.sqrt(4 * (l_vec + 2) ** 2 - 1)
    alpha_lp2 = a_lp2 * a_lp1
    beta_l = a_lp1 ** 2 + a_l ** 2 - 1 / 3
    h_matrix = (np.diag(l_vec * (l_vec + 1) - beta_l * lam, 0)
                - lam * np.diag(alpha_lp2[0:(alpha_max-1)], 1)
                - lam * np.diag(alpha_lp2[0:(alpha_max-1)], -1))

    # Find the poles as the eigenvalues of the h-matrix
    poles = np.linalg.eigvals(h_matrix)
    poles = np.sort(poles)[::-1]

    return poles
