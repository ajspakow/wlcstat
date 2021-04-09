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


def gen_conf_rouse_active(length_kuhn, num_beads, ka=1, gamma=0, b=1, num_modes=10000):
    r"""
    Generate a discrete chain based on the active-Brownian Rouse model

    Parameters
    ----------
    length_kuhn : float
        Length of the chain (in Kuhn segments)
    num_beads : int
        Number of beads in the discrete chain
    ka : float
        Active force rate constant
    fa : float
        Magnitude of the active forces
    b : float
        Kuhn length
    num_modes : int
        Number of Rouse modes in calculation

    Returns
    -------
    r_poly : (num_beads, 3) float
        Conformation of the chain subjected to active-Brownian forces

    """

    r_poly = np.zeros((num_beads, 3))
    k1 = 3 * np.pi ** 2 / (length_kuhn * (b ** 2))
    ind = np.arange(num_beads)

    for p in range(1, num_modes + 1):
        kp = k1 * p ** 2
        xp_mag = np.sqrt(1 / kp * (1 + gamma / (1 + p ** 2 / (ka * length_kuhn ** 2))))
        xp = np.random.randn(3) * xp_mag
        phi = np.sqrt(2 / length_kuhn) * np.cos(p * np.pi * ind / (num_beads - 1))
        r_poly += np.outer(phi, xp)

    return r_poly


def gen_pymol_file(r_poly, filename='r_poly.pdb', ring=False):
    r"""

    Parameters
    ----------
    r_poly
    filename
    ring

    Returns
    -------

    """

    # Open the file
    f = open(filename, 'w')

    atomname1 = "A1"    # Chain atom type
    atomname2 = "A2"    # Ribbon atom type
    atomname3 = "A3"    # Extra atom type
    atomname4 = "A4"    # Extra atom type
    atomname5 = "A5"    # Extra atom type
    atomname6 = "A6"    # Extra atom type
    atomname7 = "A7"    # Extra atom type
    atomname8 = "A8"    # Extra atom type
    resname = "SSN"     # Type of residue (UNKnown/Single Stranded Nucleotide)
    chain = "A"         # Chain identifier
    resnum = 1
    numresidues = len(r_poly[:, 0])
    descrip = "Pseudo atom representation of DNA"
    chemicalname = "Body and ribbon spatial coordinates"

    # Write the preamble to the pymol file

    f.write('HET    %3s  %1s%4d   %5d     %-38s\n' % (resname, chain, resnum, numresidues, descrip))
    f.write('HETNAM     %3s %-50s\n' % (resname, chemicalname))
    f.write('FORMUL  1   %3s    C20 N20 P21\n' % (resname))

    # Write the conformation to the pymol file

    for ind in range(numresidues):
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                (ind + 1, atomname1, resname, chain, r_poly[ind, 0], r_poly[ind, 1], r_poly[ind, 2], 1.00, 1.00))

    # Define the connectivity in the chain

    if ring:
        f.write('CONECT%5d%5d%5d\n' % (1, 2, numresidues))
    else:
        f.write('CONECT%5d%5d\n' % (1, 2))

    for ind in range(2, numresidues):
        f.write('CONECT%5d%5d%5d\n' % (ind, ind - 1, ind + 1))

    if ring:
        f.write('CONECT%5d%5d%5d\n' % (numresidues, numresidues - 1, 1))
    else:
        f.write('CONECT%5d%5d\n' % (numresidues, numresidues - 1))

    # Close the file
    f.write('END')

    f.close()

    return

