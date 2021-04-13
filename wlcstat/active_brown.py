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


def mscd_active(t, length_kuhn, delta, ka, gamma, b=1, num_modes=20000):
    r"""
    Compute mscd for two points on an active-Brownian polymer.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD
    length_kuhn : float
        Length of the chain (in Kuhn segments)
    delta : float
        Length of the chain between loci (in Kuhn segments)
    ka : float
        Active force rate constant
    gamma : float
        Magnitude of the active forces
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

    msd_coeff = 6 * length_kuhn * b ** 2 / (3 * np.pi ** 2)
    for p in range(1, num_modes+1, 2):
        msd_p = msd_coeff / (p ** 2) * (1 - np.exp(-p ** 2 * t / length_kuhn ** 2)
                                        + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                        1 - np.exp(-p ** 2 * t / length_kuhn ** 2)
                                        - p ** 2 / (ka * length_kuhn ** 2) * (1 - np.exp(-ka * t))))
        mscd += 8 * msd_p * np.sin(np.pi * p * delta / length_kuhn) ** 2

    return mscd


def mscd_plateau_active(length_kuhn, delta, ka, gamma, b=1, num_modes=20000):
    r"""
    Compute mscd plateau for two points on an active-Brownian polymer.

    Parameters
    ----------
    length_kuhn : float
        Length of the chain (in Kuhn segments)
    delta : float, array_like
        Length of the chain between loci (in Kuhn segments)
    ka : float
        Active force rate constant
    gamma : float
        Magnitude of the active forces
    b : float
        The Kuhn length (in desired length units).
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    mscd_plateau : (M,) np.array<float>
        result

    """
    mscd_plateau = np.zeros_like(delta)

    msd_coeff = 6 * length_kuhn * b ** 2 / (3 * np.pi ** 2)
    for p in range(1, num_modes+1, 2):
        msd_p = msd_coeff / (p ** 2) * (1 + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                        1 - p ** 2 / (ka * length_kuhn ** 2)))
        mscd_plateau += 8 * msd_p * np.sin(np.pi * p * delta / length_kuhn) ** 2

    return mscd_plateau


def msd_active(t, length_kuhn, ka, gamma, b=1, num_modes=20000):
    r"""
    Compute msd for the midpoint points on an active-Brownian polymer.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD
    length_kuhn : float
        Length of the chain (in Kuhn segments)
    ka : float
        Active force rate constant
    gamma : float
        Magnitude of the active forces
    b : float
        The Kuhn length (in desired length units).
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    msd : (M,) np.array<float>
        result

    """
    msd_coeff = 6 * b ** 2 / (3 * np.pi ** 2 * length_kuhn)
    msd_com = msd_coeff * ((1 + gamma) * t + gamma / ka * (np.exp(-ka * t) - 1))
    msd = msd_com

    msd_coeff = 6 * length_kuhn * b ** 2 / (3 * np.pi ** 2)
    for p in range(2, num_modes+1, 2):
        msd_p = msd_coeff / (p ** 2) * (1 - np.exp(-p ** 2 * t / length_kuhn ** 2)
                                        + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                        1 - np.exp(-p ** 2 * t / length_kuhn ** 2)
                                        - p ** 2 / (ka * length_kuhn ** 2) * (1 - np.exp(-ka * t))))
        msd += 2 * msd_p

    return msd


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
    gamma : float
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


def gen_conf_wlc(length_kuhn, num_beads, ka=1, gamma=0, b=1):
    return


def gen_pymol_file(r_poly, filename='r_poly.pdb', ring=False):
    r"""

    Parameters
    ----------
    r_poly : (num_beads, 3) float
        Conformation of the chain subjected to active-Brownian forces
    filename : str
        File name to write the pdb file
    ring : bool
        Boolean to close the polymer into a ring

    Returns
    -------

    """

    # Open the file
    f = open(filename, 'w')

    atomname1 = "A1"    # Chain atom type
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

