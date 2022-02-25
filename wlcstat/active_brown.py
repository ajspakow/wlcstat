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


def structure_factor_active(k, t, length_kuhn, ka, gamma, b=1, num_nint = 1000, num_modes=20000):
    r"""

    Parameters
    ----------
    k
    t
    length_kuhn
    ka
    gamma
    b
    num_nint
    num_modes

    Returns
    -------

    """

    n_vec = np.linspace(0, 1, num_nint)
    dn = n_vec[1] - n_vec[0]
    structure_factor = np.zeros((np.size(k), np.size(t)))

    for i_t in range(np.size(t)):
        print(i_t)
        if np.size(t) == 1:
            t_i = t
        else:
            t_i = t[i_t]

        delta_r2 = np.zeros((num_nint, num_nint))

        # Case 1: time == 0 (not confirmed)
        if t_i == 0:
            delta_n1_n2 = np.abs(n_vec * np.ones((num_nint, 1))
                                 - np.transpose(n_vec * np.ones((num_nint, 1))))

            exp_diff_mat = (np.exp(-np.pi * np.sqrt(ka) * n_vec) * np.ones((num_nint, 1))
                            - np.transpose(np.exp(-np.pi * np.sqrt(ka) * n_vec) * np.ones((num_nint, 1))))
            sum_n1_n2 = np.abs(n_vec * np.ones((num_nint, 1))
                               + np.transpose(n_vec * np.ones((num_nint, 1))))

            delta_r2 += (1 + gamma) * np.pi ** 2 * delta_n1_n2 / (2 * length_kuhn)
            delta_r2 += - gamma * np.pi / (4 * length_kuhn * np.sqrt(ka)) * (
                                 2 * np.sinh(np.pi * np.sqrt(ka) * delta_n1_n2)
                                 - exp_diff_mat ** 2 * (np.exp(np.pi * np.sqrt(ka) * sum_n1_n2) - 1))
            delta_r2 *= length_kuhn * b ** 2 / (3 * np.pi ** 2)

        # Case 2: time != 0
        else:
            cp_coef = length_kuhn * b ** 2 / (3 * np.pi ** 2)
            for p in range(1, num_modes+1):
                cp = cp_coef / p ** 2 * (np.exp(-p ** 2 * t_i / length_kuhn ** 2)
                                         + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                                 np.exp(-p ** 2 * t_i / length_kuhn ** 2)
                                                - p ** 2 / (ka * length_kuhn ** 2) * np.exp(-ka * t_i)))
                cp0 = cp_coef / p ** 2 * (1 + gamma / (1 + p ** 2 / (ka * length_kuhn ** 2)))
                phip1_2 = np.ones((num_nint, num_nint)) * np.cos(p * np.pi * n_vec) ** 2
                phip1_phip2 = np.outer(np.cos(p * np.pi * n_vec), np.cos(p * np.pi * n_vec))
                delta_r2 += cp0 * (phip1_2 + np.transpose(phip1_2)) - 2 * cp * phip1_phip2

        for i_k in range(np.size(k)):
            integrand_n1_n2 = np.exp(- k[i_k] ** 2 * delta_r2)
            integrand_n1 = dn * (np.sum(integrand_n1_n2, axis=0)
                                 - 0.5 * integrand_n1_n2[:, 0] - 0.5 * integrand_n1_n2[:, -1])
            structure_factor[i_k, i_t] = dn * (np.sum(integrand_n1)
                                               - 0.5 * integrand_n1[0] - 0.5 * integrand_n1[-1])

    return structure_factor


def flow_spec_active(k, t, length_kuhn, ka, gamma, b=1, num_nint = 1000, num_modes=20000):
    r"""

    Parameters
    ----------
    k
    t
    length_kuhn
    ka
    gamma
    b
    num_nint
    num_modes

    Returns
    -------

    """

    n_vec = np.linspace(0, 1, num_nint)
    dn = n_vec[1] - n_vec[0]
    structure_factor = np.zeros((np.size(k), np.size(t)))

    for i_t in range(np.size(t)):
        print(i_t)
        if np.size(t) == 1:
            t_i = t
        else:
            t_i = t[i_t]
        delta_r2 = np.zeros((num_nint, num_nint))
        vel_int2 = np.zeros((num_nint, num_nint))
        cp_coef = length_kuhn * b ** 2 / (3 * np.pi ** 2)
        vel_int1 = np.ones((num_nint, num_nint)) * (2 * b ** 2 / (np.pi ** 2 * length_kuhn)) * (
            (1 + gamma) * t_i + gamma / ka * (np.exp(-ka * t_i) - 1))

        for p in range(1, num_modes+1):
            cp = cp_coef / p ** 2 * (np.exp(-p ** 2 * t_i / length_kuhn ** 2)
                                     + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                             np.exp(-p ** 2 * t_i / length_kuhn ** 2)
                                             - p ** 2 / (ka * length_kuhn ** 2) * np.exp(-ka * t_i)))
            cp0 = cp_coef / p ** 2 * (1 + gamma / (1 + p ** 2 / (ka * length_kuhn ** 2)))
            phip1_2 = np.ones((num_nint, num_nint)) * np.cos(p * np.pi * n_vec) ** 2
            phip1_phip2 = np.outer(np.cos(p * np.pi * n_vec), np.cos(p * np.pi * n_vec))
            delta_r2 += cp0 * (phip1_2 - np.transpose(phip1_2)) ** 2
            vel_int1 += 4 * (cp0 - cp) * phip1_phip2
            vel_int2 += 2 * (cp0 - cp) * (phip1_2 - phip1_phip2)

        vel_int2 = vel_int2 * np.transpose(vel_int2)

        for i_k in range(np.size(k)):
            integrand_n1_n2 = (vel_int1 + vel_int2 * k[i_k] ** 2) * np.exp(- k[i_k] ** 2 * delta_r2)
            integrand_n1 = dn * (np.sum(integrand_n1_n2, axis = 0)
                                 - 0.5 * integrand_n1_n2[:, 0] - 0.5 * integrand_n1_n2[:, -1])
            structure_factor[i_k, i_t] = length_kuhn / (t_i ** 2) * dn * (np.sum(integrand_n1)
                                               - 0.5 * integrand_n1[0] - 0.5 * integrand_n1[-1])

    return structure_factor


def phi_active(t, length_kuhn, ka, gamma, b=1, num_modes=20000):
    r"""

    Parameters
    ----------
    t
    length_kuhn
    ka
    gamma
    b
    num_modes

    Returns
    -------
    phia

    """
    phia = np.zeros_like(t)

    phia_coeff = 8 / (np.pi ** 2) / (
                1 + gamma * (1 - 2.0 / (np.pi * length_kuhn * np.sqrt(ka)) * np.tanh(np.pi * length_kuhn * np.sqrt(ka) / 2)))

    for p in range(1, num_modes + 1, 2):
        phi_p = 1 / (p ** 2) * (np.exp(-p ** 2 * t / length_kuhn ** 2)
                                + gamma / (1 - p ** 4 / (ka ** 2 * length_kuhn ** 4)) * (
                                        np.exp(-p ** 2 * t / length_kuhn ** 2)
                                        - p ** 2 / (ka * length_kuhn ** 2) * np.exp(-ka * t)))
        phia += phi_p

    phia = phia_coeff * phia

    return phia


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





def gen_conf_rouse_active(length_kuhn, num_beads, ka=1, gamma=0, b=1, force_calc=False, num_modes=10000):
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

    if not force_calc:
        # Calculate the conformation without determining the active and Brownian force
        r_poly = np.zeros((num_beads, 3))

        ind = np.arange(num_beads)
        for p in range(1, num_modes + 1):
            sig_fp_tilde = np.sqrt(length_kuhn ** 2 * gamma * ka / p ** 2)
            ka_tilde = ka * length_kuhn ** 2 / p ** 2
            sig_xp = np.sqrt(1 + sig_fp_tilde ** 2 / (1 + ka_tilde))
            xp_tilde = np.random.randn(3) * sig_xp
            phi = np.sqrt(2) * np.cos(p * np.pi * ind / (num_beads - 1))
            r_poly += np.outer(phi, xp_tilde) * (np.sqrt(length_kuhn) / p) * b / np.sqrt(3 * np.pi ** 2)

        return r_poly
    else:
        # Calculate the conformation, Brownian force, and active force
        r_poly = np.zeros((num_beads, 3))
        f_active = np.zeros((num_beads, 3))

        ind = np.arange(num_beads)
        for p in range(1, num_modes + 1):
            sig_fp_tilde = np.sqrt(length_kuhn ** 2 * gamma * ka / p ** 2)
            ka_tilde = ka * length_kuhn ** 2 / p ** 2
            fp_tilde = np.random.randn(3) * sig_fp_tilde
            mu_xp = fp_tilde / (1 + ka_tilde)
            sig_xp = np.sqrt(1 + sig_fp_tilde ** 2 * ka_tilde / (1 + ka_tilde) ** 2)
            xp_tilde = np.random.randn(3) * sig_xp + mu_xp
            phi = np.sqrt(2) * np.cos(p * np.pi * ind / (num_beads - 1))
            r_poly += np.outer(phi, xp_tilde) * (np.sqrt(length_kuhn) / p) * b / np.sqrt(3 * np.pi ** 2)
            f_active += np.outer(phi, fp_tilde) * (p / np.sqrt(length_kuhn)) * np.sqrt(3 * np.pi ** 2) / b

        return r_poly, f_active


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
    atomname2 = "A2"    # Chain atom type
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

