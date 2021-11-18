
from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np
import scipy as sp


# Quadratic order density functions

def s2_wlc_diblock(k_val_vector, length_kuhn, fa, dimensions=3, alpha_max=25):
    r"""
    s2_wlc_randcopoly - Evaluate the 2-point structure factor for the wormlike chain model
    for a diblock copolymer

    Parameters
    ----------
    k_val_vector : float (array)
        The value of the Fourier vector magnitude :math:`K`
    length_kuhn : float (array)
        The length of the polymer chains in Kuhn lengths
    fa : float
        Fraction of A monomers in the random copolymer
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    s2 : float (vector)
        Structure factor for the wormlike chain model for every k_val in k_val_vector

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """
    if type(length_kuhn) == float or type(length_kuhn) == int:
        length_kuhn = np.array([length_kuhn])

    s2_aa = np.zeros((len(k_val_vector), len(length_kuhn)), dtype=type(1+1j))
    s2_ab = np.zeros((len(k_val_vector), len(length_kuhn)), dtype=type(1+1j))
    s2_bb = np.zeros((len(k_val_vector), len(length_kuhn)), dtype=type(1+1j))

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_kuhn, dimensions=3)
    k_cutoff_factor = 0.01
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    for ind_k_val in range(0, len(k_val_vector)):
        k_val = k_val_vector[ind_k_val]
        lam_max = alpha_max

        poles = eval_poles(k_val, 0, dimensions, alpha_max)
        residues = eval_residues(k_val, 0, poles, True, dimensions, lam_max, alpha_max, k_val_cutoff=10 ** -3)
        residue_zero, ddp_residue_zero = eval_residue_zero(k_val, 0, dimensions)

        for ind_length in range(0, len(length_kuhn)):
            s2_aa[ind_k_val, ind_length] = (fa * length_kuhn[ind_length] * residue_zero + ddp_residue_zero)
            s2_bb[ind_k_val, ind_length] = ((1 - fa) * length_kuhn[ind_length] * residue_zero + ddp_residue_zero)
            s2_ab[ind_k_val, ind_length] = -ddp_residue_zero

            for alpha in range(0, alpha_max):
                s2_aa[ind_k_val, ind_length] += (np.exp(poles[alpha] * fa * length_kuhn[ind_length]) *
                                                  residues[alpha] / poles[alpha] ** 2)
                s2_bb[ind_k_val, ind_length] += (np.exp(poles[alpha] * (1 - fa) * length_kuhn[ind_length]) *
                                                  residues[alpha] / poles[alpha] ** 2)
                s2_ab[ind_k_val, ind_length] += ((np.exp(poles[alpha] * length_kuhn[ind_length])
                                                 - np.exp(poles[alpha] * fa * length_kuhn[ind_length])
                                                 - np.exp(poles[alpha] * (1 - fa) * length_kuhn[ind_length])) *
                                                 residues[alpha] / poles[alpha] ** 2)

            s2_aa[ind_k_val, ind_length] *= 2 / length_kuhn[ind_length] ** 2
            s2_bb[ind_k_val, ind_length] *= 2 / length_kuhn[ind_length] ** 2
            s2_ab[ind_k_val, ind_length] *= 1 / length_kuhn[ind_length] ** 2

            # Reset the s2 value if below the cutoff
            if k_val < k_cutoff[ind_length]:
                s2_aa[ind_k_val, ind_length] = fa ** 2
                s2_bb[ind_k_val, ind_length] = (1 - fa) ** 2
                s2_ab[ind_k_val, ind_length] = fa * (1 - fa)

    return s2_aa, s2_ab, s2_bb


def s2_wlc_randcopoly(k_val_vector, length_mono_kuhn, num_mono, fa, lam=0, dimensions=3, alpha_max=25):
    r"""
    s2_wlc_randcopoly - Evaluate the 2-point structure factor for the wormlike chain model
    for a random copolymer

    Parameters
    ----------
    k_val_vector : float (array)
        The value of the Fourier vector magnitude :math:`K`
    length_mono_kuhn : float (array)
        The length of a monomer unit in Kuhn lengths
    num_mono : float
        The number of monomer units in the chain
    fa : float
        Fraction of A monomers in the random copolymer
    lam : float
        Statistical correlation parameter for random copolymer
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    s2 : float (vector)
        Structure factor for the wormlike chain model for every k_val in k_val_vector

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """
    if type(length_mono_kuhn) == float or type(length_mono_kuhn) == int:
        length_mono_kuhn = np.array([length_mono_kuhn])

    s2_aa = np.zeros((len(k_val_vector), len(length_mono_kuhn)), dtype=type(1+1j))
    s2_ab = np.zeros((len(k_val_vector), len(length_mono_kuhn)), dtype=type(1+1j))
    s2_bb = np.zeros((len(k_val_vector), len(length_mono_kuhn)), dtype=type(1+1j))

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_mono_kuhn, dimensions=3)
    k_cutoff_factor = 0.01
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    for ind_k_val in range(0, len(k_val_vector)):
        k_val = k_val_vector[ind_k_val]
        lam_max = alpha_max

        poles = eval_poles(k_val, 0, dimensions, alpha_max)
        residues = eval_residues(k_val, 0, poles, True, dimensions, lam_max, alpha_max, k_val_cutoff=10 ** -3)

        for ind_length in range(0, len(length_mono_kuhn)):
            for alpha in range(0, alpha_max):
                z0 = np.exp(poles[alpha] * length_mono_kuhn[ind_length])
                z1 = z0 * lam
                valeq = (2 * num_mono * (z0 / poles[alpha] ** 2 - 1 / poles[alpha] ** 2 -
                                         length_mono_kuhn[ind_length] / poles[alpha]))
                valne1 = (4 / poles[alpha] ** 2 * z0 *
                          (z0 ** num_mono - num_mono * z0 + num_mono - 1) /
                          (1 - z0) ** 2 * (np.cosh(poles[alpha] * length_mono_kuhn[ind_length]) - 1))
                valne2 = (4 / poles[alpha] ** 2 * z1 *
                          (z1 ** num_mono - num_mono * z1 + num_mono - 1) /
                          (1 - z1) ** 2 * (np.cosh(poles[alpha] * length_mono_kuhn[ind_length]) - 1))

                s2_aa[ind_k_val, ind_length] += (residues[alpha] * (fa * valeq +
                                                                    fa ** 2 * valne1 + fa * (1 - fa) * valne2))
                s2_ab[ind_k_val, ind_length] += (residues[alpha] * (fa * (1 - fa) * valne1 -
                                                                    fa * (1 - fa) * valne2))
                s2_bb[ind_k_val, ind_length] += (residues[alpha] * ((1 - fa) * valeq +
                                                                    (1 - fa) ** 2 * valne1 + fa * (1 - fa) * valne2))

            s2_aa[ind_k_val, ind_length] /= (length_mono_kuhn[ind_length] * num_mono) ** 2
            s2_ab[ind_k_val, ind_length] /= (length_mono_kuhn[ind_length] * num_mono) ** 2
            s2_bb[ind_k_val, ind_length] /= (length_mono_kuhn[ind_length] * num_mono) ** 2

            # Reset the s2 value if below the cutoff
#            if k_val < k_cutoff[ind_length]:
#                s2_aa[ind_k_val, ind_length] = fa ** 2
#                s2_ab[ind_k_val, ind_length] = fa * (1 - fa)
#                s2_bb[ind_k_val, ind_length] = (1 - fa) ** 2

    return s2_aa, s2_ab, s2_bb


# Cubic order density functions

def s3_wlc_diblock(k1_val_vector, k2_val_vector, length_kuhn, fa, dimensions=3, alpha_max=25):
    r"""
    s3_wlc - Evaluate the 3-point structure factor for the wormlike chain model

    Parameters
    ----------
    k1_val_vector : float (array)
        The value of the Fourier vector :math:`\vec{K}_{1}`
    k2_val_vector : float (array)
        The value of the Fourier vector :math:`\vec{K}_{2}`
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    fa : float
        Fraction of A monomers in the random copolymer
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    s3_wlc : float (vector)
        3-point Structure factor for the wormlike chain model for every k1_val and k2_val
        in k1_val_vector and k2_val_vector

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    # Assign the size of s3_wlc based on whether multiple values of k1 and length are input
    if type(length_kuhn) == float or type(length_kuhn) == int:
        length_kuhn = np.array([length_kuhn])
    s3_aaa = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))
    s3_aab = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))
    s3_abb = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))
    s3_bbb = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_kuhn, dimensions=3)
    k_cutoff_factor_zero = 0.2
    k_cutoff_zero = k_cutoff_factor_zero * np.sqrt(1 / rg2)
    k_cutoff_factor = 0.5
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    # Loop through the k values

    for ind_k_val in range(0, np.size(k1_val_vector, axis=0)):
        # Calculate the scattering vector magnitudes
        k1_val = k1_val_vector[ind_k_val, :]    # k1 vector
        k2_val = k2_val_vector[ind_k_val, :]    # k2 vector
        k12_val = k1_val + k2_val               # k1 + k2 vector
        k1_mag = np.linalg.norm(k1_val)         # k1 magnitude
        k2_mag = np.linalg.norm(k2_val)         # k2 magnitude
        k12_mag = np.linalg.norm(k12_val)       # k12 magnitude

        # Evaluate the crossover fraction for the two algorithms
        k_max = np.max([k1_mag, k2_mag, k12_mag])
        frac_algo_zero = 1 - np.exp(- k_max / k_cutoff)

        # Evaluate the angular arguments for the Legendre polynomials

        rho_12 = -np.dot(k1_val, k2_val) / (k1_mag * k2_mag)        # Cos angle between k1 and k2
        rho_112 = np.dot(k1_val, k12_val) / (k1_mag * k12_mag)      # Cos angle between k1 and k2
        rho_212 = np.dot(k2_val, k12_val) / (k2_mag * k12_mag)      # Cos angle between k1 and k2

        # Calculate the Legendre polynomials
        p_lam_12 = eval_legendre(rho_12, 0, alpha_max)
        p_lam_112 = eval_legendre(rho_112, 0, alpha_max)
        p_lam_212 = eval_legendre(rho_212, 0, alpha_max)

        # Evaluate the poles and residues of the three vectors
        lam_zero_only = False
        lam_max = alpha_max
        poles_k1 = eval_poles(k1_mag, 0, dimensions, alpha_max=alpha_max)
        residues_k1 = eval_residues(k1_mag, 0, poles_k1, lam_zero_only, dimensions,
                                    lam_max, alpha_max, k_val_cutoff=2e-4)
        residue_zero_k1, ddp_residue_zero_k1 = eval_residue_zero(k1_mag, 0, lam_zero_only, lam_max, dimensions)

        poles_k2 = eval_poles(k2_mag, 0, dimensions, alpha_max=alpha_max)
        residues_k2 = eval_residues(k2_mag, 0, poles_k2, lam_zero_only, dimensions,
                                    lam_max, alpha_max, k_val_cutoff=2e-4)
        residue_zero_k2, ddp_residue_zero_k2 = eval_residue_zero(k2_mag, 0, lam_zero_only, lam_max, dimensions)

        poles_k12 = eval_poles(k12_mag, 0, dimensions, alpha_max=alpha_max)
        residues_k12 = eval_residues(k12_mag, 0, poles_k12, lam_zero_only, dimensions,
                                     lam_max, alpha_max, k_val_cutoff=2e-4)
        residue_zero_k12, ddp_residue_zero_k12 = eval_residue_zero(k12_mag, 0, lam_zero_only, lam_max, dimensions)

        # Add the contribution from the p=0 double pole

        # Case 1: s1 < s2 < s3
        for ind_length in range(0, len(length_kuhn)):
            s3_aaa[ind_k_val, ind_length] += frac_algo_zero[ind_length] * 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + ddp_residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + residue_zero_k1[:, 0] * ddp_residue_zero_k12[:, 0] * p_lam_112)

        # Case 2: s1 < s3 < s2
        for ind_length in range(0, len(length_kuhn)):
            s3_aaa[ind_k_val, ind_length] += frac_algo_zero[ind_length] * 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                + ddp_residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                + residue_zero_k1[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_12)

        # Case 3: s2 < s1 < s3
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] += frac_algo_zero[ind_length] * 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                + ddp_residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                + residue_zero_k12[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_212)

        # Sum over the poles for the wlc green functions
        for alpha in range(0, alpha_max + 1):
            for alpha_p in range(0, alpha_max + 1):

            # Case 1: s1 < s2 < s3
                poles_vec = np.array([poles_k1[alpha], poles_k12[alpha_p], 0, 0])
                integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=(1 - frac_algo_zero))
                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k1[:, 0, alpha] * residues_k12[:, 0, alpha_p] * p_lam_112) * integ

            # Case 2: s1 < s3 < s2
                poles_vec = np.array([poles_k1[alpha], poles_k2[alpha_p], 0, 0])
                integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=(1 - frac_algo_zero))
                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k1[:, 0, alpha] * residues_k2[:, 0, alpha_p] * p_lam_12) * integ

            # Case 3: s2 < s1 < s3
                poles_vec = np.array([poles_k2[alpha], poles_k12[alpha_p], 0, 0])
                integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=(1 - frac_algo_zero))
                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k2[:, 0, alpha] * residues_k12[:, 0, alpha_p] * p_lam_212) * integ

        # Rescale the values by the factor 1 / length_kuhn ** 3
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 3.

            # Reset the s2 value if below the cutoff
            if min(k1_mag, k2_mag) < k_cutoff_zero[ind_length]:
                s3[ind_k_val, ind_length] = 1

    return s3