
from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np
import scipy as sp


# Quadratic order density functions

def s2_wlc(k_val_vector, length_kuhn, dimensions=3, alpha_max=25):
    r"""
    s2_wlc - Evaluate the 2-point structure factor for the wormlike chain model

    Parameters
    ----------
    k_val_vector : float (array)
        The value of the Fourier vector magnitude :math:`K`
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
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

    s2 = np.zeros((len(k_val_vector), len(length_kuhn)), dtype=type(1+1j))

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
            s2[ind_k_val, ind_length] = (length_kuhn[ind_length] * residue_zero + ddp_residue_zero)

            for alpha in range(0, alpha_max):
                s2[ind_k_val, ind_length] += (np.exp(poles[alpha] * length_kuhn[ind_length]) *
                                                  residues[alpha] / poles[alpha] ** 2)

            s2[ind_k_val, ind_length] *= 2 / length_kuhn[ind_length] ** 2

            # Reset the s2 value if below the cutoff
            if k_val < k_cutoff[ind_length]:
                s2[ind_k_val, ind_length] = 1

    return s2


# Cubic order density functions

def s3_wlc(k1_val_vector, k2_val_vector, length_kuhn, dimensions=3, alpha_max=25):
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
    s3 = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))

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

        rho_12 = - np.dot(k1_val, k2_val) / (k1_mag * k2_mag)        # Cos angle between k1 and k2
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
            s3[ind_k_val, ind_length] += frac_algo_zero[ind_length] * 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + ddp_residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + residue_zero_k1[:, 0] * ddp_residue_zero_k12[:, 0] * p_lam_112)

        # Case 2: s1 < s3 < s2
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] += frac_algo_zero[ind_length] * 2 * np.sum(
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


# Quartic order density functions

def s4_wlc(k1_val_vector, k2_val_vector, k3_val_vector, length_kuhn, dimensions=3, alpha_max=25):
    r"""
    s4_wlc - Evaluate the 4-point structure factor for the wormlike chain model

    Parameters
    ----------
    k1_val_vector : float (array)
        The value of the Fourier vector :math:`\vec{K}_{1}`
    k2_val_vector : float (array)
        The value of the Fourier vector :math:`\vec{K}_{2}`
    k3_val_vector : float (array)
        The value of the Fourier vector :math:`\vec{K}_{2}`
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    s4_wlc : float (vector)
        4-point Structure factor for the wormlike chain model for every k1_val, k2_val, k3_val
        in k1_val_vector, k2_val_vector, k3_val_vector

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    # Assign the size of s4_wlc based on whether multiple values of k1 and length are input
    if type(length_kuhn) == float or type(length_kuhn) == int:
        length_kuhn = np.array([length_kuhn])
    s4 = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1+1j))

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_kuhn, dimensions=3)
    k_cutoff_factor = 0.3
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_kuhn, dimensions=3)
    k_cutoff_factor_zero = 0.2
    k_cutoff_zero = k_cutoff_factor_zero * np.sqrt(1 / rg2)
    k_cutoff_factor = 1
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    # Loop through the k values

    for ind_k_val in range(0, np.size(k1_val_vector, axis=0)):

        k123_val = np.zeros((3, 3))
        k123_val[0, :] = k1_val_vector[ind_k_val, :]    # k1 vector
        k123_val[1, :] = k2_val_vector[ind_k_val, :]    # k2 vector
        k123_val[2, :] = k3_val_vector[ind_k_val, :]    # k3 vector

        # Setup the ka, kb, kc vectors for each case of s ordering
        ka_vector = np.zeros((12, 3))
        kb_vector = np.zeros((12, 3))
        kc_vector = np.zeros((12, 3))
        for i_c in range(6):
            ka_vector[i_c, :] = k123_val[0, :] + k123_val[1, :] + k123_val[2, :]
        kb_vector[0, :] = k123_val[0, :] + k123_val[1, :]
        kb_vector[1, :] = k123_val[0, :] + k123_val[1, :]
        kb_vector[2, :] = k123_val[0, :] + k123_val[2, :]
        kb_vector[3, :] = k123_val[1, :] + k123_val[2, :]
        kb_vector[4, :] = k123_val[0, :] + k123_val[2, :]
        kb_vector[5, :] = k123_val[1, :] + k123_val[2, :]
        kc_vector[0, :] = k123_val[0, :]
        kc_vector[1, :] = k123_val[1, :]
        kc_vector[2, :] = k123_val[0, :]
        kc_vector[3, :] = k123_val[1, :]
        kc_vector[4, :] = k123_val[2, :]
        kc_vector[5, :] = k123_val[2, :]
        ka_vector[6, :] = - k123_val[2, :]
        ka_vector[7, :] = - k123_val[2, :]
        ka_vector[8, :] = - k123_val[1, :]
        ka_vector[9, :] = - k123_val[0, :]
        ka_vector[10, :] = - k123_val[1, :]
        ka_vector[11, :] = - k123_val[0, :]
        kb_vector[6, :] = k123_val[0, :] + k123_val[1, :]
        kb_vector[7, :] = k123_val[0, :] + k123_val[1, :]
        kb_vector[8, :] = k123_val[0, :] + k123_val[2, :]
        kb_vector[9, :] = k123_val[1, :] + k123_val[2, :]
        kb_vector[10, :] = k123_val[0, :] + k123_val[2, :]
        kb_vector[11, :] = k123_val[1, :] + k123_val[2, :]
        kc_vector[6, :] = k123_val[0, :]
        kc_vector[7, :] = k123_val[1, :]
        kc_vector[8, :] = k123_val[0, :]
        kc_vector[9, :] = k123_val[1, :]
        kc_vector[10, :] = k123_val[2, :]
        kc_vector[11, :] = k123_val[2, :]

        # Cycle through the cases

        for i_c in range(12):
            ka_val = ka_vector[i_c, :]              # ka vector
            kb_val = kb_vector[i_c, :]              # kb vector
            kc_val = kc_vector[i_c, :]              # kc vector

            ka_mag = np.linalg.norm(ka_val)         # ka magnitude
            kb_mag = np.linalg.norm(kb_val)         # kb magnitude
            kc_mag = np.linalg.norm(kc_val)         # kc magnitude

            # Evaluate the crossover fraction for the two algorithms
            k_max = np.max([ka_mag, kb_mag, kc_mag])
            frac_algo_zero = 1 - np.exp(- k_max / k_cutoff)

            # Reset the k vectors if any are equal (need to correct double pole case)
            # if ka_mag == kb_mag and ka_mag == kc_mag:
            #     ka_val = 0.99 * ka_val              # Reset ka vector
            #     kc_val = 1.01 * kc_val              # Reset kc vector
            #     ka_mag = np.linalg.norm(ka_val)  # ka magnitude
            #     kc_mag = np.linalg.norm(kc_val)  # kc magnitude
            # elif ka_mag == kb_mag or kb_mag == kc_mag:
            #     kb_val = 1.01 * kb_val               # kb vector
            #     kb_mag = np.linalg.norm(kb_val)
            # elif ka_mag == kc_mag:
            #     kc_val = 1.01 * kc_val               # kb vector
            #     kc_mag = np.linalg.norm(kc_val)

            # Setup the polar and azimuthal angles between ka, kb, kc
            rho_ab = np.dot(ka_val, kb_val) / (ka_mag * kb_mag)     # $\cos \theta$ between ka and kb
            rho_bc = np.dot(kb_val, kc_val) / (kb_mag * kc_mag)     # $\cos \theta$ between kb and kc

            eb = kb_val / np.linalg.norm(kb_val)
            ka_perp = ka_val - np.dot(ka_val, eb) * eb
            ea_perp = ka_perp / np.linalg.norm(ka_perp)
            kc_perp = kc_val - np.dot(kc_val, eb) * eb
            ec_perp = kc_perp / np.linalg.norm(kc_perp)
            phi_abc = np.arctan2(np.dot(ec_perp, np.cross(ea_perp, eb)), np.dot(ec_perp, ea_perp))

            # Cycle through the mu values
            lam_max = alpha_max
            for mu in range(lam_max):

                # Evaluate the poles and residues of the three vectors
                lam_zero_only = False
                poles_ka = eval_poles(ka_mag, mu, dimensions, alpha_max=alpha_max)
                residues_ka = eval_residues(ka_mag, mu, poles_ka, lam_zero_only, dimensions,
                                            lam_max, alpha_max, k_val_cutoff=1e-2)
                residue_zero_ka, ddp_residue_zero_ka = eval_residue_zero(ka_mag, mu, lam_zero_only,
                                                                         lam_max, dimensions)

                poles_kb = eval_poles(kb_mag, mu, dimensions, alpha_max=alpha_max)
                residues_kb = eval_residues(kb_mag, mu, poles_kb, lam_zero_only, dimensions,
                                            lam_max, alpha_max, k_val_cutoff=1e-2, lam_cont_frac_max=50)
                residue_zero_kb, ddp_residue_zero_kb = eval_residue_zero(kb_mag, mu, lam_zero_only,
                                                                         lam_max, dimensions)

                poles_kc = eval_poles(kc_mag, mu, dimensions, alpha_max=alpha_max)
                residues_kc = eval_residues(kc_mag, mu, poles_kc, lam_zero_only, dimensions,
                                            lam_max, alpha_max, k_val_cutoff=1e-2, lam_cont_frac_max=50)
                residue_zero_kc, ddp_residue_zero_kc = eval_residue_zero(kc_mag, mu, lam_zero_only,
                                                                         lam_max, dimensions)

                # Evaluate the spherical harmonic matrix

                p_lam_ab = eval_legendre(rho_ab, mu, alpha_max)
                p_lam_bc = eval_legendre(rho_bc, mu, alpha_max)
                ylm_abc = np.outer(p_lam_ab, p_lam_bc) * np.exp(1j * mu * phi_abc)

                s4_mu = np.zeros((len(length_kuhn)), dtype=type(1+1j))

                # Add the contribution from the p=0 double pole
                for ind_length in range(0, len(length_kuhn)):
                    s4_mu[ind_length] += frac_algo_zero[ind_length] * np.sum(
                        length_kuhn[ind_length] * np.outer(residue_zero_ka[:, 0], residue_zero_kc[:, 0]) *
                        residue_zero_kb * ylm_abc +
                        np.outer(ddp_residue_zero_ka[:, 0], residue_zero_kc[:, 0]) *
                        residue_zero_kb * ylm_abc +
                        np.outer(residue_zero_ka[:, 0], ddp_residue_zero_kc[:, 0]) *
                        residue_zero_kb * ylm_abc +
                        np.outer(residue_zero_ka[:, 0], residue_zero_kc[:, 0]) *
                        ddp_residue_zero_kb * ylm_abc)

                # Sum over the poles for the wlc green functions

                for alpha_a in range(0, alpha_max - abs(mu) + 1):
                    for alpha_b in range(0, alpha_max - abs(mu) + 1):
                        for alpha_c in range(0, alpha_max - abs(mu) + 1):

                            # Calculate the contribution for the three poles
                            poles_vec = np.array([poles_ka[alpha_a], poles_kb[alpha_b], poles_ka[alpha_c], 0, 0])
                            integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=(1 - frac_algo_zero))
                            s4_mu += np.sum(np.outer(residues_ka[:, 0, alpha_a], residues_kc[:, 0, alpha_c]) *
                                            residues_kb[:, :, alpha_b] * ylm_abc) * integ

                # Add the contribute at this mu value to the 4-point function
                if mu == 0:
                    s4[ind_k_val, :] += 2 * s4_mu
                else:
                    s4[ind_k_val, :] += 2 * 2 * np.real(s4_mu)

        # Rescale the values by the factor 1 / length_kuhn ** 4
        for ind_length in range(0, len(length_kuhn)):
            s4[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 4.

            # Reset the s2 value if below the cutoff
            if min(np.linalg.norm(k123_val[0, :]),
                   np.linalg.norm(k123_val[1, :]),
                   np.linalg.norm(k123_val[2, :])) < k_cutoff_zero[ind_length]:
                s4[ind_k_val, ind_length] = 1

    return s4


def eval_legendre(rho, mu=0, alpha_max=25):
    r"""
    Evaluate the vector of Legendre polynomial values

    Parameters
    ----------
    rho : float
        Cos theta of the angle
    mu : int
        Value of the z-momentum quantum number
    alpha_max : int
        Maximum number of polynomials evaluated (default 25)

    Returns
    -------
    legendre_poly : float
        Vector of Legendre polynomials

    """

    mu = abs(mu)
    legendre_poly = np.zeros((alpha_max - abs(mu) + 1), dtype=float)

    # Define the first two Legendre polynomials
    legendre_poly[0] = (-1) ** mu * sp.special.factorial2(2 * mu - 1) * (1 - rho ** 2) ** (mu / 2)
    legendre_poly[1] = rho * (2 * mu + 1) * legendre_poly[0]

    # Define the remaining Legendre polynomials based on the recurrence relation
    for i_l in range(2, alpha_max - abs(mu) + 1):
        lam = i_l + abs(mu)
        legendre_poly[i_l] = (rho * legendre_poly[i_l - 1] * (2 * lam - 1)
                              - legendre_poly[i_l - 2] * (lam + mu - 1)) / (lam - mu)

    if mu != 0:
        for i_l in range(alpha_max - abs(mu) + 1):
            lam = i_l + abs(mu)
            legendre_poly[i_l] *= np.sqrt(sp.special.factorial(lam - mu) / sp.special.factorial(lam + mu))

    return legendre_poly


def calc_int_mag(length_kuhn, poles_vec, frac_zero=1):
    r"""
    Evaluate the magnitude of the integral for a list of poles (including repeats). This algorithm includes
    cases for single, double, and triple poles (as needed in evaluation of correlation functions)

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    poles_vec : float (array)
        Array of poles
    frac_zero : float
        Multiplicative factor for contributions from poles at p=0

    Returns
    -------
    int_mag : float
        Value of the integral over chain length for the five poles

    """

    # Determine the number of repeated poles and their order
    poles_vec_unique, poles_order = np.unique(poles_vec, return_counts=True)

    # Cycle through the poles and evaluate contribution for each order
    int_mag = np.zeros_like(length_kuhn, dtype=type(1 + 1j))

    for i_pole in range(len(poles_vec_unique)):
        # Rescale the contribution for zero poles
        if poles_vec_unique[i_pole] == 0:
            frac = frac_zero
        else:
            frac = 1

        # Algorithm for the case of a simple pole
        if poles_order[i_pole] == 1:
            f_poles = 1 / (np.prod((poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) **
                                   np.delete(poles_order, i_pole)))
            int_mag += frac * f_poles * np.exp(poles_vec_unique[i_pole] * length_kuhn)

        # Algorithm for the case of a double pole
        elif poles_order[i_pole] == 2:
            f_poles = 1 / (np.prod((poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) **
                                   np.delete(poles_order, i_pole)))
            prod_poles = np.sum(np.delete(poles_order, i_pole) /
                                (poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)))
            ddp_f_poles = - f_poles * prod_poles

            int_mag += frac * (length_kuhn * f_poles + ddp_f_poles) * np.exp(poles_vec_unique[i_pole] * length_kuhn)

        # Algorithm for the case of a triple pole
        elif poles_order[i_pole] == 3:
            f_poles = 1 / (np.prod((poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) **
                                   np.delete(poles_order, i_pole)))
            prod_poles = np.sum(np.delete(poles_order, i_pole) /
                                (poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)))
            prod_poles_2 = np.sum(np.delete(poles_order, i_pole) /
                                  (poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) ** 2)
            ddp_f_poles = - f_poles * prod_poles
            d2dp2_f_poles = f_poles * (prod_poles_2 + prod_poles ** 2)

            int_mag += frac * 0.5 * (length_kuhn ** 2 * f_poles + 2 * length_kuhn * ddp_f_poles +
                                     d2dp2_f_poles) * np.exp(poles_vec_unique[i_pole] * length_kuhn)

    return int_mag
