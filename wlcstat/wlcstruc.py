
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
    k_cutoff_factor = 1.0e-5
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    for ind_k_val in range(0, len(k_val_vector)):
        k_val = k_val_vector[ind_k_val]
        lam_max = alpha_max

        poles = eval_poles(k_val, 0, dimensions, alpha_max)
        residues = eval_residues(k_val, 0, poles, True, dimensions, lam_max, alpha_max, k_val_cutoff=10 ** -3)
        residue_zero, ddp_residue_zero = eval_residue_zero(k_val, 0, dimensions)

        for ind_length in range(0, len(length_kuhn)):

            # Use the p=0 residue if k > 10.0
            if k_val > 10.0:
                s2[ind_k_val, ind_length] = (length_kuhn[ind_length] * residue_zero + ddp_residue_zero)

                for alpha in range(0, alpha_max):
                    s2[ind_k_val, ind_length] += (np.exp(poles[alpha] * length_kuhn[ind_length]) *
                                                  residues[alpha] / poles[alpha] ** 2)

            # Use the direct integration of the Laplace inversed function if k <= 10.0
            else:
                s2[ind_k_val, ind_length] = 0.
                for alpha in range(0, alpha_max):
                    s2[ind_k_val, ind_length] += (- length_kuhn[ind_length] * poles[alpha] +
                                                  np.exp(poles[alpha] * length_kuhn[ind_length])
                                                  - 1.) * residues[alpha] / poles[alpha] ** 2

            s2[ind_k_val, ind_length] *= 2 / length_kuhn[ind_length] ** 2

            # Reset the s2 value if below the cutoff
            if k_val < k_cutoff[ind_length]:
                s2[ind_k_val, ind_length] = 1

    return s2


# Cubic order density functions

def s3_wlc(k1_val_vector, k2_val_vector, length_kuhn, dimensions=3, alpha_max=15):
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

    # Define the tolerance for poles being equal valued
    pole_tol = 1.e-5

    # Loop through the k values

    for ind_k_val in range(0, np.size(k1_val_vector, axis=0)):
        # Calculate the scattering vector magnitudes
        k1_val = k1_val_vector[ind_k_val, :]    # k1 vector
        k2_val = k2_val_vector[ind_k_val, :]    # k2 vector
        k12_val = k1_val + k2_val               # k1 + k2 vector
        k1_mag = np.linalg.norm(k1_val)         # k1 magnitude
        k2_mag = np.linalg.norm(k2_val)         # k2 magnitude
        k12_mag = np.linalg.norm(k12_val)       # k12 magnitude

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

        # Determine the algorithm for the 3 cases

        if np.sqrt(k1_mag ** 2. + k12_mag ** 2.) < 10.:
            case1_zero = True
        else:
            case1_zero = False

        if np.sqrt(k1_mag ** 2. + k2_mag ** 2.) < 10.:
            case2_zero = True
        else:
            case2_zero = False

        if np.sqrt(k2_mag ** 2. + k12_mag ** 2.) < 10.:
            case3_zero = True
        else:
            case3_zero = False

        # Add the contribution from the p=0 double pole

        # Case 1: s1 < s2 < s3
        for ind_length in range(0, len(length_kuhn)):
            if not case1_zero:
                s3[ind_k_val, ind_length] += 2 * np.sum(
                    length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                    + ddp_residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                    + residue_zero_k1[:, 0] * ddp_residue_zero_k12[:, 0] * p_lam_112)

        # Case 2: s1 < s3 < s2
        for ind_length in range(0, len(length_kuhn)):
            if not case2_zero:
                s3[ind_k_val, ind_length] += 2 * np.sum(
                    length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                    + ddp_residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                    + residue_zero_k1[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_12)

        # Case 3: s2 < s1 < s3
        for ind_length in range(0, len(length_kuhn)):
            if not case3_zero:
                s3[ind_k_val, ind_length] += 2 * np.sum(
                    length_kuhn[ind_length] * residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                    + ddp_residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                    + residue_zero_k12[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_212)

        # Sum over the poles for the wlc green functions
        for alpha in range(0, alpha_max + 1):
            for alpha_p in range(0, alpha_max + 1):

            # Case 1: s1 < s2 < s3
                pole1 = poles_k1[alpha]
                pole2 = poles_k12[alpha_p]
                if np.abs(pole1 - pole2) < pole_tol:
                    pole2 = pole1
                poles_vec = np.array([pole1, pole2, 0, 0])
                if case1_zero:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=1.0)
                else:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=0.0)

                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k1[:, 0, alpha] * residues_k12[:, 0, alpha_p] * p_lam_112) * integ

            # Case 2: s1 < s3 < s2
                pole1 = poles_k1[alpha]
                pole2 = poles_k2[alpha_p]
                if np.abs(pole1 - pole2) < pole_tol:
                    pole2 = pole1
                poles_vec = np.array([pole1, pole2, 0, 0])
                if case2_zero:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=1.0)
                else:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=0.0)

                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k1[:, 0, alpha] * residues_k2[:, 0, alpha_p] * p_lam_12) * integ

            # Case 3: s2 < s1 < s3
                pole1 = poles_k2[alpha]
                pole2 = poles_k12[alpha_p]
                if np.abs(pole1 - pole2) < pole_tol:
                    pole2 = pole1
                poles_vec = np.array([pole1, pole2, 0, 0])
                if case3_zero:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=1.0)
                else:
                    integ = calc_int_mag(length_kuhn, poles_vec, frac_zero=0.0)

                s3[ind_k_val, :] += 2 * np.sum(
                    residues_k2[:, 0, alpha] * residues_k12[:, 0, alpha_p] * p_lam_212) * integ

        # Rescale the values by the factor 1 / length_kuhn ** 3
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 3.

            # Reset the s2 value if below the cutoff
            if min(k1_mag, k2_mag) < k_cutoff_zero[ind_length]:
                s3[ind_k_val, ind_length] = 1.

    return s3


# Quartic order density functions

def s4_wlc(k1_val_vector, k2_val_vector, k3_val_vector, length_kuhn, dimensions=3, alpha_max=10):
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
    s4 = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1 + 1j))
    s4_0 = np.zeros((np.size(k1_val_vector, axis=0), len(length_kuhn)), dtype=type(1 + 1j))

    # Calculate the radius of gyration to determine a cutoff value of k
    rg2 = rg2_ave(length_kuhn, dimensions=3)
    k_cutoff_factor_zero = 0.2
    k_cutoff_zero = k_cutoff_factor_zero * np.sqrt(1 / rg2)

    # Define the tolerance for poles being equal valued
    pole_tol = 1.e-3

    # Define the conditions for switching to the zero pole algorithm
#    k_case_zero = 100.
    k_case_zero = 40.

    # Loop through the k values

    for ind_k_val in range(0, np.size(k1_val_vector, axis=0)):

        k123_val = np.zeros((3, 3))
        k123_val[0, :] = k1_val_vector[ind_k_val, :]  # k1 vector
        k123_val[1, :] = k2_val_vector[ind_k_val, :]  # k2 vector
        k123_val[2, :] = k3_val_vector[ind_k_val, :]  # k3 vector

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

        # Evaluate the residues for the 12 cases
        mu = 0
        poles_ka = np.zeros((alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))
        residues_ka = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12),
                               dtype=type(1 + 1j))
        residue_zero_ka = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))
        ddp_residue_zero_ka = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))

        poles_kc = np.zeros((alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))
        residues_kc = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12),
                               dtype=type(1 + 1j))
        residue_zero_kc = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))
        ddp_residue_zero_kc = np.zeros((alpha_max - abs(mu) + 1, alpha_max - abs(mu) + 1, 12), dtype=type(1 + 1j))
        rho_ab = np.zeros((12))
        rho_bc = np.zeros((12))
        phi_abc = np.zeros((12))
        case_zero = np.zeros((12), dtype=bool)

        # Cycle through the cases

        for i_c in range(12):
            ka_val = ka_vector[i_c, :]  # ka vector
            kb_val = kb_vector[i_c, :]  # kb vector
            kc_val = kc_vector[i_c, :]  # kc vector

            ka_mag = np.linalg.norm(ka_val)  # ka magnitude
            kb_mag = np.linalg.norm(kb_val)  # kb magnitude
            kc_mag = np.linalg.norm(kc_val)  # kc magnitude

            # Determine whether to use the zero-pole algorithm
            if np.sqrt(ka_mag ** 2. + kb_mag ** 2. + kc_mag ** 2.) < k_case_zero:
                case_zero[i_c] = True
            else:
                case_zero[i_c] = False

            # Setup the polar and azimuthal angles between ka, kb, kc
            rho_ab[i_c] = np.dot(ka_val, kb_val) / (ka_mag * kb_mag)  # $\cos \theta$ between ka and kb
            rho_bc[i_c] = np.dot(kb_val, kc_val) / (kb_mag * kc_mag)  # $\cos \theta$ between kb and kc

            if kb_mag == 0:
                eb = np.array([0, 0, 1])
            else:
                eb = kb_val / np.linalg.norm(kb_val)

            ka_perp = ka_val - np.dot(ka_val, eb) * eb
            ea_perp = ka_perp / np.linalg.norm(ka_perp)
            kc_perp = kc_val - np.dot(kc_val, eb) * eb
            ec_perp = kc_perp / np.linalg.norm(kc_perp)
            phi_abc[i_c] = np.arctan2(np.dot(ec_perp, np.cross(ea_perp, eb)), np.dot(ec_perp, ea_perp))

            # Evaluate the poles and residues of vectors a and c (mu = 0)
            lam_zero_only = False
            lam_max = alpha_max
            mu = 0
            poles_ka[:, i_c] = eval_poles(ka_mag, mu, dimensions, alpha_max=alpha_max)
            residues_ka[:, :, :, i_c] = eval_residues(ka_mag, mu, poles_ka[:, i_c], lam_zero_only, dimensions,
                                        lam_max, alpha_max, k_val_cutoff=1e-2)
            residue_zero_ka[:, :, i_c], ddp_residue_zero_ka[:, :, i_c] = eval_residue_zero(ka_mag, mu, lam_zero_only,
                                                                                           lam_max, dimensions)

            poles_kc[:, i_c] = eval_poles(kc_mag, mu, dimensions, alpha_max=alpha_max)
            residues_kc[:, :, :, i_c] = eval_residues(kc_mag, mu, poles_kc[:, i_c], lam_zero_only, dimensions,
                                        lam_max, alpha_max, k_val_cutoff=1e-2)
            residue_zero_kc[:, :, i_c], ddp_residue_zero_kc[:, :, i_c] = eval_residue_zero(kc_mag, mu, lam_zero_only,
                                                                                           lam_max, dimensions)

        # Evaluate the poles and residues for vector b at each mu
        poles_kb = np.zeros((alpha_max + 1, alpha_max + 1, 12), dtype=type(1 + 1j))
        residues_kb = np.zeros((alpha_max + 1, alpha_max + 1, alpha_max + 1, alpha_max + 1, 12),
                               dtype=type(1 + 1j))
        residue_zero_kb = np.zeros((alpha_max + 1, alpha_max + 1, alpha_max + 1, 12), dtype=type(1 + 1j))
        ddp_residue_zero_kb = np.zeros((alpha_max + 1, alpha_max + 1, alpha_max + 1, 12),
                                       dtype=type(1 + 1j))

        ylm_abc = np.zeros((alpha_max + 1, alpha_max + 1, alpha_max + 1, 12),
                           dtype=type(1 + 1j))

        # Cycle through the mu values
        for mu in range(alpha_max + 1):
            for i_c in range(12):
                kb_val = kb_vector[i_c, :]  # kb vector
                kb_mag = np.linalg.norm(kb_val)  # kb magnitude

                # Evaluate the poles and residues of the vector c at each mu value
                lam_zero_only = False
                poles_kb[mu:, mu, i_c] = eval_poles(kb_mag, mu, dimensions, alpha_max=alpha_max)
                residues_kb[mu:, mu:, mu:, mu, i_c] = \
                    eval_residues(kb_mag, mu, poles_kb[mu:, mu, i_c], lam_zero_only,
                                  dimensions,lam_max, alpha_max=alpha_max)
                residue_zero_kb[mu:, mu:, mu, i_c], ddp_residue_zero_kb[mu:, mu:, mu, i_c] = \
                    eval_residue_zero(kb_mag, mu, lam_zero_only, lam_max, dimensions)

                # Evaluate the spherical harmonic matrix

                p_lam_ab = eval_legendre(rho_ab[i_c], mu, alpha_max)
                p_lam_bc = eval_legendre(rho_bc[i_c], mu, alpha_max)
                ylm_abc[mu:, mu:, mu, i_c] = np.outer(p_lam_ab, p_lam_bc) * np.cos(mu * phi_abc[i_c])

        # Calculate the contributions for each case to s4

        # Add the contribution from the p=0 double pole
        s4_zero = np.zeros((len(length_kuhn)), dtype=type(1 + 1j))
        for mu in range(alpha_max + 1):
            if mu == 0:
                coef = 2.
            else:
                coef = 4.
            for i_c in range(12):
                #                if not case_zero[i_c]:
                s4_zero += coef * (length_kuhn * np.sum(
                        np.outer(residue_zero_ka[mu:(alpha_max + 1), mu, i_c],
                                 residue_zero_kc[mu:(alpha_max + 1), mu, i_c]) *
                        residue_zero_kb[mu:, mu:, mu, i_c] * ylm_abc[mu:, mu:, mu, i_c])
                                + np.sum(np.outer(ddp_residue_zero_ka[mu:(alpha_max + 1), mu, i_c],
                                                  residue_zero_kc[mu:(alpha_max + 1), mu, i_c]) *
                                         residue_zero_kb[mu:, mu:, mu, i_c] * ylm_abc[mu:, mu:, mu, i_c])
                                + np.sum(np.outer(residue_zero_ka[mu:(alpha_max + 1), mu, i_c],
                                                  ddp_residue_zero_kc[mu:(alpha_max + 1), mu, i_c]) *
                                         residue_zero_kb[mu:, mu:, mu, i_c] * ylm_abc [mu:, mu:, mu, i_c])
                                + np.sum(np.outer(residue_zero_ka[mu:(alpha_max + 1), mu, i_c],
                                                  residue_zero_kc[mu:(alpha_max + 1), mu, i_c]) *
                                         ddp_residue_zero_kb[mu:, mu:, mu, i_c] * ylm_abc[mu:, mu:, mu, i_c]))

        s4_0[ind_k_val, :] += s4_zero

        # Sum over the poles for the wlc green functions

        for alpha_a in reversed(range(0, alpha_max + 1)):
            for alpha_c in reversed(range(0, alpha_max + 1)):
                for alpha_b in reversed(range(0, alpha_max + 1)):
                    s4_mu = np.zeros((len(length_kuhn)), dtype=type(1 + 1j))
                    s4_mu_0 = np.zeros((len(length_kuhn)), dtype=type(1 + 1j))
                    for mu in reversed(range(0, alpha_b + 1)):
                        if mu == 0:
                            coef = 2.
                        else:
                            coef = 4.
                        for i_c in range(12):
                            # Calculate the contribution for the three poles
                            pole1 = poles_ka[alpha_a, i_c]
                            pole2 = poles_kb[alpha_b, mu, i_c]
                            pole3 = poles_kc[alpha_c, i_c]
                            # Set poles to be equal if below the tolerance
                            if np.abs(pole1 - pole2) < pole_tol:
                                pole2 = pole1
                            if np.abs(pole2 - pole3) < pole_tol:
                                pole3 = pole2
                            if np.abs(pole1 - pole3) < pole_tol:
                                pole3 = pole1
#                            if case_zero[i_c]:
#                                integ = calc_s4_int_mag(length_kuhn, pole1, pole2, pole3, frac_zero=1.0)
#                            else:
#                                integ = calc_s4_int_mag(length_kuhn, pole1, pole2, pole3, frac_zero=0.0)

                            integ = calc_s4_int_mag(length_kuhn, pole1, pole2, pole3, frac_zero=1.0)
                            s4_mu += coef * np.sum(
                                np.outer(residues_ka[mu:, 0, alpha_a, i_c],
                                         residues_kc[mu:, 0, alpha_c, i_c]) *
                                residues_kb[mu:, mu:, alpha_b, mu, i_c] * ylm_abc[mu:, mu:, mu, i_c]) * integ

                            integ = calc_s4_int_mag(length_kuhn, pole1, pole2, pole3, frac_zero=0.0)
                            s4_mu += coef * np.sum(
                                np.outer(residues_ka[mu:, 0, alpha_a, i_c],
                                         residues_kc[mu:, 0, alpha_c, i_c]) *
                                residues_kb[mu:, mu:, alpha_b, mu, i_c] * ylm_abc[mu:, mu:, mu, i_c]) * integ
                    # Add the contribute at this k value to the 4-point function
                    s4[ind_k_val, :] += s4_mu
                    s4_0[ind_k_val, :] += s4_mu_0

        # Rescale the values by the factor 1 / length_kuhn ** 4
        for ind_length in range(0, len(length_kuhn)):
            s4[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 4.
            s4_0[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 4.

            # Reset the s2 value if below the cutoff
            if min(np.linalg.norm(k123_val[0, :]),
                   np.linalg.norm(k123_val[1, :]),
                   np.linalg.norm(k123_val[2, :])) < k_cutoff_zero[ind_length]:
                s4[ind_k_val, ind_length] = 1.

    s4 = np.real(s4)

    return s4, s4_0


# Auxiliary functions for the calculations

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
    if alpha_max - abs(mu) >= 1:
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


def calc_int_mag(length_kuhn, poles_vec, frac_zero=1.0):
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
            frac = 1.0

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


def calc_s4_int_mag(length_kuhn, pole1, pole2, pole3, frac_zero=1.0):
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

    # Cycle through the poles and evaluate contribution for each order
    int_mag = np.zeros_like(length_kuhn, dtype=type(1 + 1j))

    if pole1 == pole2 and pole2 == pole3:
        if pole1 == 0:
            int_mag = length_kuhn ** 4. / 24.
        else:
            int_mag = (((pole1 * length_kuhn) ** 2. - 4. * pole1 * length_kuhn + 6.) * np.exp(pole1 * length_kuhn)
                       + frac_zero * (-2. * pole1 * length_kuhn - 6.)) / (2. * pole1 ** 4.)
    elif pole1 == pole2 and pole2 != pole3:
        int_mag = ((2. * pole3 - 3. * pole1) * np.exp(pole1 * length_kuhn) / (pole1 ** 3. * (pole1 - pole3) ** 2.)
                   + length_kuhn * np.exp(pole1 * length_kuhn) / (pole1 ** 2. * (pole1 - pole3))
                   + np.exp(pole3 * length_kuhn) / (pole3 ** 2. * (pole3 - pole1) ** 2.)
                   + frac_zero * (-length_kuhn / (pole1 ** 2. * pole3)
                                  -(pole1 + 2. * pole3) / (pole1 ** 3. * pole3 ** 2.)))
    elif pole1 == pole3 and pole2 != pole3:
        int_mag = ((2. * pole2 - 3. * pole1) * np.exp(pole1 * length_kuhn) / (pole1 ** 3. * (pole1 - pole2) ** 2.)
                   + length_kuhn * np.exp(pole1 * length_kuhn) / (pole1 ** 2. * (pole1 - pole2))
                   + np.exp(pole2 * length_kuhn) / (pole2 ** 2. * (pole2 - pole1) ** 2.)
                   + frac_zero * (-length_kuhn / (pole1 ** 2. * pole2)
                                  -(pole1 + 2. * pole2) / (pole1 ** 3. * pole2 ** 2.)))
    elif pole2 == pole3 and pole1 != pole3:
        int_mag = ((2. * pole1 - 3. * pole2) * np.exp(pole2 * length_kuhn) / (pole2 ** 3. * (pole2 - pole1) ** 2.)
                   + length_kuhn * np.exp(pole2 * length_kuhn) / (pole2 ** 2. * (pole2 - pole1))
                   + np.exp(pole1 * length_kuhn) / (pole1 ** 2. * (pole1 - pole2) ** 2.)
                   + frac_zero * (-length_kuhn / (pole2 ** 2. * pole1)
                                  -(pole2 + 2. * pole1) / (pole2 ** 3. * pole1 ** 2.)))
    else:
        int_mag = (np.exp(pole1 * length_kuhn) / (pole1 ** 2. * (pole1 - pole2) * (pole1 - pole3))
                   + np.exp(pole2 * length_kuhn) / (pole2 ** 2. * (pole2 - pole1) * (pole2 - pole3))
                   + np.exp(pole3 * length_kuhn) / (pole3 ** 2. * (pole3 - pole1) * (pole3 - pole2))
                   + frac_zero * (-length_kuhn / (pole1 * pole2 * pole3)
                                  -1. / (pole1 * pole2 * pole3 ** 2.)
                                  - 1. / (pole1 * pole3 * pole2 ** 2.)
                                  - 1. / (pole2 * pole3 * pole1 ** 2.)))

    return int_mag
