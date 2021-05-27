
from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np

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
    s2_wlc : float (vector)
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
    k_cutoff_factor = 0.3
    k_cutoff = k_cutoff_factor * np.sqrt(1 / rg2)

    # Loop through the k values

    for ind_k_val in range(0, np.size(k1_val_vector, axis=0)):
        k1_val = k1_val_vector[ind_k_val, :]    # k1 vector
        k2_val = k2_val_vector[ind_k_val, :]    # k2 vector
        k12_val = k1_val + k2_val               # k1 + k2 vector
        k1_mag = np.linalg.norm(k1_val)         # k1 magnitude
        k2_mag = np.linalg.norm(k2_val)         # k2 magnitude
        k12_mag = np.linalg.norm(k12_val)       # k12 magnitude

        # Reset the k vectors if any are equal (need to correct double pole case)
        if k1_mag == k2_mag or k1_mag == k12_mag or k2_mag == k12_mag:
            k2_val = 1.01 * k2_val               # k2 vector
            k12_val = k1_val + k2_val               # k1 + k2 vector
            k2_mag = np.linalg.norm(k2_val)         # k2 magnitude
            k12_mag = np.linalg.norm(k12_val)       # k12 magnitude

        rho_12 = -np.dot(k1_val, k2_val) / (k1_mag * k2_mag)        # Cos angle between k1 and k2
        rho_112 = np.dot(k1_val, k12_val) / (k1_mag * k12_mag)      # Cos angle between k1 and k2
        rho_212 = np.dot(k2_val, k12_val) / (k2_mag * k12_mag)      # Cos angle between k1 and k2

        # Calculate the Legendre polynomials
        p_lam_12 = eval_legendre(rho_12, alpha_max)
        p_lam_112 = eval_legendre(rho_112, alpha_max)
        p_lam_212 = eval_legendre(rho_212, alpha_max)

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
            s3[ind_k_val, ind_length] += 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + ddp_residue_zero_k1[:, 0] * residue_zero_k12[:, 0] * p_lam_112
                + residue_zero_k1[:, 0] * ddp_residue_zero_k12[:, 0] * p_lam_112)

        # Case 2: s1 < s3 < s2
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] += 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                + ddp_residue_zero_k1[:, 0] * residue_zero_k2[:, 0] * p_lam_12
                + residue_zero_k1[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_12)

        # Case 3: s2 < s1 < s3
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] += 2 * np.sum(
                length_kuhn[ind_length] * residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                + ddp_residue_zero_k12[:, 0] * residue_zero_k2[:, 0] * p_lam_212
                + residue_zero_k12[:, 0] * ddp_residue_zero_k2[:, 0] * p_lam_212)

        # Calculate the residues for the poles (case 1 - double pole, case 2 - unequal poles)
        lam_max = alpha_max

        # Case 1: s1 < s2 < s3
        if k1_mag == k12_mag:
            residues_double_k1 = eval_residues_double_pole(k1_mag, 0, poles_k1,
                                                        lam_zero_only, lam_max, alpha_max, dimensions)
        else:
            residues_other_pole_k112 = eval_residues_other_pole(k1_mag, 0, poles_k12,
                                                              lam_zero_only, lam_max, alpha_max, dimensions)
            residues_other_pole_k121 = eval_residues_other_pole(k12_mag, 0, poles_k1,
                                                               lam_zero_only, lam_max, alpha_max, dimensions)

        # Case 2: s1 < s3 < s2
        if k1_mag == k2_mag:
            residues_double_k1 = eval_residues_double_pole(k1_mag, 0, poles_k1,
                                                        lam_zero_only, lam_max, alpha_max, dimensions)
        else:
            residues_other_pole_k12 = eval_residues_other_pole(k1_mag, 0, poles_k2,
                                                              lam_zero_only, lam_max, alpha_max, dimensions)
            residues_other_pole_k21 = eval_residues_other_pole(k2_mag, 0, poles_k1,
                                                              lam_zero_only, lam_max, alpha_max, dimensions)

        # Case 3: s2 < s1 < s3
        if k2_mag == k12_mag:
            residues_double_k2 = eval_residues_double_pole(k2_mag, 0, poles_k2,
                                                        lam_zero_only, lam_max, alpha_max, dimensions)
        else:
            residues_other_pole_k212 = eval_residues_other_pole(k2_mag, 0, poles_k12,
                                                              lam_zero_only, lam_max, alpha_max, dimensions)
            residues_other_pole_k122 = eval_residues_other_pole(k12_mag, 0, poles_k2,
                                                               lam_zero_only, lam_max, alpha_max, dimensions)

        # Sum over the poles for the wlc green functions
        for alpha in range(0, alpha_max + 1):
            # Case 1: s1 < s2 < s3
            if k1_mag == k12_mag:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_double_k1[:, 0, alpha] * p_lam_112 * np.exp(poles_k1[alpha] * length_kuhn[ind_length])
                        + residues_k1[:, 0, alpha] ** 2 * p_lam_112 * length_kuhn[ind_length] *
                        np.exp(poles_k1[alpha] * length_kuhn[ind_length])) / poles_k1[alpha] ** 2
            else:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_k1[:, 0, alpha] * residues_other_pole_k121[:, 0, alpha] * p_lam_112) * (
                            np.exp(poles_k1[alpha] * length_kuhn[ind_length]) / poles_k1[alpha] ** 2) + 2 * np.sum(
                        residues_k12[:, 0, alpha] * residues_other_pole_k112[:, 0, alpha] * p_lam_112) * (
                            np.exp(poles_k12[alpha] * length_kuhn[ind_length]) / poles_k12[alpha] ** 2)

            # Case 2: s1 < s3 < s2
            if k1_mag == k2_mag:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_double_k1[:, 0, alpha] * p_lam_12 * np.exp(poles_k1[alpha] * length_kuhn[ind_length])
                        + residues_k1[:, 0, alpha] ** 2 * p_lam_12 * length_kuhn[ind_length] *
                        np.exp(poles_k1[alpha] * length_kuhn[ind_length])) / poles_k1[alpha] ** 2
            else:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_k1[:, 0, alpha] * residues_other_pole_k21[:, 0, alpha] * p_lam_12) * (
                            np.exp(poles_k1[alpha] * length_kuhn[ind_length]) / poles_k1[alpha] ** 2) + 2 * np.sum(
                        residues_k2[:, 0, alpha] * residues_other_pole_k12[:, 0, alpha] * p_lam_12) * (
                            np.exp(poles_k2[alpha] * length_kuhn[ind_length]) / poles_k2[alpha] ** 2)

            # Case 3: s2 < s1 < s3
            if k2_mag == k12_mag:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_double_k2[:, 0, alpha] * p_lam_212 * np.exp(poles_k2[alpha] * length_kuhn[ind_length])
                        + residues_k2[:, 0, alpha] ** 2 * p_lam_212 * length_kuhn[ind_length] *
                        np.exp(poles_k2[alpha] * length_kuhn[ind_length])) / poles_k2[alpha] ** 2
            else:
                for ind_length in range(0, len(length_kuhn)):
                    s3[ind_k_val, ind_length] += 2 * np.sum(
                        residues_k2[:, 0, alpha] * residues_other_pole_k122[:, 0, alpha] * p_lam_212) * (
                            np.exp(poles_k2[alpha] * length_kuhn[ind_length]) / poles_k2[alpha] ** 2) + 2 * np.sum(
                        residues_k12[:, 0, alpha] * residues_other_pole_k212[:, 0, alpha] * p_lam_212) * (
                            np.exp(poles_k12[alpha] * length_kuhn[ind_length]) / poles_k12[alpha] ** 2)

        # Rescale the values by the factor 1 / length_kuhn ** 3
        for ind_length in range(0, len(length_kuhn)):
            s3[ind_k_val, ind_length] /= length_kuhn[ind_length] ** 3.

            # Reset the s2 value if below the cutoff
            if min(k1_mag, k2_mag) < k_cutoff[ind_length]:
                s3[ind_k_val, ind_length] = 1

    return s3


def eval_legendre(rho, alpha_max=25):
    r"""
    Evaluate the vector of Legendre polynomial values

    Parameters
    ----------
    rho : float
        Cos theta of the angle
    alpha_max : int
        Maximum number of polynomials evaluated (default 25)

    Returns
    -------
    legendre_poly : float
        Vector of Legendre polynomials

    """

    legendre_poly = np.zeros((alpha_max + 1), dtype=float)

    # Define the first two Legendre polynomials
    legendre_poly[0] = 1.
    legendre_poly[1] = rho

    # Define the remaining Legendre polynomials based on the recurrence relation
    for i_l in range(2, alpha_max + 1):
        legendre_poly[i_l] = (rho * legendre_poly[i_l-1] * (2 * i_l - 1) / i_l
                              - legendre_poly[i_l - 2] * (i_l - 1) / i_l)

#    for i_l in range(0, alpha_max + 1):
#        legendre_poly[i_l] *= 1 / (2 * i_l + 1)

    return legendre_poly



