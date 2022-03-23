"""
wlc_lcpoly

Module containing functions for evaluating the
self-consistent field theory for polymer solutions, including nematic
liquid crystallinity

"""
import numpy as np
import scipy as sp
from numba import jit
import matplotlib.pyplot as plt


def q_lcpoly(length_kuhn, lam, alpha_max=25):
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
        q_val = np.zeros_like(length_kuhn, dtype=type(1+1j))
        q_val += 1
    elif lam < 50 and lam != 0:
        poles = eval_poles_lcpoly(lam, 0, alpha_max)
        residues = eval_residues_lcpoly(lam, 0, poles, l_zero_only=True, l_max=alpha_max,
                                        alpha_max=alpha_max, l_cont_frac_max=50)
        q_val = 0
        for ind_alpha in range(0, alpha_max + 1):
            q_val += np.exp(poles[ind_alpha] * length_kuhn) * residues[ind_alpha]

    # Evaluate the partition function using the large-lam expansion
    else:
        q_val = np.exp(length_kuhn * (2 / 3 * lam - 2 * np.sqrt(lam) + 1 + 1 / 4 / np.sqrt(lam)
                                      + 1 / 4 / lam + 23 / 64 / lam ** (3 / 2)
                                      + 41 / 64 / lam ** 2 + 681 / 512 / lam ** (5 / 2)))
    return q_val


def m_lcpoly(length_kuhn, lam, alpha_max=25, l_cont_frac_max=50):
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

    if lam < 1000:

        # Find the poles and residues for the calculation of q and m
        poles = eval_poles_lcpoly(lam, 0, alpha_max)
        residues = eval_residues_lcpoly(lam, 0, poles, l_zero_only=False, l_max=alpha_max,
                                        alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)

        # Reset poles by subtracting zeroth pole
        poles -= poles[0]

        # Setup the l_selection matrix
        a_l_m = eval_a_l_m(np.arange(0, alpha_max + 2), 0)

        l_eq_l0 = 1.5 * (a_l_m[0:(alpha_max + 1)] ** 2 + a_l_m[1:(alpha_max + 2)] ** 2)
        l_eq_l0p2 = 1.5 * (a_l_m[1:alpha_max] * a_l_m[2:(alpha_max + 1)])
        p20_select_mat = np.diag(l_eq_l0) + np.diag(l_eq_l0p2, 2) + np.diag(l_eq_l0p2, -2)

        # Calculate the partition function and m value
        q_val = np.zeros_like(length_kuhn, dtype=type(1+1j))
        m_val = np.zeros_like(length_kuhn, dtype=type(1+1j))
        for ind_alpha in range(0, alpha_max + 1):
            q_val += np.exp(poles[ind_alpha] * length_kuhn) * residues[0, 0, ind_alpha]
            for ind_alpha_p in range(0, alpha_max + 1):
                select_mag = np.real(np.dot(residues[:, :, ind_alpha_p],
                                            np.dot(p20_select_mat, residues[:, :, ind_alpha]))[0, 0])
                if ind_alpha_p == ind_alpha:
                    int_mag = np.exp(poles[ind_alpha] * length_kuhn)
                else:
                    int_mag = ((np.exp(poles[ind_alpha] * length_kuhn) - np.exp(poles[ind_alpha_p] * length_kuhn)) /
                               (poles[ind_alpha] - poles[ind_alpha_p])) / length_kuhn
                m_val += select_mag * int_mag
        m_val = np.real(m_val / q_val) - 1 / 2

    # Asymptotic solution based on spheroidal harmonic ground-state analysis
    else:
        m_val = (1 - 3 / 2 / np.sqrt(lam) - 3 / 16 / lam ** (3 / 2) - 3 / 8 / lam ** 2
                 - 207 / 256 / lam ** (5 / 2) - 123 / 64 / lam ** 3 - 10215 / 2048 / lam ** (7 / 2))

    return m_val


def r_2_lcpoly(length_kuhn, lam, alpha_max=25, l_cont_frac_max=50):
    r"""
    Calculate the mean-square end-to-end distance for a liquid crystal polymer

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
    r_2_par : float
        mean-square distance in the parallel direction
    r_2_perp : float
        mean-square distance in the perpendicular direction
    xi_par : float
        parallel correlation length
    xi_perp : float
        perpendicular correlation length

    """

    poles_m0 = eval_poles_lcpoly(lam, m=0, alpha_max=alpha_max)
    resi_m0 = eval_residues_lcpoly(lam, m=0, poles=poles_m0, l_zero_only=False, l_max=alpha_max,
                                   alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)
    poles_m1 = eval_poles_lcpoly(lam, m=1, alpha_max=alpha_max)
    resi_m1 = eval_residues_lcpoly(lam, m=1, poles=poles_m1, l_zero_only=False, l_max=alpha_max,
                                   alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)

    # Reset poles by subtracting zeroth pole
    max_pole = poles_m0[0]
    poles_m0 -= max_pole
    poles_m1 -= max_pole

    # Setup the l_selection matrices
    a_l_m0 = np.zeros(alpha_max + 2)
    a_l_m0[1::] = eval_a_l_m(np.arange(1, alpha_max + 2), 0)
    b_l_m0 = np.zeros(alpha_max + 2)
    b_l_m0[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), 0)
    b_l_m1 = np.zeros(alpha_max + 2)
    b_l_m1[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), 1)

    l_eq_l0p1 = a_l_m0[1:(alpha_max + 1)]
    uz_select_mat = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0p1, -1)

    l_eq_l0p1 = -1 / np.sqrt(2) * b_l_m0[1:(alpha_max + 1)]
    l_eq_l0m1 = 1 / np.sqrt(2) * b_l_m1[1:(alpha_max + 1)]
    ux_select_mat = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0m1, -1)

    # Calculate the partition function and m value
    q_val = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    r_2_par = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    r_2_perp = np.zeros_like(length_kuhn, dtype=type(1 + 1j))

    for ind_alpha in range(0, alpha_max + 1, 2):
        # Contributions to q_val
        q_val += np.exp(poles_m0[ind_alpha] * length_kuhn) * resi_m0[0, 0, ind_alpha]
        for ind_alpha_p in range(1, alpha_max + 1, 2):
            for ind_alpha_pp in range(0, alpha_max + 1, 2):
                # Contributions to the r_2_par average
                if ind_alpha_p <= alpha_max:
                    select_mag = np.real(
                        np.dot(resi_m0[:, :, ind_alpha_pp], np.dot(np.transpose(uz_select_mat),
                        np.dot(resi_m0[:, :, ind_alpha_p], np.dot(uz_select_mat,
                        resi_m0[:, :, ind_alpha]))))[0, 0])
                    poles_vec = np.array([poles_m0[ind_alpha], poles_m0[ind_alpha_p], poles_m0[ind_alpha_pp]])
                    int_mag = calc_int_mag(length_kuhn, poles_vec)
                    r_2_par += 2 * select_mag * int_mag
                # Contributions to the r_2_perp average
                if ind_alpha_p <= (alpha_max - 1):
                    select_mag = np.real(
                        np.dot(resi_m0[:, :, ind_alpha_pp],  np.dot(np.transpose(ux_select_mat),
                        np.dot(resi_m1[:, :, ind_alpha_p - 1], np.dot(ux_select_mat,
                        resi_m0[:, :, ind_alpha]))))[0, 0])
                    poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p - 1], poles_m0[ind_alpha_pp]])
                    int_mag = calc_int_mag(length_kuhn, poles_vec)
                    r_2_perp += 2 * select_mag * int_mag

    # Finalize the averages
    q_val = np.real(q_val)
    r_2_par = np.real(r_2_par / q_val)
    r_2_perp = np.real(r_2_perp / q_val)
    xi_par = - 1 / poles_m0[1]
    xi_perp = - 1 / poles_m1[1]

    return r_2_par, r_2_perp, xi_par, xi_perp


def elastic_lcpoly(length_kuhn, lam, alpha_max=25, l_cont_frac_max=200):
    r"""
    Calculate the Frank elastic constants for a polymer liquid crystal solution

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    lam : float
        The value of the quadrupole field :math:`\lambda`
    alpha_max : int
        Maximum number of poles evaluated (default 50)
    l_cont_frac_max : int
        Number of levels included in the evaluation of residues from continued fraction

    Returns
    -------
    q_val : float
        Single-chain partition function
    m_val : float
        Nematic order parameter
    y21_y21 : float
        y21-y21 correlation function
    y21_ux_ux_y21 : float
        y21-y21 correlation function with x-x end-to-end distance squared
    y21_uy_uy_y21 : float
        y21-y21 correlation function with y-y end-to-end distance squared
    y21_uz_uz_y21 : float
        y21-y21 correlation function with z-z end-to-end distance squared

    """

    poles_m0 = eval_poles_lcpoly(lam, m=0, alpha_max=alpha_max)
    resi_m0 = eval_residues_lcpoly(lam, m=0, poles=poles_m0, l_zero_only=False, l_max=alpha_max,
                                   alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)
    poles_m1 = eval_poles_lcpoly(lam, m=1, alpha_max=alpha_max)
    resi_m1 = eval_residues_lcpoly(lam, m=1, poles=poles_m1, l_zero_only=False, l_max=alpha_max,
                                   alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)
    poles_m2 = eval_poles_lcpoly(lam, m=2, alpha_max=alpha_max)
    resi_m2 = eval_residues_lcpoly(lam, m=2, poles=poles_m2, l_zero_only=False, l_max=alpha_max,
                                   alpha_max=alpha_max, l_cont_frac_max=l_cont_frac_max)

    # Reset poles by subtracting zeroth pole
    max_pole = poles_m0[0]
    poles_m0 -= max_pole
    poles_m1 -= max_pole
    poles_m2 -= max_pole

    # Setup the l_selection matrices
    a_l_m0 = np.zeros(alpha_max + 2)
    a_l_m0[1::] = eval_a_l_m(np.arange(1, alpha_max + 2), 0)
    a_l_m1 = np.zeros(alpha_max + 2)
    a_l_m1[1::] = eval_a_l_m(np.arange(1, alpha_max + 2), 1)
    b_l_m0 = np.zeros(alpha_max + 2)
    b_l_m0[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), 0)
    b_l_m1 = np.zeros(alpha_max + 2)
    b_l_m1[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), 1)
    b_l_mm1 = np.zeros(alpha_max + 2)
    b_l_mm1[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), -1)
    b_l_m2 = np.zeros(alpha_max + 2)
    b_l_m2[1::] = eval_b_l_m(np.arange(1, alpha_max + 2), 2)

    l_eq_l0 = a_l_m0[0:(alpha_max + 1)] ** 2 + a_l_m0[1:(alpha_max + 2)] ** 2
    l_eq_l0p2 = a_l_m0[1:alpha_max] * a_l_m0[2:(alpha_max + 1)]
    y20_select_mat = np.diag(l_eq_l0) + np.diag(l_eq_l0p2, 2) + np.diag(l_eq_l0p2, -2)

    l_eq_l0 = (b_l_m1[1:(alpha_max + 2)] * a_l_m1[1:(alpha_max + 2)]
               - b_l_m0[0:(alpha_max + 1)] * a_l_m1[0:(alpha_max + 1)])
    l_eq_l0m2 = b_l_m1[1:alpha_max] * a_l_m1[2:(alpha_max + 1)]
    l_eq_l0p2 = - b_l_m0[2:(alpha_max + 1)] * a_l_m1[1:alpha_max]
    y2l_select_mat = np.diag(l_eq_l0) + np.diag(l_eq_l0p2, 2) + np.diag(l_eq_l0m2, -2)

    l_eq_l0p1 = a_l_m1[1:(alpha_max + 1)]
    uz_select_mat = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0p1, -1)

    l_eq_l0m1 = 0.5 * b_l_m2[1:(alpha_max + 1)]
    l_eq_l0p1 = -0.5 * b_l_mm1[1:(alpha_max + 1)]
    ux_select_mat_m2 = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0m1, -1)
    uy_select_mat_m2 = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0m1, -1)

    l_eq_l0m1 = -1 / np.sqrt(2) * b_l_m0[1:(alpha_max + 1)]
    l_eq_l0p1 = 1 / np.sqrt(2) * b_l_m1[1:(alpha_max + 1)]
    ux_select_mat_m0 = np.diag(l_eq_l0p1, 1) + np.diag(l_eq_l0m1, -1)

    # Calculate the partition function and m value
    q_val = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    m_val = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    y21_y21 = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    y21_uz_uz_y21 = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    y21_ux_ux_y21 = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    y21_uy_uy_y21 = np.zeros_like(length_kuhn, dtype=type(1 + 1j))
    for ind_alpha in range(0, alpha_max + 1):
        # Contributions to q_val
        q_val += np.exp(poles_m0[ind_alpha] * length_kuhn) * resi_m0[0, 0, ind_alpha]
        for ind_alpha_p in range(0, alpha_max + 1):
            # Contributions to m_val
            select_mag = np.real(np.dot(resi_m0[:, :, ind_alpha_p],
                                        np.dot(y20_select_mat, resi_m0[:, :, ind_alpha]))[0, 0])
            poles_vec = np.array([poles_m0[ind_alpha], poles_m0[ind_alpha_p]])
            int_mag = calc_int_mag(length_kuhn, poles_vec)
            m_val += select_mag * int_mag / length_kuhn

            for ind_alpha_pp in range(0, alpha_max + 1):
                # Contributions to the y21_y21 average
                if ind_alpha_p <= (alpha_max - 1):
                    select_mag = np.real(
                        np.dot(resi_m0[:, :, ind_alpha_pp], np.dot(np.transpose(y2l_select_mat),
                        np.dot(resi_m1[:, :, ind_alpha_p], np.dot(y2l_select_mat,
                        resi_m0[:, :, ind_alpha]))))[0, 0])
                    poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p], poles_m0[ind_alpha_pp]])
                    int_mag = calc_int_mag(length_kuhn, poles_vec)
                    y21_y21 += 2 * select_mag * int_mag / length_kuhn ** 2

                for ind_alpha_ppp in range(0, alpha_max + 1):
                    for ind_alpha_pppp in range(0, alpha_max + 1):

                        # Contributions to the y21_uz_uz_y21 average
                        if (ind_alpha_p <= (alpha_max - 1) and ind_alpha_pp <= (alpha_max - 1)
                                and ind_alpha_ppp <= (alpha_max - 1)):
                            select_mag = np.real(
                                np.dot(resi_m0[:, :, ind_alpha_pppp], np.dot(np.transpose(y2l_select_mat),
                                np.dot(resi_m1[:, :, ind_alpha_ppp], np.dot(np.transpose(uz_select_mat),
                                np.dot(resi_m1[:, :, ind_alpha_pp], np.dot(uz_select_mat,
                                np.dot(resi_m1[:, :, ind_alpha_p], np.dot(y2l_select_mat,
                                resi_m0[:, :, ind_alpha]))))))))[0, 0])
                            poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p],
                                                  poles_m1[ind_alpha_pp], poles_m1[ind_alpha_ppp],
                                                  poles_m0[ind_alpha_pppp]])
                            int_mag = calc_int_mag(length_kuhn, poles_vec)
                            y21_uz_uz_y21 += 4 * select_mag * int_mag / length_kuhn ** 2

                        # Contributions to the y21_uy_uy_y21 average
                        if (ind_alpha_p <= (alpha_max - 1) and ind_alpha_pp <= (alpha_max - 2)
                                and ind_alpha_ppp <= (alpha_max - 1)):
                            select_mag = np.real(
                                np.dot(resi_m0[:, :, ind_alpha_pppp], np.dot(np.transpose(y2l_select_mat),
                                np.dot(resi_m1[:, :, ind_alpha_ppp], np.dot(np.transpose(uy_select_mat_m2),
                                np.dot(resi_m2[:, :, ind_alpha_pp], np.dot(uy_select_mat_m2,
                                np.dot(resi_m1[:, :, ind_alpha_p], np.dot(y2l_select_mat,
                                resi_m0[:, :, ind_alpha]))))))))[0, 0])
                            poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p],
                                                  poles_m2[ind_alpha_pp], poles_m1[ind_alpha_ppp],
                                                  poles_m0[ind_alpha_pppp]])
                            int_mag = calc_int_mag(length_kuhn, poles_vec)
                            y21_uy_uy_y21 += 4 * select_mag * int_mag / length_kuhn ** 2

                        # Contributions to the y21_ux_ux_y21 average
                        if (ind_alpha_p <= (alpha_max - 1) and ind_alpha_pp <= (alpha_max - 2)
                                and ind_alpha_ppp <= (alpha_max - 1)):
                            select_mag = np.real(
                                np.dot(resi_m0[:, :, ind_alpha_pppp], np.dot(np.transpose(y2l_select_mat),
                                np.dot(resi_m1[:, :, ind_alpha_ppp], np.dot(np.transpose(ux_select_mat_m2),
                                np.dot(resi_m2[:, :, ind_alpha_pp], np.dot(ux_select_mat_m2,
                                np.dot(resi_m1[:, :, ind_alpha_p], np.dot(y2l_select_mat,
                                resi_m0[:, :, ind_alpha]))))))))[0, 0])
                            poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p],
                                                  poles_m2[ind_alpha_pp], poles_m1[ind_alpha_ppp],
                                                  poles_m0[ind_alpha_pppp]])
                            int_mag = calc_int_mag(length_kuhn, poles_vec)
                            y21_ux_ux_y21 += 4 * select_mag * int_mag / length_kuhn ** 2

                        if (ind_alpha_p <= (alpha_max - 1) and ind_alpha_ppp <= (alpha_max - 1)):
                            select_mag = np.real(
                                np.dot(resi_m0[:, :, ind_alpha_pppp], np.dot(np.transpose(y2l_select_mat),
                                np.dot(resi_m1[:, :, ind_alpha_ppp], np.dot(np.transpose(ux_select_mat_m0),
                                np.dot(resi_m0[:, :, ind_alpha_pp], np.dot(ux_select_mat_m0,
                                np.dot(resi_m1[:, :, ind_alpha_p], np.dot(y2l_select_mat,
                                resi_m0[:, :, ind_alpha]))))))))[0, 0])
                            poles_vec = np.array([poles_m0[ind_alpha], poles_m1[ind_alpha_p],
                                                  poles_m0[ind_alpha_pp], poles_m1[ind_alpha_ppp],
                                                  poles_m0[ind_alpha_pppp]])
                            int_mag = calc_int_mag(length_kuhn, poles_vec)
                            y21_ux_ux_y21 += 4 * select_mag * int_mag / length_kuhn ** 2

    # Finalize the averages
    q_val = np.real(q_val)
    m_val = np.real(1.5 * m_val / q_val) - 1 / 2
    y21_y21 = np.real(y21_y21 / q_val) * 15 / (8 * np.pi)
    y21_ux_ux_y21 = np.real(y21_ux_ux_y21 / q_val) * 15 / (8 * np.pi)
    y21_uy_uy_y21 = np.real(y21_uy_uy_y21 / q_val) * 15 / (8 * np.pi)
    y21_uz_uz_y21 = np.real(y21_uz_uz_y21 / q_val) * 15 / (8 * np.pi)

    return q_val, m_val, y21_y21, y21_ux_ux_y21, y21_uy_uy_y21, y21_uz_uz_y21


def elastic_rr(length_kuhn, lam):
    r"""
    Calculate the Frank elastic constants for a polymer liquid crystal solution

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    lam : float
        The value of the quadrupole field :math:`\lambda`

    Returns
    -------
    q_val : float
        Single-chain partition function
    m_val : float
        Nematic order parameter
    y21_y21 : float
        y21-y21 correlation function
    y21_ux_ux_y21 : float
        y21-y21 correlation function with x-x end-to-end distance squared
    y21_uy_uy_y21 : float
        y21-y21 correlation function with y-y end-to-end distance squared
    y21_uz_uz_y21 : float
        y21-y21 correlation function with z-z end-to-end distance squared

    """

    a_val = np.outer(lam, length_kuhn)

    # Calculate the single-chain partition function
    z_val = np.sqrt(a_val)
    q_val = np.sqrt(np.pi) * sp.special.erfi(z_val) * np.exp(-a_val) / z_val
    q_val[np.isnan(q_val)] = 1 / a_val[np.isnan(q_val)]     # Reset to the asymptotic form for instances of NaN

    # Calculate rho averages
    rho2_ave = (1 / a_val / q_val - 1 / (2 * a_val))
    rho4_ave = ((2 * a_val - 3) / (2 * a_val ** 2) / q_val + 3 / (4 * a_val ** 2))
    rho6_ave = ((4 * a_val ** 2 - 10 * a_val + 15) / (4 * a_val ** 3) / q_val
                - 15 / (8 * a_val ** 3))

    # Determine the order parameter and the correlation functions
    n2_mat = np.outer(np.ones(len(lam)), length_kuhn ** 2) / 6

    m_val = 3 * rho2_ave / 2 - 1 / 2
    y21_y21 = 15 / np.pi * (rho2_ave - rho4_ave) / 2
    y21_ux_ux_y21 = n2_mat * 15 / np.pi * (rho2_ave - 2 * rho4_ave + rho6_ave) * 3 / 8
    y21_uy_uy_y21 = y21_ux_ux_y21 / 3
    y21_uz_uz_y21 = n2_mat * 15 / np.pi * (rho4_ave - rho6_ave) / 2

    return q_val, m_val, y21_y21, y21_ux_ux_y21, y21_uy_uy_y21, y21_uz_uz_y21


def eval_poles_lcpoly(lam, m=0, alpha_max=25):
    r"""
    eval_poles_lcpoly - Evaluate the poles for given :math:`\lambda` and :math:`\mu`
    using the matrix method for intermediate :math:`K`

    Parameters
    ----------
    lam : float
        The value of the nematic field :math:`\lambda`
    m : int
        Value of the mu parameter
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    poles : float
        Evaluated poles for the given :math:`\lambda` and :math:`\mu`

    Notes
    -----

    """

    # Determine the number of poles evaluated noting conjugate pairing
    # Evaluate 4 times the number needed to achieve accuracy

    num_total = int(4 * 2 * (np.ceil(alpha_max / 2.0)))  # Eigenvalues come in pairs

    # Build the h-matrix used to evaluate the poles of the Green function
    l = m + np.arange(0, num_total)
    a_l = eval_a_l_m(l, m)
    a_lp1 = eval_a_l_m(l + 1, m)
    a_lp2 = eval_a_l_m(l + 2, m)

    diagonal = l * (l + 1) - lam * (a_lp1 ** 2 + a_l ** 2 - 1 / 3)
    diagonal_plus2 = - lam * a_lp1[0:(num_total - 2)] * a_lp2[0:(num_total - 2)]

    h_matrix = np.diag(diagonal) + np.diag(diagonal_plus2, 2) + np.diag(diagonal_plus2, -2)

    # Find the poles as the eigenvalues of the h-matrix
    poles_total = -1 * np.linalg.eigvals(h_matrix)
    poles_total = np.sort(poles_total)[::-1]

    poles = poles_total[0:(alpha_max - abs(m) + 1)]

    # Reset poles using asymptotic behavior for Spheroidal harmonic eigenvalues (see Meixner and Schafke, page 319)

    for ind in range(0, (alpha_max - abs(m)), 2):
        l_val = ind + abs(m)
        if abs(poles[ind] - poles[ind + 1]) < 1e-10:
            q = l_val + 1
            p = int(0.5 * (q - m - 1))
            delta_poles = np.exp(- 2 * np.sqrt(lam)) * (
                    2 * (4 * np.sqrt(lam)) ** (q + 1) / (sp.math.factorial(p) * sp.math.factorial(m + p))) * (
                    1 + (m ** 2 - (q + 1) * (3 * q + 1)) / (8 * np.sqrt(lam)))
            poles[ind + 1] = poles[ind] - delta_poles

    return poles


def eval_residues_lcpoly(lam, m, poles, l_zero_only=True, l_max=25, alpha_max=25, l_cont_frac_max=50):
    r"""
    eval_residues_lcpoly -
    Evaluate the residues for the Green's function of a nematic polymer

    Parameters
    ----------
    lam : float
        The value of the nematic field :math:`\lambda`
    m : int
        Value of the mu parameter
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu`
    l_zero_only : boolean
        Indicates whether the residues will be evaluated over the range of :math:`\lambda` and :math:`lambda_{0}`
    l_max : int
        Maximum lambda value evaluated
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    l_cont_frac_max : int
        Maximum :math:`\lambda` value in the continued fraction evaluation

    Returns
    -------
    residues : complex float
        Evaluated residues for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2008]_ for intermediate-k algorithms

    """

    # Initialize the residue array based on whether lam_zero_only

    l_cont_frac_max = l_cont_frac_max + np.mod(l_cont_frac_max, 2) + np.mod(m, 2)
    l_max = max(l_max, alpha_max)

    if l_zero_only:
        residues = np.zeros((alpha_max - abs(m) + 1), dtype=type(1+1j))
    else:
        residues = np.zeros((l_max + 1, l_max + 1, alpha_max - abs(m) + 1), dtype=type(1+1j))

    for alpha in range(abs(m), alpha_max + 1):
        ind_alpha = alpha - abs(m)

        # Build the continued fractions
        j_plus = np.zeros((l_cont_frac_max - abs(m) + 1), dtype=type(1+1j))
        djdp_plus = np.zeros((l_cont_frac_max - abs(m) + 1), dtype=type(1+1j))
        j_minus = np.zeros((l_cont_frac_max - abs(m) + 1), dtype=type(1+1j))
        djdp_minus = np.zeros((l_cont_frac_max - abs(m) + 1), dtype=type(1+1j))
        a_l_m = eval_a_l_m(np.arange(abs(m), l_cont_frac_max + 1), m)

        # l contributions to the residue matrix with same parity as m

        l = l_cont_frac_max
        ind_l = l - abs(m)
        j_plus[ind_l] = poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 - 1 / 3)
        djdp_plus[ind_l] = 1

        for l in reversed(range(abs(m), l_cont_frac_max - 1, 2)):
            ind_l = l - abs(m)
            j_plus[ind_l] = (poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
                             - (a_l_m[ind_l + 1] * a_l_m[ind_l + 2] * lam) ** 2 / j_plus[ind_l + 2])
            djdp_plus[ind_l] = (1 + (a_l_m[ind_l + 1] * a_l_m[ind_l + 2] * lam /
                                     j_plus[ind_l + 2]) ** 2 * djdp_plus[ind_l + 2])

        l = abs(m)
        ind_l = l - abs(m)
        j_minus[0] = poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l + 1] ** 2 - 1 / 3)
        djdp_minus[0] = 1
        for l in range(abs(m) + 2, max(l_max, alpha_max) + 1, 2):
            ind_l = l - abs(m)
            j_minus[ind_l] = (poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
                              - (a_l_m[ind_l] * a_l_m[ind_l - 1] * lam) ** 2 / j_minus[ind_l - 2])
            djdp_minus[ind_l] = 1 + (a_l_m[ind_l] * a_l_m[ind_l - 1] * lam /
                                     j_minus[ind_l - 2]) ** 2 * djdp_minus[ind_l - 2]

        # l contributions to the residue matrix with opposite parity as m

        l = l_cont_frac_max - 1
        ind_l = l - abs(m)
        j_plus[ind_l] = poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
        djdp_plus[ind_l] = 1

        for l in reversed(range(abs(m) + 1, l_cont_frac_max - 2, 2)):
            ind_l = l - abs(m)
            j_plus[ind_l] = (poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
                             - (a_l_m[ind_l + 1] * a_l_m[ind_l + 2] * lam) ** 2 / j_plus[ind_l + 2])
            djdp_plus[ind_l] = (1 + (a_l_m[ind_l + 1] * a_l_m[ind_l + 2] * lam /
                                     j_plus[ind_l + 2]) ** 2 * djdp_plus[ind_l + 2])

        l = abs(m) + 1
        ind_l = l - abs(m)
        j_minus[1] = poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
        djdp_minus[1] = 1
        for l in range(abs(m) + 3, max(l_max, alpha_max), 2):
            ind_l = l - abs(m)
            j_minus[ind_l] = (poles[ind_alpha] + l * (l + 1) - lam * (a_l_m[ind_l] ** 2 + a_l_m[ind_l + 1] ** 2 - 1 / 3)
                              - (a_l_m[ind_l] * a_l_m[ind_l - 1] * lam) ** 2 / j_minus[ind_l - 2])
            djdp_minus[ind_l] = 1 + (a_l_m[ind_l] * a_l_m[ind_l - 1] * lam /
                                     j_minus[ind_l - 2]) ** 2 * djdp_minus[ind_l - 2]

        # Construct the residue matrix

        if l_zero_only:
            if np.mod(ind_alpha, 2) == 0:
                if ind_alpha == 0:
                    residues[ind_alpha] = 1 / djdp_plus[0]
                else:
                    w_alpha = 1 / (djdp_plus[ind_alpha]
                                   + (a_l_m[ind_alpha] * a_l_m[ind_alpha - 1] * lam /
                                      j_minus[ind_alpha - 2]) ** 2 * djdp_minus[ind_alpha - 2])
                    w_prod_left = np.prod(lam * a_l_m[1:(ind_alpha + 1):2] * a_l_m[2:(ind_alpha + 2):2] /
                                          j_minus[0:ind_alpha:2])
                    residues[ind_alpha] = w_prod_left ** 2 * w_alpha
        else:
            if np.mod(ind_alpha, 2) == 0:
                if ind_alpha == 0:
                    w_alpha = 1 / djdp_plus[0]
                else:
                    w_alpha = 1 / (djdp_plus[ind_alpha] + (a_l_m[ind_alpha] * a_l_m[ind_alpha - 1] * lam /
                                                           j_minus[ind_alpha - 2]) ** 2 * djdp_minus[ind_alpha - 2])

                w_prod_left = np.flip(np.cumprod(np.flip(lam * a_l_m[1:(ind_alpha + 1):2] * a_l_m[2:(ind_alpha + 2):2] /
                                                         j_minus[0:ind_alpha:2])))
                w_prod_right = np.cumprod(lam * a_l_m[(ind_alpha + 1):(l_max - abs(m)):2] *
                                          a_l_m[(ind_alpha + 2):(l_max - abs(m) + 1):2] /
                                          j_plus[(ind_alpha + 2):(l_max - abs(m) + 1):2])
                w_prod = np.concatenate((w_prod_left, np.ones(1), w_prod_right))
                residues[(abs(m) + 0)::2, (abs(m) + 0)::2, ind_alpha] = np.outer(w_prod, w_prod) * w_alpha
            elif np.mod(ind_alpha, 2) == 1:
                if ind_alpha == 1:
                    w_alpha = 1 / djdp_plus[1]
                else:
                    w_alpha = 1 / (djdp_plus[ind_alpha] + (a_l_m[ind_alpha] * a_l_m[ind_alpha - 1] * lam /
                                                           j_minus[ind_alpha - 2]) ** 2 * djdp_minus[ind_alpha - 2])

                w_prod_left = np.flip(np.cumprod(np.flip(lam * a_l_m[2:(ind_alpha + 1):2] * a_l_m[3:(ind_alpha + 2):2] /
                                                         j_minus[1:ind_alpha:2])))
                w_prod_right = np.cumprod(lam * a_l_m[(ind_alpha + 1):(l_max - abs(m)):2] *
                                          a_l_m[(ind_alpha + 2):(l_max - abs(m) + 1):2] /
                                          j_plus[(ind_alpha + 2):(l_max - abs(m) + 1):2])
                w_prod = np.concatenate((w_prod_left, np.ones(1), w_prod_right))
                residues[(abs(m) + 1)::2, (abs(m) + 1)::2, ind_alpha] = np.outer(w_prod, w_prod) * w_alpha

    return residues


def eval_a_l_m(l, m):
    r"""
    eval_a_l_m - Evaluate the coefficient from a ladder operation :math:`cos \theta Y_{\lambda;\mu}`
    on the spherical harmonic

    Parameters
    ----------
    l : int (array)
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    m : int
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    and Arfken (1999) (Ref [Arfken1999]_)
    """
    a_l_m = np.sqrt((l - m) * (l + m) /
                    ((2 * l - 1) * (2 * l + 1)))

    return a_l_m


def eval_b_l_m(l, m):
    r"""
    eval_a_l_m - Evaluate the coefficient from a ladder operation :math:`e^{i \phi} sin \theta Y_{\lambda;\mu}`
    on the spherical harmonic

    Parameters
    ----------
    l : int (array)
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    m : int
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    and Arfken (1999) (Ref [Arfken1999]_)
    """
    b_l_m = np.sqrt((l + m) * (l + m - 1) /
                    ((2 * l - 1) * (2 * l + 1)))

    return b_l_m


def calc_int_mag(length_kuhn, poles_vec):
    r"""
    Evaluate the magnitude of the integral for a list of poles (including repeats). This algorithm includes
    cases for single, double, and triple poles (as needed in evaluation of correlation functions)

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    poles_vec : float (array)
        Array of poles

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

        # Algorithm for the case of a simple pole
        if poles_order[i_pole] == 1:
            f_poles = 1 / (np.prod((poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) **
                                   np.delete(poles_order, i_pole)))
            int_mag += f_poles * np.exp(poles_vec_unique[i_pole] * length_kuhn)

        # Algorithm for the case of a double pole
        elif poles_order[i_pole] == 2:
            f_poles = 1 / (np.prod((poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)) **
                                   np.delete(poles_order, i_pole)))
            prod_poles = np.sum(np.delete(poles_order, i_pole) /
                                (poles_vec_unique[i_pole] - np.delete(poles_vec_unique, i_pole)))
            ddp_f_poles = - f_poles * prod_poles

            int_mag += (length_kuhn * f_poles + ddp_f_poles) * np.exp(poles_vec_unique[i_pole] * length_kuhn)

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

            int_mag += 0.5 * (length_kuhn ** 2 * f_poles + 2 * length_kuhn * ddp_f_poles +
                              d2dp2_f_poles) * np.exp(poles_vec_unique[i_pole] * length_kuhn)

    return int_mag
