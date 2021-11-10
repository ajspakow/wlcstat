import numpy as np


"""
Functions for evaluation of vertex functions
"""


def eval_residue_zero(k_val, mu = 0, lam_zero_only=True, lam_max=25, dimensions=3, lam_cont_frac_max=500):
    r"""
    eval_residues_zeros -
    Evaluate the residue at p = 0 using the intermediate-k algorithm provided in Ref. [Mehraeen2008]_

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    lam_cont_frac_max : int
        Maximum :math:`\lambda` value in the continued fraction evaluation

    Returns
    -------
    residue_zero, ddp_residue_zero : complex float
        Evaluated residues for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k algorithms

    """

    # Build the continued fractions
    j_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
    djdp_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
    a_lam_mu = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))

    lam = lam_cont_frac_max
    ind_lam = lam - abs(mu)
    j_plus[ind_lam] = lam * (lam + dimensions - 2)
    djdp_plus[ind_lam] = 1
    a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
    for lam in reversed(range(abs(mu), lam_cont_frac_max)):
        ind_lam = lam - abs(mu)
        if lam != abs(mu):
            a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
        j_plus[ind_lam] = (lam * (lam + dimensions - 2) + (a_lam_mu[ind_lam + 1] * k_val) ** 2 / j_plus[ind_lam + 1])
        djdp_plus[ind_lam] = 1 - (a_lam_mu[ind_lam + 1] * k_val / j_plus[ind_lam + 1]) ** 2 * djdp_plus[ind_lam + 1]

    if lam_zero_only:
        residue_zero = 1 / j_plus[0]
        ddp_residue_zero = - djdp_plus[0] / j_plus[0] ** 2.
    else:
        # Need to test this to verify unequal lambda values
        w_alpha = 1 / j_plus[0]
        ddp_w_alpha = - djdp_plus[0] / j_plus[0] ** 2.

        w_prod_right = np.cumprod(1j * k_val * a_lam_mu[1:(lam_max - abs(mu) + 1)] /
                                  j_plus[1:(lam_max - abs(mu) + 1)])
        w_prod = np.concatenate((np.ones(1), w_prod_right))

        ddp_w_prod_right = - np.cumsum(djdp_plus[1:(lam_max - abs(mu) + 1)] / j_plus[1:(lam_max - abs(mu) + 1)])
        ddp_w_prod = np.concatenate((np.zeros(1), ddp_w_prod_right))

        residue_zero = np.outer(w_prod, w_prod) * w_alpha
        ddp_residue_zero = (np.outer(w_prod, w_prod) * ddp_w_alpha
                            + np.outer(w_prod, w_prod * ddp_w_prod) * w_alpha
                            + np.outer(w_prod * ddp_w_prod, w_prod) * w_alpha)

    return residue_zero, ddp_residue_zero


def eval_residues_double_pole(k_val, mu, poles, lam_zero_only=True, lam_max=25, alpha_max=25,
                                     dimensions=3, lam_cont_frac_max=500):
    r"""
    eval_residues_double_pole -
    Evaluate the residues when a double pole occurs using the intermediate-k algorithm provided in Ref. [Mehraeen2008]_

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu`
    lam_zero_only : boolean
        Indicates whether the residues will be evaluated over the range of :math:`\lambda` and :math:`lambda_{0}`
    lam_max : int
        Maximum lambda value evaluated
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    lam_cont_frac_max : int
        Maximum :math:`\lambda` value in the continued fraction evaluation

    Returns
    -------
    residues : complex float
        Evaluated residues for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k algorithms

    """

    # Initialize the residue array based on whether lam_zero_only
    if lam_zero_only:
        residues_double = np.zeros((alpha_max - abs(mu) + 1), dtype=type(1+1j))
    else:
        residues_double = np.zeros((lam_max - abs(mu) + 1, lam_max - abs(mu) + 1,
                                    alpha_max - abs(mu) + 1), dtype=type(1+1j))

    for alpha in range(abs(mu), alpha_max + 1):
        ind_alpha = alpha - abs(mu)

        # Build the continued fractions
        j_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
        djdp_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
        d2jdp2_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1 + 1j))
        j_minus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
        djdp_minus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
        d2jdp2_minus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))
        a_lam_mu = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1+1j))

        lam = lam_cont_frac_max
        ind_lam = lam - abs(mu)
        j_plus[ind_lam] = poles[ind_alpha] + lam * (lam + dimensions - 2)
        djdp_plus[ind_lam] = 1
        d2jdp2_plus[ind_lam] = 0
        a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
        for lam in reversed(range(abs(mu), lam_cont_frac_max)):
            ind_lam = lam - abs(mu)
            if lam != abs(mu):
                a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
            j_plus[ind_lam] = (poles[ind_alpha] + lam * (lam + dimensions - 2)
                               + (a_lam_mu[ind_lam + 1] * k_val) ** 2 / j_plus[ind_lam + 1])
            djdp_plus[ind_lam] = 1 - ((a_lam_mu[ind_lam + 1] * k_val) ** 2 /
                                      j_plus[ind_lam + 1] ** 2 * djdp_plus[ind_lam + 1])
            d2jdp2_plus[ind_lam] = (2 * (a_lam_mu[ind_lam + 1] * k_val) ** 2 /
                                    j_plus[ind_lam + 1] ** 3 * djdp_plus[ind_lam + 1] ** 2
                                    - (a_lam_mu[ind_lam + 1] * k_val) ** 2 /
                                    j_plus[ind_lam + 1] ** 2 * d2jdp2_plus[ind_lam + 1])
        lam = abs(mu)
        j_minus[0] = poles[ind_alpha] + lam * (lam + dimensions - 2)
        djdp_minus[0] = 1
        d2jdp2_minus[0] = 0
        for lam in range(abs(mu) + 1, max(lam_max, alpha_max) + 1):
            ind_lam = lam - abs(mu)
            j_minus[ind_lam] = (poles[ind_alpha] + lam * (lam + dimensions - 2)
                                + (a_lam_mu[ind_lam] * k_val) ** 2 / j_minus[ind_lam - 1])
            djdp_minus[ind_lam] = 1 - ((a_lam_mu[ind_lam] * k_val) ** 2 /
                                       j_minus[ind_lam - 1] ** 2 * djdp_minus[ind_lam - 1])
            d2jdp2_minus[ind_lam] = (2 * (a_lam_mu[ind_lam] * k_val) ** 2 /
                                     j_minus[ind_lam - 1] ** 3 * djdp_minus[ind_lam - 1] ** 2
                                     - (a_lam_mu[ind_lam] * k_val) ** 2 /
                                     j_minus[ind_lam - 1] ** 2 * d2jdp2_minus[ind_lam - 1])

        if lam_zero_only:
            if ind_alpha == 0:
                w_alpha = 1 / djdp_plus[0]
                d2jdp2_alpha = d2jdp2_plus[0]
                residues_double[ind_alpha] = - d2jdp2_alpha * w_alpha ** 3
            else:
                w_alpha = 1 / (djdp_plus[ind_alpha] - (a_lam_mu[ind_alpha] * k_val /
                                                       j_minus[ind_alpha - 1]) ** 2 * djdp_minus[ind_alpha - 1])
                w_prod_left = np.prod(1j * k_val * a_lam_mu[1:(ind_alpha + 1)] / j_minus[0:ind_alpha])
                d2jdp2_alpha = d2jdp2_plus[ind_alpha] + (2 * (a_lam_mu[ind_alpha] * k_val) ** 2 /
                                     j_minus[ind_alpha - 1] ** 3 * djdp_minus[ind_alpha - 1] ** 2
                                     - (a_lam_mu[ind_alpha] * k_val) ** 2 /
                                     j_minus[ind_alpha - 1] ** 2 * d2jdp2_minus[ind_alpha - 1])

                residues_double[ind_alpha] = - w_prod_left ** 2 * d2jdp2_alpha * w_alpha ** 3
        else:
            if ind_alpha == 0:
                w_alpha = 1 / djdp_plus[0]
                d2jdp2_alpha = d2jdp2_plus[0]
            else:
                w_alpha = 1 / (djdp_plus[ind_alpha] -
                               (a_lam_mu[ind_alpha] * k_val / j_minus[ind_alpha - 1]) ** 2 * djdp_minus[ind_alpha - 1])
                d2jdp2_alpha = d2jdp2_plus[ind_alpha] + (
                        2 * (a_lam_mu[ind_alpha] * k_val) ** 2 / j_minus[ind_alpha - 1] ** 3 *
                        djdp_minus[ind_alpha - 1] ** 2
                        - (a_lam_mu[ind_alpha] * k_val) ** 2 / j_minus[ind_alpha - 1] ** 2 *
                        d2jdp2_minus[ind_alpha - 1])

            w_prod_left = np.flip(np.cumprod(np.flip(1j * k_val * a_lam_mu[1:(ind_alpha + 1)] / j_minus[0:ind_alpha])))
            w_prod_right = np.cumprod(1j * k_val * a_lam_mu[(ind_alpha + 1):(lam_max - abs(mu) + 1)] /
                                      j_plus[(ind_alpha + 1):(lam_max - abs(mu) + 1)])
            w_prod = np.concatenate((w_prod_left, np.ones(1), w_prod_right))

            ddp_w_prod_right = - np.cumsum(
                djdp_plus[(ind_alpha + 1):(lam_max - abs(mu) + 1)] / j_plus[(ind_alpha + 1):(lam_max - abs(mu) + 1)])
            ddp_w_prod_left = - np.flip(np.cumsum(np.flip(djdp_minus[0:ind_alpha] / j_minus[0:ind_alpha])))
            ddp_w_prod = np.concatenate((ddp_w_prod_left, np.zeros(1), ddp_w_prod_right))

            residues_double[:, :, ind_alpha] = (-np.outer(w_prod ** 2, w_prod ** 2) * d2jdp2_alpha * w_alpha ** 3
                                                + 2 * np.outer(w_prod ** 2, w_prod ** 2 * ddp_w_prod) * w_alpha ** 2
                                                + 2 * np.outer(w_prod ** 2 * ddp_w_prod, w_prod ** 2) * w_alpha ** 2)

    return residues_double


def eval_residues_other_pole(k_val, mu, poles, lam_zero_only=True, lam_max=25, alpha_max=25,
                              dimensions=3, lam_cont_frac_max=500):
    r"""
    eval_residues_double_pole -
    Evaluate the residues at a different pole occurs using the intermediate-k algorithm provided in Ref. [Mehraeen2008]_

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu`
    lam_zero_only : boolean
        Indicates whether the residues will be evaluated over the range of :math:`\lambda` and :math:`lambda_{0}`
    lam_max : int
        Maximum lambda value evaluated
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    lam_cont_frac_max : int
        Maximum :math:`\lambda` value in the continued fraction evaluation

    Returns
    -------
    residues : complex float
        Evaluated residues for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k algorithms

    """

    # Initialize the residue array based on whether lam_zero_only
    if lam_zero_only:
        residues_other_pole = np.zeros((alpha_max - abs(mu) + 1), dtype=type(1 + 1j))
    else:
        residues_other_pole = np.zeros((lam_max - abs(mu) + 1, lam_max - abs(mu) + 1,
                                    alpha_max - abs(mu) + 1), dtype=type(1 + 1j))

    for alpha in range(abs(mu), alpha_max + 1):
        ind_alpha = alpha - abs(mu)

        # Build the continued fractions
        j_plus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1 + 1j))
        j_minus = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1 + 1j))
        a_lam_mu = np.zeros((lam_cont_frac_max - abs(mu) + 1), dtype=type(1 + 1j))

        lam = lam_cont_frac_max
        ind_lam = lam - abs(mu)
        j_plus[ind_lam] = poles[ind_alpha] + lam * (lam + dimensions - 2)
        a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
        for lam in reversed(range(abs(mu), lam_cont_frac_max)):
            ind_lam = lam - abs(mu)
            if lam != abs(mu):
                a_lam_mu[ind_lam] = eval_a_lam_mu(lam, mu, dimensions)
            j_plus[ind_lam] = (poles[ind_alpha] + lam * (lam + dimensions - 2)
                               + (a_lam_mu[ind_lam + 1] * k_val) ** 2 / j_plus[ind_lam + 1])
        lam = abs(mu)
        j_minus[0] = poles[ind_alpha] + lam * (lam + dimensions - 2)
        for lam in range(abs(mu) + 1, max(lam_max, alpha_max) + 1):
            ind_lam = lam - abs(mu)
            if j_minus[ind_lam - 1] != 0:
                j_minus[ind_lam] = (poles[ind_alpha] + lam * (lam + dimensions - 2)
                                    + (a_lam_mu[ind_lam] * k_val) ** 2 / j_minus[ind_lam - 1])

        if lam_zero_only:
            if ind_alpha == 0:
                w_alpha = 1 / j_plus[0]
                residues_other_pole[ind_alpha] = w_alpha
            else:
                w_alpha = 1 / (j_plus[ind_alpha] + (a_lam_mu[ind_alpha] * k_val) ** 2 / j_minus[ind_alpha - 1])

                w_prod_left = np.prod(1j * k_val * a_lam_mu[1:(ind_alpha + 1)] / j_minus[0:ind_alpha])
                residues_other_pole[ind_alpha] = w_prod_left ** 2 * w_alpha
        else:
            if ind_alpha == 0:
                w_alpha = 1 / j_plus[0]
            else:
                w_alpha = 1 / (j_plus[ind_alpha] + (a_lam_mu[ind_alpha] * k_val) ** 2 / j_minus[ind_alpha - 1])

            w_prod_left = np.flip(np.cumprod(np.flip(1j * k_val * a_lam_mu[1:(ind_alpha + 1)] / j_minus[0:ind_alpha])))
            w_prod_right = np.cumprod(1j * k_val * a_lam_mu[(ind_alpha + 1):(lam_max - abs(mu) + 1)] /
            j_plus[(ind_alpha + 1):(lam_max - abs(mu) + 1)])
            w_prod = np.concatenate((w_prod_left, np.ones(1), w_prod_right))
            residues_other_pole[:, :, ind_alpha] = np.outer(w_prod, w_prod) * w_alpha

    return residues_other_pole


def eval_a_lam_mu(lam, mu, dimensions=3):
    r"""
    eval_a_lam_mu - Evaluate the coefficient from a ladder operation :math:`cos \theta Y_{\lambda;\mu}`
    on the spherical harmonic

    Parameters
    ----------
    lam : int (array)
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    mu : int
        The angular kinetic energy quantum index of the spherical harmonic :math:`Y_{\lambda;\mu}`

    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    and Arfken (1999) (Ref [Arfken1999]_)
    """
    a_lam_mu = np.sqrt((lam - mu) * (lam + mu + dimensions - 3) /
                       ((2 * lam + dimensions - 2) * (2 * lam + dimensions - 4)))

    return a_lam_mu
