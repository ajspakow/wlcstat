import numpy as np


"""
Functions for evaluation of vertex functions
"""


def eval_residue_zero(k_val, dimensions=3, lam_cont_frac_max=500):
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
    j_plus = np.zeros((lam_cont_frac_max + 1), dtype=type(1+1j))
    djdp_plus = np.zeros((lam_cont_frac_max + 1), dtype=type(1+1j))
    a_lam_mu = np.zeros((lam_cont_frac_max + 1), dtype=type(1+1j))

    lam = lam_cont_frac_max
    ind_lam = lam
    j_plus[ind_lam] = lam * (lam + dimensions - 2)
    djdp_plus[ind_lam] = 1
    a_lam_mu[ind_lam] = eval_a_lam_mu(lam, 0, dimensions)
    for lam in reversed(range(0, lam_cont_frac_max)):
        ind_lam = lam
        if lam != 0:
            a_lam_mu[ind_lam] = eval_a_lam_mu(lam, 0, dimensions)
        j_plus[ind_lam] = (lam * (lam + dimensions - 2) + (a_lam_mu[ind_lam + 1] * k_val) ** 2 / j_plus[ind_lam + 1])
        djdp_plus[ind_lam] = 1 - (a_lam_mu[ind_lam + 1] * k_val / j_plus[ind_lam + 1]) ** 2 * djdp_plus[ind_lam + 1]

    residue_zero = 1 / j_plus[0]
    ddp_residue_zero = - djdp_plus[0] / j_plus[0] ** 2.

    return residue_zero, ddp_residue_zero


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
