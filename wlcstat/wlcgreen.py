import numpy as np
from wlcstat.util.wlc_poles_residues import *


def eval_poles_and_residues(k_val, mu, dimensions=3):
    r"""
    eval_poles_and_residues - Evaluate the poles and the residues for a given value of the
    Fourier vector magnitude :math:`K`

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    poles : float (complex)

    residues : float (complex)

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    residues = np.zeros(1)
    poles = eval_poles(k_val, mu, dimensions)

    return poles, residues


r"""
    mu_zero_only : boolean
        Determines whether poles and residues are determined for non-zero :math:`\mu` (default True)
    lam_zero_only : boolean
        Determines whether residues are determined for non_zero :math:`\lambda` (default True)
"""
