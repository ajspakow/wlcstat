
from wlcstat.util.wlc_poles_residues import *
import numpy as np
import scipy.special as sp


def eval_poles_and_residues(k_val, mu, lam_zero_only=True, dimensions=3):
    r"""
    eval_poles_and_residues - Evaluate the poles and the residues for a given value of the
    Fourier vector magnitude :math:`K`

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter (:math:`z`-component of the angular momentum)
    lam_zero_only : boolean
        Determines whether residues are determined for non_zero :math:`\lambda` (default True)
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu`
    residues : complex float
        Evaluated residues for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    poles = eval_poles(k_val, mu, dimensions)
    residues = eval_residues(k_val, mu, poles, lam_zero_only, dimensions)

    return poles, residues


def gwlc_r(r_val, length_kuhn, dimensions=3, alpha_max=25, k_val_max=1e5, delta_k_val_max=0.1):
    r"""
    gwlc_r - Evaluate the orientation-independent Green's function for the wormlike chain model
    
    Parameters
    ----------
    r_val : float (array)
        The values of the end-to-end distance :math:`r = R/L` to be evaluated
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    k_val_max : float
        Cutoff value of :math:`K` for numerical integration
    delta_k_val_max : float
        Maximum value of the integration step size

    Returns
    -------
    gwlc : float
        The orientation-independent Green's function [size len(r_val) x len(length_kuhn)]

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    delta_k_val = min(delta_k_val_max, 2 * np.pi / np.max(length_kuhn) / 10)

    # Eliminate 0 and 1 from the r_val array
    r_val[r_val == 0] = 1e-10
    r_val[r_val == 1] = 1-1e-10

    # Initialize the Green's function
    if type(length_kuhn) == float or type(length_kuhn) == int:
        gwlc = np.zeros((len(r_val)), dtype=type(1+1j))
    else:
        gwlc = np.zeros((len(r_val), len(length_kuhn)), dtype=type(1+1j))

    tolerance = 1e-15

    k_val_output = k_val_max / 100
    k_val = delta_k_val
    int_count = 0
    contains_nan = False
    while k_val <= k_val_max and not contains_nan:
        int_count += 1

        poles = eval_poles(k_val, 0, dimensions, alpha_max)
        residues = eval_residues(k_val, 0, poles, True, dimensions, alpha_max, alpha_max)

        if int_count == 1:
            int_coef = 55 / 24
        elif int_count == 2:
            int_coef = -1 / 6
        elif int_count == 3:
            int_coef = 11 / 8
        else:
            int_coef = 1

        for alpha in range(0, alpha_max + 1):
            gkwlc_kval = residues[alpha] * np.exp(poles[alpha] * length_kuhn)
            if type(length_kuhn) == float or type(length_kuhn) == int:
                integrand = (int_coef * k_val ** (dimensions / 2)
                             * sp.jv(dimensions / 2 - 1, k_val * r_val * length_kuhn) * gkwlc_kval)
            else:
                integrand = (int_coef * k_val ** (dimensions / 2)
                             * sp.jv(dimensions / 2 - 1, k_val * np.outer(r_val, length_kuhn))
                             * np.outer(np.ones((len(r_val)), dtype=type(1+1j)), gkwlc_kval))
            if not contains_nan:
                contains_nan = np.isnan(integrand).any()
            if not contains_nan:
                gwlc += integrand
            else:
                print("Encountered NaN at k_val = " + str(k_val))

        k_val += delta_k_val
        if k_val >= k_val_output:
            k_val_output += k_val_max / 100
            print("Current k_val = " + str(k_val) + " with delta_k_val = " + str(delta_k_val)
                  + " and k_val_max = " + str(k_val_max))
            print(np.max(abs(np.exp(poles * np.min(length_kuhn)))))

        if np.min(abs(np.exp(poles * np.min(length_kuhn)))) < tolerance and alpha_max > 10 and k_val > 2840:
            alpha_max -= 1
            print("Reducing alpha_max to " + str(alpha_max) + " at k_val = "
                  + str(k_val) + " with delta_k_val = " + str(delta_k_val))
            print(np.max(abs(np.exp(poles * np.min(length_kuhn)))))

        if np.max(abs(np.exp(poles * np.min(length_kuhn)))) < tolerance:
            print("Achieved accuracy of " + str(tolerance) + " at k_val = "
                  + str(k_val) + " with delta_k_val = " + str(delta_k_val))
            print(np.max(abs(np.exp(poles * np.min(length_kuhn)))))
            k_val = k_val_max + delta_k_val

    gwlc *= delta_k_val / (2 * np.pi) ** (dimensions / 2) * np.outer(r_val ** (-dimensions / 2 + 1),
                                                                     length_kuhn ** (dimensions / 2 + 1))

    return gwlc
