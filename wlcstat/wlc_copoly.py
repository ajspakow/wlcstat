
from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np
import scipy as sp


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