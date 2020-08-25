
from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
import numpy as np


def eval_structure_factor(k_val_vector, length_kuhn, dimensions=3, alpha_max=25):
    r"""
    eval_structure factor - Evaluate the structure factor for the wormlike chain model

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
    structure_factor : float (vector)
        Structure factor for the wormlike chain model for every k_val in k_val_vector

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """
    if type(length_kuhn) == float or type(length_kuhn) == int:
        structure_factor = np.zeros((len(k_val_vector)), dtype=type(1+1j))
    else:
        structure_factor = np.zeros((len(k_val_vector), len(length_kuhn)), dtype=type(1+1j))

    for ind_k_val in range(0, len(k_val_vector)):
        k_val = k_val_vector[ind_k_val]

        poles = eval_poles(k_val, 0, dimensions, alpha_max)
        residues = eval_residues(k_val, 0, poles, True, dimensions, alpha_max)
        residue_zero, ddp_residue_zero = eval_residue_zero(k_val, dimensions)

        if type(length_kuhn) == float or type(length_kuhn) == int:
            structure_factor[ind_k_val] = (length_kuhn * residue_zero + ddp_residue_zero)
            for alpha in range(0, alpha_max):
                structure_factor[ind_k_val] += (np.exp(poles[alpha] * length_kuhn) *
                                                residues[alpha] / poles[alpha] ** 2)
            structure_factor[ind_k_val] *= 2 / length_kuhn ** 2
        else:
            for ind_length in range(0, len(length_kuhn)):
                structure_factor[ind_k_val, ind_length] = (length_kuhn[ind_length] * residue_zero + ddp_residue_zero)

                for alpha in range(0, alpha_max):
                    structure_factor[ind_k_val, ind_length] += (np.exp(poles[alpha] * length_kuhn[ind_length]) *
                                                                residues[alpha] / poles[alpha] ** 2)

                structure_factor[ind_k_val, ind_length] *= 2 / length_kuhn[ind_length] ** 2

    return structure_factor


