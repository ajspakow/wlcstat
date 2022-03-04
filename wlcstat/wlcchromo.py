from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np
import scipy as sp


# Component of structure factor between two monomers

def s2_wlc_monomers(k_val_vector, delta, epsilon=1, length_kuhn, dimensions=3, alpha_max=25):
    r"""
    s2_wlc_monomers - Evaluate the component of the 2-point structure factor between two monomers for the wormlike chain model

    Parameters
    ----------
    k_val_vector : float (array)
        The value of the Fourier vector magnitude :math:`K`
    delta : float (array)
        The arc length separation between the two monomers in Kuhn lengths. Note that either delta >= epsilon or delta = 0.
    epsilon : float
        The monomer length in Kuhn lengths (default 1)
    length_kuhn : float
        The polymer length in Kuhn lengths, used to determine cutoff k value
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    s2 : float (vector)
        Component of the structure factor from these two monomers for the wormlike chain model for every k_val in k_val_vector
    """
    
    if type(delta) == float or type(delta) == int:
        delta = np.array([delta])

    s2 = np.zeros((len(k_val_vector), len(delta)), dtype=type(1+1j))

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

        for ind_delta in range(0, len(delta)):
            if delta[ind_delta] == 0:
                s2[ind_k_val, ind_delta] = (epsilon * residue_zero + ddp_residue_zero)

                for alpha in range(0, alpha_max):
                    s2[ind_k_val, ind_delta] += (np.exp(poles[alpha] * epsilon) *
                                                      residues[alpha] / poles[alpha] ** 2)

                s2[ind_k_val, ind_delta] *= 2 / epsilon ** 2

                # Reset the s2 value if below the cutoff
                if k_val < k_cutoff[ind_length]:
                    s2[ind_k_val, ind_length] = 1
            elif delta[ind_delta] >= epsilon:
                for alpha in range(0, alpha_max):
                    s2[ind_k_val, ind_delta] += residues[alpha]*np.exp(poles[alpha] * delta[ind_delta]) * (np.cosh(poles[alpha] * delta[ind_delta]) - 1) / poles[alpha] ** 2
                
                s2[ind_k_val, ind_delta] *= 2 / epsilon ** 2
                
                # Reset the s2 value if below the cutoff
                if k_val < k_cutoff[ind_length]:
                

    return s2
