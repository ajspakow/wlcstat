from wlcstat.util.wlc_poles_residues import *
from wlcstat.util.wlc_vertex import *
from wlcstat.wlcave import *
import numpy as np
import scipy as sp


# Component of structure factor between two monomers

def s2_wlc_monomers(k_val_vector, delta, length_kuhn, epsilon=1, dimensions=3, alpha_max=25):
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
        
    if type(k_val_vector) == float or type(k_val_vector) == int:
        k_val_vector = np.array([k_val_vector])

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
                if k_val < k_cutoff:
                    s2[ind_k_val, ind_delta] = 1
            elif delta[ind_delta] >= epsilon:
                for alpha in range(0, alpha_max):
                    s2[ind_k_val, ind_delta] += residues[alpha] * (np.exp(poles[alpha]*(delta[ind_delta]+epsilon)) + np.exp(poles[alpha]*(delta[ind_delta]-epsilon)) - 2*np.exp(poles[alpha]*delta[ind_delta])) / poles[alpha] ** 2
                
                s2[ind_k_val, ind_delta] *= 1 / epsilon ** 2
                
                # Reset the s2 value if below the cutoff
                if k_val < k_cutoff:
                    s2[ind_k_val, ind_delta] = 1
    return s2


# Two-point structure factors weighted by protein binding to marks

def s2_wlc_marked(k_val_vector, N, M, exp_sigma, exp_sigma_squared, epsilon=1, dimensions=3, alpha_max=25, delta_max = np.Inf):
    r"""
    s2_wlc_marked - Evaluate the 2-point structure factor for the wormlike chain model weighted by protein binding for various types

    Parameters
    ----------
    k_val_vector : float (array)
        The value of the Fourier vector magnitude :math:`K`
    N : int
        The number of monomers in the polymer.
    M : int
        The number of mark types.
    exp_sigma : float (N-by-M array)
        The expected number of proteins of each type (M) bound to each monomer (N) of the polymer.
    exp_sigma_squared : float (N-by-M array)
        The expected squared number <sigma^2> of proteins of each type (M) bound to each monomer (N) of the polymer.
    epsilon : float
        The monomer length in Kuhn lengths (default 1)
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    delta_max : int
        Maximum interaction distance along polymer chain that contributes to structure factor (default Infinity)

    Returns
    -------
    s2List : List of float (vector)
        First element: sum over all monomer pairs of S_{ij}<sigma_j> for each value of k and mark type
        Second element: sum over all monomer pairs of S_{ij}<sigma_i^\alpha><sigma_j^\beta> for each value of k and pairs of mark types (1,1; 1,2; 1,3; ... 1,M; 2,2; ... M,M)
    """
    
    if type(k_val_vector) == float or type(k_val_vector) == int:
        k_val_vector = np.array([k_val_vector])
    
    if type(exp_sigma) == float or type(exp_sigma) == int:
        exp_sigma = exp_sigma*np.ones((N,M))
        
    if type(exp_sigma_squared) == float or type(exp_sigma_squared) == int:
        exp_sigma_squared = exp_sigma_squared*np.ones((N,M))
        
    deltas = np.arange(0, np.min([N, delta_max]))
    s_monos = s2_wlc_monomers(k_val_vector, deltas*epsilon, N*epsilon, epsilon, dimensions, alpha_max)
    
    one_mark_coeffs = np.zeros((len(deltas),M))
    two_mark_coeffs = np.zeros((len(deltas),M,M))
    
    for ind_mark1 in range(0,M):
        one_mark_coeffs[0,ind_mark1] = 2*np.sum(exp_sigma[:,ind_mark1])
        for ind_mark2 in range(ind_mark1,M):
            if ind_mark1 == ind_mark2:
                two_mark_coeffs[0,ind_mark1,ind_mark2] = np.sum(exp_sigma_squared[:,ind_mark1])
            else:
                two_mark_coeffs[0,ind_mark1,ind_mark2] = np.dot(exp_sigma[:,ind_mark1],exp_sigma[:,ind_mark2])
    
    for ind_delta in range(1,len(deltas)):
        for ind_mark1 in range(0,M):
            one_mark_coeffs[ind_delta,ind_mark1] = one_mark_coeffs[ind_delta-1,ind_mark1] - (exp_sigma[ind_delta-1,ind_mark1] + exp_sigma[-ind_delta,ind_mark1])
            for ind_mark2 in range(ind_mark1,M):
                two_mark_coeffs[ind_delta,ind_mark1,ind_mark2] = np.dot(exp_sigma[:-ind_delta,ind_mark1], exp_sigma[ind_delta:,ind_mark2]) + np.dot(exp_sigma[ind_delta:,ind_mark1], exp_sigma[:-ind_delta,ind_mark2])
    
    one_mark_coeffs[0,:] /= 2
        
    s2_one_mark = np.zeros((len(k_val_vector), M), dtype=type(1+1j))
    s2_two_marks = np.zeros((len(k_val_vector), int(np.round(M*(M+1)/2))), dtype=type(1+1j))
    
    ind_marks = 0
    for ind_mark1 in range(0,M):
        s2_one_mark[:,ind_mark1] = np.tensordot(s_monos,one_mark_coeffs[:,ind_mark1],1)
        for ind_mark2 in range(ind_mark1,M):
            s2_two_marks[:,ind_marks] = np.tensordot(s_monos,two_mark_coeffs[:,ind_mark1,ind_mark2],1)
            ind_marks += 1
 
    s2_one_mark /= N ** 2
    s2_two_marks /= N ** 2
    
    return [s2_one_mark, s2_two_marks]


# Calculate expected number of reader proteins bound to each nucleosome

def exp_sigma_calc(marks, mu, e_bind, J, phi_p, v_p=np.pi/6):
    r"""
    exp_sigma_calc - Evaluate <sigma> and <sigma^2> as a function of marks

    Parameters
    ----------
    marks : N-by-M int array
        The number of tails marked with each mark (M) on each nucleosome (N)
    mu : M-by-1 float array
        Chemical potential of each reader protein (units k_BT)
    e_bind : M-by-1 float array
        The binding energy for each mark (units k_BT)
    J : M-by-M float array
        Symmetric array with coupling energy per volume between reader proteins
    phi_p : Float, 0 < phi_p < 1
        Volume fraction taken up by polymer vs. solvent
    v_p : Float
        Volume of a nucleosome in Kuhn lengths^3 (default = pi/6)

    Returns
    -------
    sigmaList : list of two N-by-M float arrays
        <sigma> and <sigma^2> for each mark type on each nucleosome
    """
    
    N = marks.shape[0]
    
    if len(marks.shape) == 1:
        marks = marks.reshape((N,1))
        
    M = marks.shape[1]
    
    # Numerically solve mean protein densities rho
    def density_func(rho):
        sm = np.zeros(M)
        for i in range(0,N):            
            kappa = mu+np.multiply(e_bind,marks[i,:])+np.sum(np.multiply(J,rho),0)
            sm += np.divide((2*np.exp(kappa) + 2*np.exp(2*kappa)), (1+2*np.exp(kappa) + np.exp(2*kappa)))/N
        return sm - rho*v_p/phi_p

    # Avoid wells by starting search at both maximum and minimum rho
    # One of these will often trigger a warning, but the other should be fine
    rho_hi = sp.optimize.fsolve(density_func,np.ones(M)*2*phi_p/v_p)
    rho_lo = sp.optimize.fsolve(density_func,np.zeros(M))
    
    rho = rho_lo
    if np.linalg.norm(density_func(rho_hi)) + 1e-10 < np.linalg.norm(density_func(rho_lo)):
        rho = rho_hi
    print("Error = " + str(np.linalg.norm(density_func(rho))))
        
    sigma = np.zeros((N,M))
    sigma_squared = np.zeros((N,M))
    
    for i in range(0,N):
        kappa = mu+np.multiply(e_bind,marks[i,:])+np.sum(np.multiply(J,rho),0)
        sigma[i,:] = np.divide((2*np.exp(kappa) + 2*np.exp(2*kappa)), (1+2*np.exp(kappa) + np.exp(2*kappa)))
        sigma_squared[i,:] = np.divide((2*np.exp(kappa) + 4*np.exp(2*kappa)), (1+2*np.exp(kappa) + np.exp(2*kappa)))
        
    return [sigma,sigma_squared]
