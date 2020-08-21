import numpy as np


def eval_poles(k_val, mu, dimensions=3, alpha_max=25, k_val_cutoff=8000):
    r"""
    eval_poles - Evaluate the poles for the wormlike chain Green's function for a given :math:`K`
    and :math:`z`-component quantum index :math:`\mu`

    Parameters
    ----------

    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)
    k_val_cutoff : float
        Cutoff value of :math:`K` for crossover from intermediate-k algorithm to large-k algorithm

    Returns
    -------
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k and large-k algorithms

    """

    # Size of 'poles' set by the pole index running from :math:`\mu` to :math:`\alpha_{max}`
    poles = np.zeros((alpha_max + 1 - abs(mu)), dtype='complex')

    # Use the intermediate-k or large-k algorithm based on k_val_cutoff
    if abs(k_val) < k_val_cutoff:
        poles = eval_poles_intermediate_k_val(k_val, mu, dimensions, alpha_max)
    else:
        poles = eval_poles_large_k_val(k_val, mu, dimensions, alpha_max)

    return poles


def eval_poles_intermediate_k_val(k_val, mu, dimensions=3, alpha_max=25):
    r"""
    eval_poles_intermediate_k_val - Evaluate the poles for given :math:`K` and :math:`\mu`
    using the matrix method for intermediate :math:`K`

    Parameters
    ----------
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu` using the intermediate-k algorithm

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k algorithms

    """

    # Determine the number of poles evaluated noting conjugate pairing
    # Evaluate 4 times the number needed to achieve accuracy

    num_total = int(4 * 2 * (np.ceil(alpha_max / 2.0)))  # Eigenvalues come in pairs

    h_matrix = eval_h_matrix(num_total, k_val, mu, dimensions)

    poles_total = -1 * np.linalg.eigvals(h_matrix)
    poles_total = np.sort(poles_total)[::-1]

    # Fix sorting
    for ii in range(0, (num_total - 1), 2):
        poles_total[ii] = np.real(poles_total[ii]) + 1j * abs(np.imag(poles_total[ii]))
        poles_total[ii + 1] = np.real(poles_total[ii + 1]) - 1j * abs(np.imag(poles_total[ii + 1]))

    if k_val > 1:
        poles_total = poles_total * k_val

    poles = poles_total[0:(alpha_max - abs(mu)+1)]

    return poles


def eval_h_matrix(num_total, k_val, mu, dimensions=3):
    r"""
    eval_h_matrix - Build the h-matrix used to evaluate the poles of the Green function

    Parameters
    ----------
    num_total : int
        The rank of the h-matrix that determines the number of poles to evaluate by eigenvalue decomposition
    mu : int
        Value of the mu parameter
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    h_matrix : complex float (num_total x num_total)
        h_matrix for evaluation of the poles for the given :math:`K` and :math:`\mu`

    Notes
    -----
    See [Mehraeen2009]_ for intermediate-k algorithms

    """
    h_matrix = np.zeros((num_total, num_total), dtype=type(1 + 1j))
    left = (0 + 0j) * np.NaN

    for row in range(0, num_total):
        j = row + 1
        if k_val <= 1:
            # for above diagonal
            right = -1j * k_val * eval_a_lam_mu(j + mu, mu, dimensions)
            if row > 0:
                # for below diagonal
                left = -1j * k_val * eval_a_lam_mu(j + mu - 1, mu, dimensions)
            # diagonal element
            diag = (mu + j - 1) * (mu + j + dimensions - 3)
            if row == 0:
                h_matrix[row, 0:2] = [diag, right]
            elif row == num_total - 1:
                h_matrix[row, row-1:row+1] = [left, diag]
            else:
                h_matrix[row, row-1:row+2] = [left, diag, right]
        else:
            # for above diagonal
            right = -1j * eval_a_lam_mu(j + mu, mu, dimensions)
            if row > 0:
                # for below diagonal
                left = -1j * eval_a_lam_mu(j + mu - 1, mu, dimensions)
            # diagonal element
            diag = (mu + j - 1) * (mu + j + dimensions - 3.0) / k_val
            if row == 0:
                h_matrix[row, 0:2] = [diag, right]
            elif row == num_total - 1:
                h_matrix[row, row-1:row+1] = [left, diag]
            else:
                h_matrix[row, row-1:row+2] = [left, diag, right]

    return h_matrix


def eval_poles_large_k_val(k_val, mu, dimensions=3, alpha_max=25):
    r"""
    eval_poles_large_k_val - Evaluate the poles for given :math:`K` and :math:`\mu`
    using the expansion method (perturbation theory) for large :math:`K`

    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter
    dimensions : int
        The number of dimensions (default to 3 dimensions)
    alpha_max : int
        Maximum number of poles evaluated (default 25)

    Returns
    -------
    poles : complex float
        Evaluated poles for the given :math:`K` and :math:`\mu` using the large-k algorithm

    Notes
    -----
    See [Mehraeen2009]_ for large-k algorithms

    """
    # k, ORDEig, mu, d=3):

    poles_total = np.zeros(alpha_max + 1, dtype='complex') * np.NaN
    for l_val in range(0, alpha_max, 2):
        poles_total[l_val] = 1j * k_val - mu * (mu + dimensions - 2) - eval_epsilon(l_val / 2.0, dimensions, k_val, mu)
        poles_total[l_val + 1] = np.conj(poles_total[l_val])

    poles = poles_total[0:(alpha_max - abs(mu)+1)]

    return poles


def eval_epsilon(l_val, d, k_val, mu):
    r"""
    eval_epsilon - Evaluation of the large-k expansion coefficients for asymptotic analysis of poles

    Parameters
    ----------
    l_val : float
    d : int
        Dimensions
    k_val : float
        The value of the Fourier vector magnitude :math:`K`
    mu : int
        Value of the mu parameter

    Returns
    -------
    epsilon_sum : complex float
        The value of the large-k expansion from the quantum mechanical perturbation theory

    Notes
    -----
    See [Mehraeen2009]_ for large-k algorithms

    """
    alpha = 1.0 / np.sqrt(8 * k_val)
    beta = -np.sqrt(2) / 4 * (1+1j)
    m = mu + (d - 3) / 2.0
    n = 2 * l_val + m + 1

    epsilon_0 = (-1.0/2/beta)**(-1)*(n/2.0)
    epsilon_1 = (-1.0/2/beta)**(0)*(-1.0/8*(n**2+3-3*m**2)-m*(m+1))
    epsilon_2 = (-1.0/2/beta)**(1)*(-1.0/2**5*n*(n**2+3-9*m**2))
    epsilon_3 = (-1.0/2/beta)**(2)*(-1.0/2**8*(5*n**4+34*n**2+9)-(102*n**2+42)*m**2+33*m**4)
    epsilon_4 = (-1.0/2/beta)**(3)*(-1.0/2**11*n*(33*n**4+410*n**2+405)-(1230*n**2+1722)*m**2+813*m**4)
    epsilon_5 = (-1.0/2/beta)**(4)*(-1.0/2**12*9*(7*n**6+140*n**4+327*n**2+54
                                    - (420*n**4+1350*n**2+286)*m**2+(495*n**2+314)*m**4-82*m**6))
    epsilon_sum = epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha**2+epsilon_4*alpha**3+epsilon_5*alpha**4

    return epsilon_sum



"""

Functions for residue evaluation

"""


def eval_residues(k_val, eig, lam, mu, nlam=10, dimensions=3, lam_max=500, cutoff=10**-11):
    r"""

    if k_val < 10 ** -3:
        res = np.zeros((nlam, nlam), dtype=type(1+1j))*np.NaN
        for lam in range(abs(mu),nlam):
            for lam0 in range(abs(mu),nlam):
                smallK = SmallAysmpRes(k_val, l, lam, lam0, mu, dimensions)
                res[lam0, lam] = smallK

    res = largeKResidues(K, eig, mu, nlam, lam_max)

    for lam in range(abs(mu),nlam):
        for lam0 in range(abs(mu),nlam):
            smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
            if abs(smallK)<cutoff:
                res[lam0,lam]=smallK
    return res
    """

    pass


def eval_a_lam_mu(lam, mu, dimensions=3):
    r"""
    eval_a_lam_mu - Evaluate the coefficient from a ladder operation :math:`cos \theta Y_{\lambda;\mu}`
    on the spherical harmonic

    Parameters
    ----------
    lam : int
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
                       float((2 * lam + dimensions - 2) * (2 * lam + dimensions - 4)))

    return a_lam_mu
