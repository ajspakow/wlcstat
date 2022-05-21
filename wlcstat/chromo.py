r"""Kinked WLC with twist and fixed ends

This module calculates statistics for a series of worm-like chains with twist (DNA linkers) connected
by kinks imposed by nucleosomes. Calculations include R^2, Kuhn length, propogator matrices,
full Green's function for the end-to-end distance of the polymer, and looping statistics.

Code written by Bruno Beltran, Deepti Kannan, and Andy Spakowitz

"""
import numpy as np
import scipy as sp
from numba import jit
import matplotlib.pyplot as plt
from wlcstat.util import wignerD as wd

# Constants for bare DNA properties

naked_dna_length_per_base = 0.332               # length of DNA per base pair (nm)
default_lt = 100 / naked_dna_length_per_base    # twist persistance length of DNA in bp: 100 nm, or ~301 bp
default_lp = 50 / naked_dna_length_per_base     # bend persistance length of DNA in bp: 50 nm, or ~150 bp
default_w_in = 73                               # number of bound base pairs in the entry side of dyad axis base in bp
default_w_out = 73                              # number of bound base pairs in the exit side of dyad axis base in bp
default_Lw = default_w_in + default_w_out

naked_dna_twist_density = 10.5
naked_dna_length_per_base = 0.332
naked_dna_thickness = 2.0
nucleosome_bound_dna_twist_density = 10.17

dna_params = {'tau_n': 2 * np.pi / nucleosome_bound_dna_twist_density,
              'tau_d': 2 * np.pi / naked_dna_twist_density,
              'lpb': naked_dna_length_per_base,
              'r_dna': naked_dna_thickness / 2}

helix_params_reported_richmond_davey_2003 = {
        'T': 1.67, 'b': 133.6, 'c': 41.181, 'r': 41.9,
        'phi': 0, 'psi0': 0, 'theta': 0, 'x0': 0, 'y0': 0, 'z0': 0}

helix_params_best = {'r': 4.1899999999999995, 'c': 4.531142964071856, 'T': 1.8375, 'b': 147}
#helix_params_best = {'r': 4.119206148106671, 'c': 5.1680650307039455, 'T': 1.9833333333333336, 'b': 147}

# Parameters for Fourier inversions

Klin = np.linspace(0, 10**5, 20000)
Klog = np.logspace(-3, 5, 10000)
Kvals = np.unique(np.concatenate((Klin, Klog)))
#convert to little k -- units of inverse bp (this results in kmax = 332)
kvals = Kvals / (2 * default_lp)
"""Good values to use for integrating our Green's functions. If the lp of the
bare chain under consideration changes, these should change."""
TOL = 10e-14


def r2_kinked_twlc(links, lt=default_lt, lp=default_lp, kd_unwrap=None, w_ins=default_w_in,
                   w_outs=default_w_out, tau_d=dna_params['tau_d'], tau_n=dna_params['tau_n'],
                   lmax=2, helix_params=helix_params_best, unwraps=None, random_phi=False):
    r"""Calculate the mean squared end-to-end distance, or :math:`\langle{R^2}\rangle` of a kinked WLC with a given set of
     linkers and unwrapping amounts.

    Parameters
    ----------
    links : (L,) array-like
        linker length in bp
    w_ins : float or (L+1,) array_like
        amount of DNA wrapped on entry side of central dyad base in bp
    w_outs : float or (L+1,) array_like
        amount of DNA wrapped on exit side of central dyad base in bp
    tau_n : float
        twist density of nucleosome-bound DNA in rad/bp
    tau_d : float
        twist density of naked DNA in rad/bp
    lt : float
        twist persistence length in bp
    lp : float
        DNA persistence length in bp
    lmax : int
        maximum eigenvalue l for which to compute wigner D' (default lmax = 2)

    Returns
    -------
    r2 : (L,) array-like
        mean square end-to-end distance of kinked chain as a function of chain length in nm^2
    ldna : (L,) array-like
        mean square end-to-end distance of kinked chain as a function of chain length in nm
    kuhn : float
        Kuhn length as defined by :math:`\langle{R^2}\rangle / R_{max}` in long chain limit
    """
    b = helix_params['b']
    num_linkers = len(links)
    num_nucleosomes = num_linkers + 1
    # resolve kd_unwrap
    if kd_unwrap is not None:
        sites_unbound_left = scipy.stats.binom(7, kd_unwrap).rvs(num_nucleosomes)
        sites_unbound_right = scipy.stats.binom(7, kd_unwrap).rvs(num_nucleosomes)
        w_ins, w_outs = resolve_wrapping_params(sites_unbound_left + sites_unbound_right,
                w_ins, w_outs, num_nucleosomes, unwrap_is='sites')
    else:
        w_ins, w_outs = resolve_wrapping_params(unwraps, w_ins, w_outs, num_nucleosomes)

    # calculate unwrapping amounts based on w_ins and w_outs
    mu_ins = (b - 1)/2 - w_ins
    mu_outs = (b - 1)/2 - w_outs
    # only need one g matrix per linker length, no need to recalculate each time
    # perhaps we tabulate all g's and all M's and then mix and match to grow chain?
    # for now, build dictionary of (link, wrapping) -> [B0, B1, B2]
    bmats = {}

    # B0-2curr will keep track of the B matrices as they propogate along the chain
    # initialize based on very first linker in chain
    link = mu_outs[0] + links[0] + mu_ins[1]
    wrap = w_outs[0] + w_ins[1]
    key = (link, wrap)
    R = OmegaE2E(wrap, tau_n=tau_n)
    # recall that our OmegaE2E matrix is designed to be applied from the right
    # so in order to add an arbitrary twist *before* the action of the
    # nucleosome (as if from changing the linker length) then we should apply
    # Rz to the left of R so that when the combined R is applied on the *right*
    # then the extra Rz is applied "first".
    # in this code, we use "(-gamma, -beta, -alpha)" from the left as a proxy
    # from right multiplication in build_B_matrices_for_R2
    if random_phi:
        R = Rz(2*np.pi*np.random.rand()) @ R
    alpha, beta, gamma = zyz_from_matrix(R)
    bmats[key] = build_B_matrices_for_R2(link, alpha, beta, gamma, lt, lp, tau_d, lmax)
    B0curr, B1curr, B2curr = bmats[key]

    # calculate R^2 as a function of number of nucleosomes (r2[0] is 0 nucleosomes)
    r2 = np.zeros((num_linkers,))
    lengthDNA = np.zeros_like(r2)
    r2[0] = 3 * np.real(B2curr[0,0]/B0curr[0, 0])
    lengthDNA[0] = link

    # recursively calculate Nth propagator using B matrices
    for i in range(1, num_linkers):
        # add up the effective linker lengths including unwrapping
        link = mu_outs[i] + links[i] + mu_ins[i+1]
        # w_ins[i+1] because the ith linker is between i, and i+1 nucs
        wrap = w_outs[i] + w_ins[i+1]
        key = (link, wrap)
        # update dictionary for this linker and wrapping amount, if necessary
        if key not in bmats:
            R = OmegaE2E(wrap, tau_n=tau_n)
            if random_phi:
                R = Rz(2*np.pi*np.random.rand()) @ R
            alpha, beta, gamma = zyz_from_matrix(R)
            bmats[key] = build_B_matrices_for_R2(link, alpha, beta, gamma, lt, lp, tau_d, lmax)
        B0next, B1next, B2next = bmats[key]

        # propogate B0curr, B1curr, B2curr matrices by a linker
        B0temp = B0next@B0curr
        B1temp = B1next@B0curr + B0next@B1curr
        B2temp = B2next@B0curr + 2*B1next@B1curr + B0next@B2curr
        B0curr = B0temp
        B1curr = B1temp
        B2curr = B2temp

        # Rz^2 is B2[0,0], so multiply by 3 to get R^2
        # divide by B0[0,0] to ensure normalization is OK (shouldn't matter, since B0[0,0] is 1)
        r2[i] = 3*np.real(B2curr[0,0]/B0curr[0, 0])
        lengthDNA[i] = lengthDNA[i-1] + link

    # take absolute value of R2, convert to nm^2
    r2 = np.abs(r2) * (dna_params['lpb'])**2
    lengthDNA = lengthDNA * dna_params['lpb']

    # Find scaling of R2 with length of chain at larger length scales to calculate Kuhn length
    try:
        min_i = np.round(len(lengthDNA)*5/6).astype(int)
        kuhn = stats.linregress(lengthDNA[min_i:], r2[min_i:])[0]
    except:
        kuhn = np.nan

    return r2, lengthDNA, kuhn, w_ins, w_outs



def build_B_matrices_for_R2(link, alpha, beta, gamma, lt=default_lt, lp=default_lp, tau_d=dna_params['tau_d'], lmax=2):
    r"""Helper function to construct propogator B matrices for a single linker length link with kink
    rotation given by alpha, beta, gamma. Returns the following three matrices:

    .. math::

        B^{(n)} = \lim_{k \to 0}\frac{d^n B}{dk^n}

    for n = 0, 1, and 2. where B is defined as

    .. math::

        B^{l_f j_f}_{l_0 j_0} = \sqrt{\frac{8\pi^2}{2l_f+1}} {\mathcal{D}}^{j_f j_0}_{l_f}(-\gamma, -\beta, -\alpha) g^{j_0}_{l_f l_0}
        B^{(n)}[I_f, I_0] = M[I_f, I_0] * g^{(n)}[I_f, I_0]

    The g matrix is used for the linker propogator, and the M matrix represents the rotation due to the kink.
    All matrices are super-indexed by B[If, I0] where :math:`I(l, j) = l^2 + l + j`. :math:`I` can take on :math:`l^2+2l+1` possible values.

    Notes
    -----
    Andy adds 1 to the above formula for :math:`I` since his script is in Matlab, which is
    one-indexed.

    Parameters
    ----------
    link : float
        linker length in bp
    lt : float
        twist persistence length in bp
    lp : float
        DNA persistence length in bp
    tau_d : float
        twist density of naked DNA in rad/bp
    lmax : int
        maximum eigenvalue l for which to compute wigner D' (default lmax = 2)

    Returns
    -------
    [B0, B1, B2] : (3,) list
        list of 3 matrices, each of dimension :math:`l_{max}^2 + 2l_{max} + 1`
    """

    # corresponds to If, I0 indices in Andy's notes: for every l, (2l+1) possible values of j
    ntot = lmax**2 + 2*lmax + 1
    # so matrix elements look like g[If,I0] where I0 = l0**2 + l0 + j0 and If = lf**2 + lf + jf
    # Note that for linker propogators (g matrices), jf always equals j0 (only perturbs l)
    # and for kink propogators (M matrices), lf always equals l0 (only perturbs j)
    # NOTE: for python, need to subtract 1 from Andy's formulas for indexing to work

    # g0 represents the 0th derivative of g with respect to k in the limit as k goes to 0
    g0 = np.zeros((ntot, ntot), 'complex')
    g1 = np.zeros_like(g0)
    g2 = np.zeros_like(g0)
    M  = np.zeros_like(g0)
    mywd = wd.wigner_d_vals()

    # define useful lambda functions of l and j used to compute eigenvalues, matrix elements
    I   = lambda l, j: l**2 + l + j # indexing
    al  = lambda l, j: np.sqrt((l-j)*(l+j)/(4*l**2 - 1)) # ladder coefficients alpha
    lam = lambda l, j: (l*(l+1))/(2*lp) + 0.5*((1/lt)-(1/lp))*j**2 - 1j*tau_d*j # eigenvalue of H0

    # build g and M matrices by looping over l0 and j0
    for l0 in range(lmax+1):
        for j0 in range(-l0, l0+1):
            # for this particular tuple (l0, j0), compute the relevant index in the g matrix:
            I0 = I(l0, j0)

            # Compute relevant values of lambda_lj and alpha_lj to construct g0, g1, g2
            laml0   = lam(l0, j0)
            laml0p1 = lam(l0+1, j0)
            laml0m1 = lam(l0-1, j0)
            laml0p2 = lam(l0+2, j0)
            laml0m2 = lam(l0-2, j0)
            all0    = al(l0, j0)
            all0p1  = al(l0+1, j0)
            # NOTE: this will produce nans for (l0, j0) = (2, 2) and (2, -2), but
            # this quantity isn't used for those values of l0, j0 so no problem
            if (l0 != 2 or abs(j0) != 2):
                all0m1  = al(l0-1, j0)
            all0p2  = al(l0+2, j0)

            ### Construct g0 matrix###

            # answer for g0 says l = l0 due to delta function, so g0 is diagonal
            g0[I0, I0] = np.exp(-laml0*link)

            #### Construct g1 matrix###

            # first consider case where l = l0-1
            l = l0 - 1
            If = I(l, j0)
            # check out of bounds for l, ensure l does not exceed j0 (because j = j0)
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)):
                g1[If, I0] = (1j*all0/(laml0m1 - laml0))*(np.exp(-laml0*link) - np.exp(-laml0m1*link))

            # next consider the case where l = l0+1
            l = l0 + 1
            If = I(l, j0)
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)):
                g1[If, I0] = (1j*all0p1/(laml0p1 - laml0))*(np.exp(-laml0*link) - np.exp(-laml0p1*link))

            #### Construct g2 matrix###

            # Case 1: l = l0 + 2
            l = l0 + 2
            If = I(l, j0)
            # only valid case is when (l0, j0) = (0, 0)
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)):
                g2[If, I0] = -2*all0p1*all0p2*(np.exp(-laml0*link)/((laml0p1 - laml0)*(laml0p2 - laml0)) +
                                              np.exp(-laml0p1*link)/((laml0 - laml0p1)*(laml0p2 - laml0p1)) +
                                              np.exp(-laml0p2*link)/((laml0 - laml0p2)*(laml0p1 - laml0p2)))

            # Case 2: l = l0 --- diagonal entries of g2 matrix
            l = l0
            If = I(l, j0)
            # Case 2A: terms with l0 + 1
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)):
                g2[If, I0] = -2*all0p1**2*(link*np.exp(-laml0*link)/(laml0p1 - laml0) +
                                          (np.exp(-laml0p1*link) - np.exp(-laml0*link))/(laml0 - laml0p1)**2)
            # Case 2B: terms with l0 - 1
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)+1):
                g2[If, I0] += -2*all0**2*(link*np.exp(-laml0*link)/(laml0m1 - laml0) +
                                          (np.exp(-laml0m1*link) - np.exp(-laml0*link))/(laml0 - laml0m1)**2)

            # Case 3: l = l0 - 2
            l = l0 - 2
            If = I(l, j0)
            # only valid case is when (l0, j0) = (2, 0)
            if (l >= 0) and (l <= lmax) and (l >= np.abs(j0)):
                g2[If, I0] = -2*all0*all0m1*(np.exp(-laml0*link)/((laml0m1 - laml0)*(laml0m2 - laml0)) +
                                              np.exp(-laml0m1*link)/((laml0 - laml0m1)*(laml0m2 - laml0m1)) +
                                              np.exp(-laml0m2*link)/((laml0 - laml0m2)*(laml0m1 - laml0m2)))

            # Next build M matrix
            for jf in range(-l0, l0+1):
                If = I(l0, jf)
                M[If, I0] = mywd.get(l0, jf, j0, -gamma, -beta, -alpha) / mywd.normalize(l0, jf, j0)

    B0 = M@g0
    B1 = M@g1
    B2 = M@g2
    return [B0, B1, B2]


def minimum_energy_no_sterics_linker_only(links, *, w_ins=default_w_in,
        w_outs=default_w_out, tau_n=dna_params['tau_n'],
        tau_d=dna_params['tau_d'], lpb=dna_params['lpb'],
        helix_params=helix_params_best, unwraps=None, random_phi=None):
    """Calculate the theoretical conformation of a chain of nucleosomes
    connected by straight linkers (ground state) with no steric exclusion.

    Parameters
    ----------
    links : (L,) array_like
        lengths of each linker segment
    tau_n : float
        twist density of nucleosome-bound DNA in rad/bp
    tau_d : float
        twist density of linker DNA in rad/bp
    w_ins : float or (L+1,) array_like
        amount of DNA wrapped on entry side of central dyad base
    w_outs : float or (L+1,) array_like
        amount of DNA wrapped on exit side of central dyad base

    Returns
    -------
    entry_rots : (L+1,3,3) np.ndarray
        the orientation of the material normals (as columns) of the DNA at the
        entry site
    entry_pos : (L+1,3) np.ndarray
        the position of the entry site of each nucleosome (where the first bp
        is bound)
    """
    b = helix_params['b']
    num_linkers = len(links)
    num_nucleosomes = num_linkers + 1
    w_ins, w_outs = resolve_wrapping_params(unwraps, w_ins, w_outs, num_nucleosomes)
    # initialize to valid orientation matrix
    entry_rots = np.tile(np.identity(3), (num_nucleosomes, 1, 1))
    # and to start at the origin
    entry_pos = np.zeros((num_nucleosomes, 3))
    for i in range(num_linkers):
        # w_ins[i+1] because the ith linker is between i, and i+1 nucs
        Onext = OmegaNextEntry(links[i], tau_n=tau_n, tau_d=tau_d,
                w_in=w_ins[i+1], w_out=w_outs[i], helix_params=helix_params)
        if random_phi is not None:
            Onext = Rz(random_phi*np.random.rand()) @ Onext
        entry_rots[i+1] = entry_rots[i] @ Onext
        exit_u = entry_rots[i+1,:,2]
        exit_u = exit_u/np.linalg.norm(exit_u)
        mu_out = (b - 1)/2 - w_outs[i]
        mu_in = (b - 1)/2 - w_ins[i+1]
        #no translation included in this treatment for comparison to theory
        entry_pos[i+1] = entry_pos[i] + exit_u*(mu_out + links[i] + mu_in)*lpb
    return entry_rots, entry_pos


def resolve_wrapping_params(unwraps, w_ins=None, w_outs=None, N=None, unwrap_is='bp'):
    """Allow the user to specify either one value (and tile appropriately) or
    an array of values for the number of base pairs bound to the nucleosome.
    Also allow either the number of base pairs bound on each side of the dyad to be
    specified or the number of base pairs unwrapped relative to the crystal
    structure (in which case the unwrapping is split evenly on either side of
    the dyad for simplicity.

    Parameters
    ----------
    unwraps : float or (N,) array_like
        total amount of unwrapped DNA on both sides
    w_ins (optional): float or (N,) array_like
        wrapped DNA on entry side of dyad axis
    w_outs (optional): float or (N,) array_like
        wrapped DNA on exit side of dyad axis
    N (optional): int
        output size, if other params are not array_like
    unwrap_is : string
        'bp' or 'sites', to specify whether we're counting the number of bp
        bound or the number of nucleosome-to-dna contacts (respectively)

    Returns
    -------
    w_in : (N,) np.ndarray
        N output wrapping lengths on entry side of dyad axis
    w_out : (N,) np.ndarray
        N output wrapping lengths on exit side of dyad axis

    All functions take w_in, w_out directly now."""
    if (w_ins is None) != (w_outs is None):
        raise ValueError("Either none or both of w_in and w_out must be specified.")

    if unwraps is not None:
        unwraps = np.atleast_1d(unwraps)
        w_ins, w_outs = zip(*map(functools.partial(resolve_unwrap, unwrap_is=unwrap_is), unwraps))
    w_ins = np.atleast_1d(w_ins)
    w_outs = np.atleast_1d(w_outs)
    if len(w_ins) == 1 and N is not None:
        w_ins = np.tile(w_ins, (N,))
        w_outs = np.tile(w_outs, (N,))
    return w_ins, w_outs


def OmegaNextAssoc(Ll, *, w_in=default_w_in, w_out=default_w_out,
              tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
              helix_params=helix_params_best):
    r"""Nucleosome entry to exit (including unwrapped DNA on both sides)."""
    b = helix_params['b']
    mu_in = (b - 1)/2 - w_in
    mu_out = (b - 1)/2 - w_out
    return Rz(mu_in*tau_d) @ OmegaE2E(w_in + w_out, tau_n=tau_n,
            helix_params=helix_params) @ Rz((Ll + mu_out)*tau_d)

def OmegaNextExit(Ll, *, w_in=default_w_in, w_out=default_w_out,
        tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
        helix_params=helix_params_best):
    r"""Nucleosome exit to exit."""
    b = helix_params['b']
    mu_in = (b - 1)/2 - w_in
    mu_out = (b - 1)/2 - w_out
    OE2E = OmegaE2E(w_in + w_out, tau_n=tau_n, helix_params=helix_params)
    return Rz((mu_out + Ll + mu_in)*tau_d) @ OE2E

def OmegaNextEntry(Ll, *, w_in=default_w_in, w_out=default_w_out,
        tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
        helix_params=helix_params_best):
    r"""Nucleosome entry to entry."""
    b = helix_params['b']
    mu_in = (b - 1)/2 - w_in
    mu_out = (b - 1)/2 - w_out
    OE2E = OmegaE2E(w_in + w_out, tau_n=tau_n, helix_params=helix_params)
    return OE2E @ Rz((mu_out + Ll + mu_in)*tau_d)


def OmegaE2E(Lw=default_Lw, *, tau_n=dna_params['tau_n'], helix_params=helix_params_best):
    r"""Nucleosome core entry to exit rotation matrix

    :math:`\Omega_C \cdot R_z(L_w \left[ \tau_n - \tau_H \right])`

    Should be multiplied against entry orientation from the right to get the
    exit orientation.

    Example
    -------
    To get the exit orientation from the nucleosome assuming the bound DNA
    adopts a configuration with zero intrinsic twist is just
        >>> import nuc_chain.geometry as ncg
        >>> Omega0 = np.identity(3)
        >>> Omega0Exit = Omega0 @ ncg.OmegaE2E(tau_n=0)

    Parameters
    ----------
    tau_n : float
        the twist density of nucleosome-bound dna in radians/bp
    Lw : float
        Length of DNA bound to the nucleosome (#bound bp - 1)
    helix_params : Dict[str, float]
        parameters to :py:func:`H` family of function

    Returns
    -------
    OmegaE2E : (3,3) np.ndarray[float64]
        Entry to exit rotation matrix
    """
    # N basepairs -> N-1 inter-base-pair segments
    # in other words, no +1 since e.g. if w_in=w_out=1, then there are three
    # base pairs but the dna "length" is 2, not 3.
    theta = Lw * (tau_n - Htau(**helix_params))
    Oc = OmegaC(Lw, helix_params=helix_params)

    return Oc @ Rz(theta)


def OmegaC(Lw, helix_params=helix_params_best):
    """inv(H_rot(0)) @ H_rot(Lw)

    The non-twist-corrected rotation matrix that takes the n, b, and u vectors
    of the helix from the start of the helix to the exit.

    Parameters
    ----------
    Lw : float
        Length of DNA bound to the nucleosome (#bound bp - 1)
    helix_params : (optional) Dict[str, float]

    Returns
    -------
    Omega_next : (3, 3) np.ndarray[float64]
        Rotation matrix taking the input to output orientations of the normal,
        binormal, tangent triad of the helix.
    """
    I = np.identity(3)
    # since the helix has uniform torsion and curvature, this calculation is
    # invariant to translation of the parametric coordinate (could just as
    # easily use x, x + Lw for any x)
    entry_inv = np.linalg.solve(H_rot(0, **helix_params), I)
    exit_rotation = H_rot(Lw, **helix_params)
    return entry_inv@exit_rotation


def H(i, r, c, T, b):
    """Parametric description of a simple left-handed helix aligned with
    :math:`\hat{z}`.

    The helix is parametrized on [0,1]. Has height c, radius r, and undergoes T
    twists in that range. Centerline is the positive z-axis and starting point
    is (r,0,0).

    Parameters
    ----------
    ${helix_doc}

    Returns
    -------
    X : (3,N) np.ndarray
        Position coordinate for each i.
    """
    if i is None:
        i = np.arange(b)
    i = np.atleast_1d(i)
    phase = (2*np.pi*i*T)/b
    return np.array([r*np.cos(-phase), r*np.sin(-phase), c*(i/b)])

def Hu(i, r, c, T, b, normed=True):
    r"""Tangent vector to :py:func:`H`, :math:`\frac{dH}{di}`.

    Parameters
    ----------
    ${helix_doc}
    ${normed_doc}

    Returns
    -------
    U : (3,N) np.ndarray
        Tangent vector for each i.
    """
    if i is None:
        i = np.arange(b)
    i = np.atleast_1d(i)
    p = 2*np.pi*T/b
    U = np.array([p*r*np.sin(-p*i),
                 -p*r*np.cos(-p*i),
                  c/b*np.ones_like(i)])
    if normed:
        return normalize_HX(U)
    return U

def Hn(i, r, c, T, b, normed=True):
    r"""Normal vector to :py:func:`H`. The unit vector in the direction of
    :math:`\frac{d^2 H}{di^2}`.

    Parameters
    ----------
    ${helix_doc}
    ${normed_doc}

    Returns
    -------
    N : (3,N) np.ndarray
        Normal vector for each i.
    """
    if i is None:
        i = np.arange(b)
    i = np.atleast_1d(i)
    p = 2*np.pi*T/b
    N = np.array([-p*p*r*np.cos(-p*i),
                  -p*p*r*np.sin(-p*i),
                  np.zeros_like(i)])
    if normed:
        return normalize_HX(N)
    return N

def Hb(*args, **kwargs):
    r"""Bi-normal vector to :py:func:`H`. Unit vector in the direction of
    :math:`\frac{dH}{di}\times\frac{d^2
    H}{di^2}`.

    Parameters
    ----------
    ${helix_doc}
    ${normed_doc}

    Returns
    -------
    B : (3,N) np.ndarray
        Bi-normal vector for each i.
    """
    U = Hu(*args, **kwargs)
    N = Hn(*args, **kwargs)
    B = np.cross(U, N, axis=0)
    if 'normed' in kwargs and kwargs['normed']:
        return normalize_HX(B)
    return B

def Hk(r, c, T, b):
    r"""Curvature of our standard helix.

    Parameters
    ----------
    ${helix_params_doc}

    Returns
    -------
    k : float64
        Curvature of the helix.

    Notes
    -----
    The curvature of a helix is a constant, in our parameterization it works
    out to :math:`\frac{r}{\gamma}\left(\frac{2\pi T}{b}\right)^2`, where
    :math:`\gamma = (1/b)\sqrt{4\pi^2r^2T^2 + c^2}`.
    """
    psi = np.power(2*np.pi*T, 2)
    c2 = c*c
    r2 = r*r
    return r*psi/np.sqrt(psi*r2 + c2)/b

def Htau(r, c, T, b):
    r"""Torsion of our standard helix.

    Parameters
    ----------
    ${helix_params_doc}

    Returns
    -------
    tau : float64
        Torsion of the helix.

    Notes
    -----
    The torsion of a helix is a constant, in our parameterization it works out
    to :math:`\frac{-2\pi c T}{b} \left[4\pi^2r^2T^2 + c^2\right]^{-1/2}`.
    """
    return -2*np.pi*c*T/np.sqrt(np.power(2*np.pi*r*T, 2) + c*c)/b

def H_rot(i, r, c, T, b):
    """OmegaH(s). Given OmegaH(0) = [Hn, Hb, Hu].

    A rotation element, so det(H_rot) == 1."""
    return np.concatenate([Hn(i, r, c, T, b),
            Hb(i, r, c, T, b), Hu(i, r, c, T, b)], axis=1)

def H_rot_oriented(*args, **kwargs):
    """OmegaH(s) given OmegaH(0) arbitrary.

    A rotation element, so det(H_rot) == 1.
    For some reason, tends to produce matrices that are numerically off of
    H_rot by a factor of like 1-2% in each element for equivalent parameters?
    Off by one somewhere? That's the order of magnitude of a bp."""
    return np.concatenate([Hn_oriented(*args, **kwargs),
            Hb_oriented(*args, **kwargs), Hu_oriented(*args, **kwargs)], axis=1)

def H_oriented(*args, **kwargs):
    """Position a helix at entry_pos and entry_rot, as returned by e.g.
    minimum_energy_no_sterics.

    Parameters
    ----------
    ${helix_doc}
    ${orientation_doc}

    Returns
    -------
    X : (3,N) np.ndarray
        Position coordinate for each i.
    """
    return HX_at_pos_rot_(H, *args, **kwargs)

def Hu_oriented(*args, **kwargs):
    """Tangent to :py:func:`H_oriented` at entry_rot.

    Parameters
    ----------
    ${helix_params_doc}
    ${entry_rot_doc}

    Returns
    -------
    U : (3,N) np.ndarray
        Tangent vector for each i.
    """
    return HX_at_pos_rot_(Hu, *args, **kwargs)

def Hn_oriented(*args, **kwargs):
    """Normal to :py:func:`H_oriented` at entry_rot.

    Parameters
    ----------
    ${helix_params_doc}
    ${entry_rot_doc}

    Returns
    -------
    n : (3,N) np.ndarray
        normal vector for each i.
    """
    return HX_at_pos_rot_(Hn, *args, **kwargs)


def normalize_HX(X):
    """Normalize a list of 3-vectors of shape (N,3).

    Parameters
    ----------
    X : (N, 3) array_like
        N input vectors to be normalized.

    Returns
    -------
    U : (N, 3) np.ndarray
        N outputs, such that np.linalg.norm(U[i,:]) == 1 for i in range N.
    """
    return X/np.linalg.norm(X, axis=0)[None,:]


def Rx(theta):
    r"""Rotation matrix about the x axis with angle theta.

    Notes
    -----
    ..math::

        \frac{1}{\sqrt{3}}
        \begin{bmatrix}
            1 &             0 &              0\\
            0 & np.cos(theta) & -np.sin(theta)\\
            0 & np.sin(theta) &  np.cos(theta)
        \end{bmatrix}
    """
    return np.array([[1,             0,              0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta),  np.cos(theta)]])


def Ry(theta):
    r"""Rotation matrix about the z axis with angle theta.

    Notes
    -----
    ..math::

        \frac{1}{\sqrt{3}}
        \begin{bmatrix}
            np.cos(theta) & 0 & -np.sin(theta) \\
                        0 & 1 &              0 \\
            np.sin(theta) & 0 &  np.cos(theta)
        \end{bmatrix}
    """
    return np.array([[np.cos(theta), 0, -np.sin(theta)],
                     [            0, 1,              0],
                     [np.sin(theta), 0,  np.cos(theta)]])


def Rz(theta):
    r"""Rotation matrix about the z axis with angle theta.

    Notes
    -----
    ..math::

        \frac{1}{\sqrt{3}}
        \begin{bmatrix}
            np.cos(theta) & -np.sin(theta) & 0 \\
            np.sin(theta) &  np.cos(theta) & 0 \\
                        0 &             0  & 1
        \end{bmatrix}
    """
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta),  np.cos(theta), 0],
                    [            0,             0,  1]])


def zyz_from_matrix(R):
    """Convert from a rotation matrix to (extrinsic zyz) Euler angles.

    Parameters
    ----------
    R : (3, 3) array_like
        Rotation matrix to convert

    Returns
    -------
    alpha : float
        The final rotation about z''
    beta : float
        The second rotation about y'
    gamma : float
        the initial rotation about z

    Notes
    -----
    See pg. 13 of Bruno's (NCG) nucleosome geometry notes for a detailed
    derivation of this formula. Easy to test correctness by simply
    reversing the process R == Rz(alpha)@Ry(beta)@Rz(gamma).
    """
    beta = np.arccos(R[2,2])
    if np.abs(1 - np.abs(R[2,2])) < TOL:
        beta = 0
        gamma = 0
        Rz_alpha = R
    else:
        gamma = np.arctan2(-R[2,1]/np.sin(beta), R[2,0]/np.sin(beta))
        Rz_alpha = R @ np.linalg.inv(Rz(gamma)) @ np.linalg.inv(Ry(beta))
    alpha = np.arccos(Rz_alpha[0,0])
    # # couldn't get this part of the formula to work for some reason
    # alpha = np.arctan2(R[0,2]/np.sin(beta), -R[1,2]/np.sin(beta))
    return alpha, beta, gamma

def gen_chromo_conf(links, lt=default_lt, lp=default_lp, kd_unwrap=None, w_ins=default_w_in,
             w_outs=default_w_out, tau_d=dna_params['tau_d'], tau_n=dna_params['tau_n'],
             lpb=dna_params['lpb'], r_dna=dna_params['r_dna'],
             helix_params=helix_params_best, unwraps=None, random_phi=False):
    """
    Generate DNA and nucleosome conformation based on chain growth algorithm

    Parameters
    ----------
    links : (L,) array-like
        linker length in bp
    w_ins : float or (L+1,) array_like
        amount of DNA wrapped on entry side of central dyad base in bp
    w_outs : float or (L+1,) array_like
        amount of DNA wrapped on exit side of central dyad base in bp
    tau_n : float
        twist density of nucleosome-bound DNA in rad/bp
    tau_d : float
        twist density of naked DNA in rad/bp
    lt : float
        twist persistence length in bp
    lp : float
        DNA persistence length in bp

    Returns
    -------

    """


    # Put the parameters into an appropriate form
    nbpnuc = helix_params['b']
    num_linkers = len(links)
    num_nuc = num_linkers + 1
    hnuc = helix_params['c'] / 2
    rnuc = helix_params['r']
    ltnuc = np.sqrt(hnuc ** 2 + (2 * np.pi * rnuc) ** 2)
    eps = lp / lpb
    epst = lt / lpb
    om = tau_d
    omnuc = tau_n
    omdna = 0.75 * np.pi

    # resolve kd_unwrap
    if kd_unwrap is not None:
        sites_unbound_left = scipy.stats.binom(7, kd_unwrap).rvs(num_nuc)
        sites_unbound_right = scipy.stats.binom(7, kd_unwrap).rvs(num_nuc)
        w_ins, w_outs = resolve_wrapping_params(sites_unbound_left + sites_unbound_right,
                w_ins, w_outs, num_nucleosomes, unwrap_is='sites')
    else:
        w_ins, w_outs = resolve_wrapping_params(unwraps, w_ins, w_outs, num_nuc)
    bounds = w_ins + w_outs + 1

    # Initialize the conformation
    num_bp_total = np.sum(links) + np.sum(bounds)
    r = np.zeros((num_bp_total, 3))
    rdna1 = np.zeros((num_bp_total, 3))
    rdna2 = np.zeros((num_bp_total, 3))
    t1 = np.zeros((num_bp_total, 3))
    t2 = np.zeros((num_bp_total, 3))
    t3 = np.zeros((num_bp_total, 3))
    rn = np.zeros((num_nuc, 3))
    un = np.zeros((num_nuc, 3))

    # Generate the conformation using the chain-growth algorithm
    t10 = np.array([1, 0, 0])
    t30 = np.array([0, 0, 1])
    t20 = np.cross(t30, t10)
    r0 = np.array([0, 0, 0])

    count = 0
    for inuc in range(num_nuc):
        # Generate the nucleosomal dna
        bound = bounds[inuc]
        n1 = -t10
        n3 = (2 * np.pi * rnuc / ltnuc) * t20 + (hnuc / ltnuc) * t30
        n2 = np.cross(n3, n1)

        delta = - rnuc * n1 + (hnuc * (bound - 1) * lpb / (2 * ltnuc)) * n3
        rn[inuc, :] = r0 + delta
        un[inuc, :] = n3

        for i in range(bound):
            s = i * lpb
            r[count, :] = (rnuc * np.cos(2 * np.pi * s / ltnuc) * n1 +
                           rnuc * np.sin(2 * np.pi * s / ltnuc) * n2 +
                           (hnuc * s / ltnuc) * n3 + r0 - rnuc * n1)
            if i == 0:
                t1[count, :] = t10
                t2[count, :] = t20
                t3[count, :] = t30
            else:
                t3[count, :] = (- (2 * np.pi * rnuc / ltnuc) * np.sin(2 * np.pi * s / ltnuc) * n1
                                + (2 * np.pi * rnuc / ltnuc) * np.cos(2 * np.pi * s / ltnuc) * n2
                                + hnuc / ltnuc * n3)
                t3[count, :] /= np.linalg.norm(t3[count, :])
                th = np.arccos(np.dot(t3[count, :], t3[count - 1, :]))
                phi = np.arctan2(np.dot(t3[count, :], t2[count - 1, :]), np.dot(t3[count, :], t1[count - 1, :]))
                psi = -phi + omnuc

                t1p = (np.cos(th) * np.cos(phi) * t1[count - 1, :]
                       + np.cos(th) * np.sin(phi) * t2[count - 1, :]
                       - np.sin(th) * t3[count - 1, :])
                t1p -= np.dot(t3[count, :], t1p) * t3[count, :]
                t1p /= np.linalg.norm(t1p)
                t2p = np.cross(t3[count, :], t1p)
                t1[count, :] = np.cos(psi) * t1p + np.sin(psi) * t2p
                t2[count, :] = np.cross(t3[count, :], t1[count, :])

            rdna1[count, :] = r[count, :] + t1[count, :] * r_dna
            rdna2[count, :] = r[count, :] + r_dna * (np.cos(omdna) * t1[count, :] +
                                                     np.sin(omdna) * t2[count, :])
            count += 1

        # Calculate the position and orientation heading into the next linker

        th = np.arccos(1 / eps * np.log(
            np.random.uniform() * 2 * np.sinh(eps) + np.exp(-eps)))
        phi = 2 * np.pi * np.random.uniform()
        psi = -phi + om + np.random.normal() / np.sqrt(epst)

        t1p = (np.cos(th) * np.cos(phi) * t1[count - 1, :]
               + np.cos(th) * np.sin(phi) * t2[count - 1, :]
               - np.sin(th) * t3[count - 1, :])
        t3p = (np.sin(th) * np.cos(phi) * t1[count - 1, :]
               + np.sin(th) * np.sin(phi) * t2[count - 1, :]
               + np.cos(th) * t3[count - 1, :])
        t3p /= np.linalg.norm(t3p)
        t1p -= np.dot(t3p, t1p) * t3p
        t1p /= np.linalg.norm(t1p)
        t2p = np.cross(t3p, t1p)

        t10 = np.cos(psi) * t1p + np.sin(psi) * t2p
        t30 = t3p
        t20 = np.cross(t30, t10)
        r0 = r[count - 1, :] + t30 * lpb

        # Generate the linker dna
        if inuc < (num_nuc - 1):
            link = links[inuc]

            t1[count, :] = t10
            t2[count, :] = t20
            t3[count, :] = t30
            r[count, :] = r0
            rdna1[count, :] = r[count, :] + t1[count, :] * r_dna
            rdna2[count, :] = r[count, :] + r_dna * (np.cos(omdna) * t1[count, :] +
                                                     np.sin(omdna) * t2[count, :])
            count += 1

            for i in range(1, link):
                th = np.arccos(1 / eps * np.log(
                    np.random.uniform() * 2 * np.sinh(eps) + np.exp(-eps)))
                phi = 2 * np.pi * np.random.uniform()
                psi = -phi + om + np.random.normal() / np.sqrt(epst)

                t1p = (np.cos(th) * np.cos(phi) * t1[count - 1, :]
                       + np.cos(th) * np.sin(phi) * t2[count - 1, :]
                       - np.sin(th) * t3[count - 1, :])
                t3p = (np.sin(th) * np.cos(phi) * t1[count - 1, :]
                       + np.sin(th) * np.sin(phi) * t2[count - 1, :]
                       + np.cos(th) * t3[count - 1, :])
                t3p /= np.linalg.norm(t3p)
                t1p -= np.dot(t3p, t1p) * t3p
                t1p /= np.linalg.norm(t1p)
                t2p = np.cross(t3p, t1p)

                t1[count, :] = np.cos(psi) * t1p + np.sin(psi) * t2p
                t3[count, :] = t3p
                t2[count, :] = np.cross(t3[count, :], t1[count, :])

                r[count, :] = r[count - 1, :] + t3[count, :] * lpb
                rdna1[count, :] = r[count, :] + t1[count, :] * r_dna
                rdna2[count, :] = r[count, :] + r_dna * (np.cos(omdna) * t1[count, :] +
                                                     np.sin(omdna) * t2[count, :])
                count += 1

            # Calculate the position and orientation heading into the nucleosome
            th = np.arccos(1 / eps * np.log(
                np.random.uniform() * 2 * np.sinh(eps) + np.exp(-eps)))
            phi = 2 * np.pi * np.random.uniform()
            psi = -phi + om + np.random.normal() / np.sqrt(epst)

            t1p = (np.cos(th) * np.cos(phi) * t1[count - 1, :]
                   + np.cos(th) * np.sin(phi) * t2[count - 1, :]
                   - np.sin(th) * t3[count - 1, :])
            t3p = (np.sin(th) * np.cos(phi) * t1[count - 1, :]
                   + np.sin(th) * np.sin(phi) * t2[count - 1, :]
                   + np.cos(th) * t3[count - 1, :])
            t3p /= np.linalg.norm(t3p)
            t1p -= np.dot(t3p, t1p) * t3p
            t1p /= np.linalg.norm(t1p)
            t2p = np.cross(t3p, t1p)

            t10 = np.cos(psi) * t1p + np.sin(psi) * t2p
            t30 = t3p
            t20 = np.cross(t30, t10)

            r0 = r[count - 1, :] + t30 * lpb

    return r, rdna1, rdna2, rn, un


def gen_chromo_pymol_file(r, rdna1, rdna2, rn, un, filename='r_poly.pdb', ring=False):
    """

    Parameters
    ----------
    r
    rdna1
    rdna2
    rn
    un
    filename
    ring

    Returns
    -------

    """
    # Setup the parameters for imaging nucleosome array
    hnuc = 2.265571482035928

    # Open the file
    f = open(filename, 'w')

    atomname1 = "A1"    # Chain atom type
    atomname2 = "A2"    # Chain atom type
    atomname3 = "A3"    # Chain atom type
    atomname4 = "A4"    # Chain atom type
    resname = "SSN"     # Type of residue (UNKnown/Single Stranded Nucleotide)
    chain = "A"         # Chain identifier
    resnum = 1
    numdna = len(r[:, 0])
    numnuc = len(rn[:, 0])
    descrip = "Pseudo atom representation of DNA"
    chemicalname = "Body and ribbon spatial coordinates"

    # Write the preamble to the pymol file

    f.write('HET    %3s  %1s%4d   %5d     %-38s\n' % (resname, chain, resnum, numdna, descrip))
    f.write('HETNAM     %3s %-50s\n' % (resname, chemicalname))
    f.write('FORMUL  1   %3s    C20 N20 P21\n' % (resname))

    # Write the conformation to the pymol file

    # Define the dna centerline positions
    count = 1
    for ind in range(numdna):
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                    (count, atomname1, resname, chain, r[ind, 0], r[ind, 1], r[ind, 2], 1.00, 1.00))
        count += 1

    # Define the dna strand 1 positions
    for ind in range(numdna):
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                    (count, atomname2, resname, chain, rdna1[ind, 0], rdna1[ind, 1], rdna1[ind, 2], 1.00, 1.00))
        count += 1

    # Define the dna strand 2 positions
    for ind in range(numdna):
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                    (count, atomname3, resname, chain, rdna2[ind, 0], rdna2[ind, 1], rdna2[ind, 2], 1.00, 1.00))
        count += 1

    # Define the nucleosome positions
    for ind in range(numnuc):
        rnucind = rn[ind, :] + 0.5 * hnuc *un[ind, :]
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                    (count, atomname4, resname, chain, rnucind[0], rnucind[1], rnucind[2], 1.00, 1.00))
        count += 1
        rnucind = rn[ind, :] - 0.5 * hnuc *un[ind, :]
        f.write('ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n' %
                    (count, atomname4, resname, chain, rnucind[0], rnucind[1], rnucind[2], 1.00, 1.00))
        count += 1



    # Define the connectivity in the chain

    # Connectivity for the center beads
    count = 1
    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count + 1, count - 1 + numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count + 1))
    count += 1

    for ind in range(2, numdna):
        f.write('CONECT%5d%5d%5d\n' % (count , count - 1, count + 1))
        count += 1

    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count - 1, 1 + count - numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count - 1))
    count += 1

    # Connectivity for the dna chain 1 beads
    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count + 1, count - 1 + numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count + 1))
    count += 1

    for ind in range(2, numdna):
        f.write('CONECT%5d%5d%5d\n' % (count , count - 1, count + 1))
        count += 1

    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count - 1, 1 + count - numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count - 1))
    count += 1

    # Connectivity for the dna chain 2 beads
    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count + 1, count - 1 + numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count + 1))
    count += 1

    for ind in range(2, numdna):
        f.write('CONECT%5d%5d%5d\n' % (count , count - 1, count + 1))
        count += 1

    if ring:
        f.write('CONECT%5d%5d%5d\n' % (count, count - 1, 1 + count - numdna))
    else:
        f.write('CONECT%5d%5d\n' % (count, count - 1))
    count += 1

    # Connectivity for the nucleosome positions
    for ind in range(numnuc):
        f.write('CONECT%5d%5d\n' % (count, count + 1))
        count += 1
        f.write('CONECT%5d%5d\n' % (count, count - 1))
        count += 1

    # Close the file
    f.write('END')

    f.close()

    return