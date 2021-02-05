"""Rouse polymer, analytical results.

Notes
-----
There are two parameterizations of the "Rouse" polymer that are commonly used,
and they use the same variable name for two different things.

In one, N is the number of Kuhn lengths, and in the other, N is the number of
beads, each of which can represent an arbitrary number of Kuhn lengths.
"""
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
# import mpmath

from functools import lru_cache
from pathlib import Path
import os


@jit
def rouse_mode(p, n, N=1):
    """Eigenbasis for Rouse model.

    Indexed by p, depends only on position n/N along the polymer of length N.
    N=1 by default.

    Weber, Phys Rev E, 2010 (Eq 14)"""
    p = np.atleast_1d(p)
    phi = np.sqrt(2)*np.cos(p*np.pi*n/N)
    phi[p == 0] = 1
    return phi


@jit(nopython=True)
def rouse_mode_coef(p, b, N, kbT=1):
    """k_p: Weber Phys Rev E 2010, after Eq. 18."""
    # alternate: k*pi**2/N * p**2, i.e. k = 3kbT/b**2
    return 3*np.pi**2*kbT/(N*b**2)*p**2


@jit(nopython=True)
def kp_over_kbt(p : float, b : float, N : float):
    """k_p/(k_B T) : "non-dimensionalized" k_p is all that's needed for most
    formulas, e.g. MSD."""
    return (3*np.pi*np.pi)/(N*b*b) * p*p


@jit(nopython=True)
def linear_mid_msd(t, b, N, D, num_modes=20000):
    """
    modified from Weber Phys Rev E 2010, Eq. 24.
    """
    msd = np.zeros_like(t)
    for p in range(1, num_modes+1):
        k2p_norm = kp_over_kbt(2*p, b, N)
        msd += (1 / k2p_norm) * (1 - np.exp(-k2p_norm * (D / N) * t))
    return 12 * msd + 6 * D * t / N


@jit(nopython=True)
def gaussian_G(r, N, b):
    """Green's function of a Gaussian chain at N Kuhn lengths of separation,
    given a Kuhn length of b"""
    r2 = np.power(r, 2)
    return np.power(3/(2*np.pi*b*b*N), 3/2)*np.exp(-(3/2)*r2/(N*b*b))


@jit(nopython=True)
def gaussian_Ploop(a, N, b):
    """Looping probability for two loci on a Gaussian chain N kuhn lengths
    apart, when the Kuhn length is b, and the capture radius is a"""
    Nb2 = N*b*b
    return spycial.erf(a*np.sqrt(3/2/Nb2)) - a*np.sqrt(6/np.pi/Nb2)/np.exp(3*a*a/2/Nb2)


@jit(nopython=True)
def _cart_to_sph(x, y, z):
    r = np.sqrt(x*x + y*y + z*z)
    if r == 0.0:
        return 0.0, 0.0, 0.0
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)
    return r, phi, theta


def confined_G(r, rp, N, b, a, n_max=100, l_max=50):
    # first precompute the zeros of the spherical bessel functions, since our
    # routine to do so is pretty slow
    if l_max >= 86: # for l >= 86, |m| >= 85 returns NaN from sph_harm
        raise ValueError("l_max > 85 causes NaN's from scipy.special.sph_harm")
    if confined_G.zl_n is None or n_max > confined_G.zl_n.shape[1] \
            or l_max > confined_G.zl_n.shape[0]:
        confined_G.zl_n = spherical_jn_zeros(l_max, n_max)
    # some aliases
    spherical_jn = scipy.special.spherical_jn
    Ylm = scipy.special.sph_harm
    zl_n = confined_G.zl_n[:l_max+1,:n_max]
    # convert to spherical coordinates
    r = np.array(r)
    if r.ndim == 1:
        x, y, z = r
        xp, yp, zp = rp
    # elif r.ndim == 2:
    #     x = r[:,0]
    #     y = r[:,1]
    #     z = r[:,2]
    #     xp = rp[:,0]
    #     yp = rp[:,1]
    #     zp = rp[:,2]
    else:
        raise ValueError("Don't understand your R vectors")
    r, phi, theta = _cart_to_sph(x, y, z)
    rp, phip, thetap = _cart_to_sph(xp, yp, zp)
    # compute terms based on indexing (ij in meshgrid to match zl_n)
    l, n = np.meshgrid(np.arange(l_max+1), np.arange(n_max), indexing='ij')
    ln_term = 2*np.exp(-b**2/6 * (zl_n/a)**2 * N)
    ln_term = ln_term/(a**3 * spherical_jn(l+1, zl_n)**2)
    ln_term = ln_term*spherical_jn(l, zl_n/a*r)*spherical_jn(l, zl_n/a*rp)
    # l,m terms
    m = np.arange(-l_max, l_max+1)
    l, m = np.meshgrid(np.arange(l_max+1), np.arange(-l_max, l_max+1))
    lm_term = Ylm(m, l, phi, theta)*Ylm(m, l, phip, thetap)
    lm_mask = np.abs(m) <= l
    lm_term[~lm_mask] = 0
    # now broadcast and sum
    G = np.sum(ln_term[None,:,:] * lm_term[:,:,None])
    return G
confined_G.zl_n = None


def linear_mscd(t, D, Ndel, N, b=1, num_modes=20000):
    r"""
    Compute mscd for two points on a linear polymer.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    Ndel : float
        Distance from the last linkage site to the measured site. This ends up
        being (1/2)*separation between the loci (in Kuhn lengths).
    N : float
        The full lengh of the linear polymer (in Kuhn lengths).
    b : float
        The Kuhn length (in desired length units).
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    mscd : (M,) np.array<float>
        result
    """
    mscd = np.zeros_like(t)

    k1 = 3 * np.pi ** 2 / (N * (b ** 2))
    sum_coeff = 48 / k1
    exp_coeff = k1 * D / N
    sin_coeff = np.pi * Ndel / N

    for p in range(1, num_modes+1, 2):
        mscd += (1/p**2) * (1 - np.exp(-exp_coeff * (p ** 2) * t)) \
                * np.sin(sin_coeff*p)**2

    return sum_coeff * mscd


def ring_mscd(t, D, Ndel, N, b=1, num_modes=20000):
    r"""
    Compute mscd for two points on a ring.

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSCD.
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    Ndel : float
        (1/2)*separation between the loci on loop (in Kuhn lengths)
    N : float
        full length of the loop (in Kuhn lengths)
    b : float
        The Kuhn length, in desired output length units.
    num_modes : int
        How many Rouse modes to include in the sum.

    Returns
    -------
    mscd : (M,) np.array<float>
        result
    """
    mscd = np.zeros_like(t)

    k1 = 12 * np.pi ** 2 / (N * (b ** 2))
    sum_coeff = 48 / k1
    exp_coeff = k1 * D / N
    sin_coeff = 2 * np.pi * Ndel / N

    for p in range(1, num_modes+1):
        mscd += (1 / p ** 2) * (1 - np.exp(-exp_coeff * p ** 2 * t)) \
                * np.sin(sin_coeff * p) ** 2
    return sum_coeff * mscd


def linear_msd_confine(t, D, Nmono, N, b=1, k_conf=1, num_modes=20000):
    r"""
    Compute msd for a Rouse polymer confined within a harmonic confining potential

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSD
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    Nmono : float
        Monomer position of the tagged locus (in Kuhn length).
    N : float
        The full lengh of the linear polymer (in Kuhn lengths).
    b : float
        The Kuhn length (in desired length units).
    k_conf : float
        Confinement strength
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    msd : (M,) np.array<float>
        result
    """
    msd = np.zeros_like(t)

    cos_coeff = np.pi * Nmono / N
    tau_p = N / D / (N * k_conf)
    msd += 6 / (N * k_conf) * (1 - np.exp(-t / tau_p))

    for p in range(1, num_modes+1):
        kp = 3 * np.pi ** 2 * p ** 2 / (N * (b ** 2))
        tau_p = N / D / (kp + N * k_conf)
        msd += 12 / (kp + N * k_conf) * (1 - np.exp(-t / tau_p)) * np.cos(cos_coeff * p) ** 2

    return msd


def linear_msd_confine_plateau(Nmono, N, b=1, k_conf=1, num_modes=20000):
    r"""
    Compute msd for a Rouse polymer confined within a harmonic confining potential

    Parameters
    ----------
    t : (M,) float, array_like
        Times at which to evaluate the MSD
    D : float
        Diffusion coefficient, (in desired output length units). Equal to
        :math:`k_BT/\xi` for :math:`\xi` in units of "per Kuhn length".
    Nmono : float
        Monomer position of the tagged locus (in Kuhn length).
    N : float
        The full lengh of the linear polymer (in Kuhn lengths).
    b : float
        The Kuhn length (in desired length units).
    k_conf : float
        Confinement strength
    num_modes : int
        how many Rouse modes to include in the sum

    Returns
    -------
    msd : (M,) np.array<float>
        result
    """
    msd_plateau = np.zeros_like(k_conf)

    cos_coeff = np.pi * Nmono / N
    msd_plateau += 6 / (N * k_conf)

    for p in range(1, num_modes+1):
        kp = 3 * np.pi ** 2 * p ** 2 / (N * (b ** 2))
        msd_plateau += 12 / (kp + N * k_conf) * np.cos(cos_coeff * p) ** 2

    return msd_plateau


def end_to_end_corr(t, D, N, num_modes=10000):
    """Doi and Edwards, Eq. 4.35"""
    mscd = np.zeros_like(t)
    tau1 = N ** 2 / (3 * np.pi * np.pi * D)
    for p in range(1, num_modes+1, 2):
        mscd += 8 / p / p / np.pi / np.pi * np.exp(- t * p ** 2 / tau1)
    return N * mscd



# Parameter values for chromosome V of S. cerevisiae

chrv_length = 577  # Mb
ura_loc = 116
centr_loc = 152
ura_locus_frac = (chrv_length - ura_loc)/chrv_length
chrv_centromere_frac = (chrv_length - centr_loc)/chrv_length
col_width = 3.405
golden_ratio = (1 + np.sqrt(5))/2
chrv_size_bp = 576874
location_ura_bp = np.mean([116167, 116970])

chrv_size_kuhn = 1165
kuhn_length = 0.015  # um
chrv_size_effective_um = chrv_size_kuhn * kuhn_length
location_ura_effective_um = location_ura_bp * (chrv_size_effective_um / chrv_size_bp)

nuc_radius_um = 1.3  # Average of het5 msd convex hull distribution
sim_nuc_radius_um = 1
sim_D = 0.02  # um^2/s


def draw_cells(linkages, min_y=0.05, max_y=0.95, locus_frac=ura_locus_frac, centromere_frac=chrv_centromere_frac,
               chr_size = chrv_size_effective_um):
    r"""
    Render the model of homologous chromosomes with linkages.

    Parameters
    ----------
    linkages : float array
        List of the link positions between the homologous chromosomes


    Returns
    -------

    """

    all_closest_links = []
    linkages = [np.sort((chr_size - linkages[i]) / chr_size) for i in range(len(linkages))]

    n_cells = len(linkages)
    fig = plt.figure(figsize=(col_width, col_width/golden_ratio))
    ax = fig.add_axes([0, 0, 1, 1])
    plt.axis('off')
    n_fences = n_cells + 1
    fence_posts = np.linspace(0, 1, n_fences)
    width_per_cell = np.diff(fence_posts)[0]
    cell_centers = (fence_posts[1:] + fence_posts[:-1]) / 2
    width_to_chr_center = width_per_cell / 5
    cell_width = 15
    for i, x in enumerate(cell_centers):
        for dx in [width_to_chr_center, -width_to_chr_center]:
            plt.plot([x + dx, x + dx],
                     [min_y, max_y], transform=ax.transAxes, linewidth=cell_width,
                     solid_capstyle='round', color=[197/255, 151/255, 143/255])
            plt.scatter([x + dx], [min_y + (max_y - min_y)*centromere_frac],
                        zorder=10, transform=ax.transAxes, s=200, color='k')
            plt.scatter([x + dx], [min_y + (max_y - min_y)*locus_frac],
                        zorder=15, transform=ax.transAxes, s=500, color='g', marker='*', edgecolors='k')
        for linkage in linkages[i]:
            plt.plot([x - width_to_chr_center, x + width_to_chr_center],
                     [min_y + (max_y - min_y)*linkage, min_y + (max_y - min_y)*linkage],
                     color=(0, 0, 1), transform=ax.transAxes,
                     linewidth=5, solid_capstyle='round')
        num_linkages = len(linkages[i])
        j = np.searchsorted(linkages[i], locus_frac)
        closest_links = []
        if j != 0:
            closest_links.append(linkages[i][j - 1])
        if j != num_linkages:
            closest_links.append(linkages[i][j])
        closest_links = np.array(closest_links)
        if len(closest_links) > 0:
            linewidths = 1.2*np.ones_like(closest_links)
            closestest_link = np.argmin(np.abs(closest_links - locus_frac))
            linewidths[closestest_link] = 3.5
        for k, linkage in enumerate(closest_links):
            plt.plot([x - width_to_chr_center, x - width_to_chr_center,
                      x + width_to_chr_center, x + width_to_chr_center],
                     [min_y + (max_y - min_y)*locus_frac, min_y + (max_y - min_y)*linkage,
                      min_y + (max_y - min_y)*linkage, min_y + (max_y - min_y)*locus_frac],
                     color=(1, 1, 1), transform=ax.transAxes,
                     linewidth=linewidths[k], linestyle='--', solid_capstyle='butt', zorder=100)
        all_closest_links.append(closest_links)

    return ax, all_closest_links


#chrv_size_bp = 576874
#location_ura_bp = np.mean([116167, 116970])

#chrv_size_kuhn = 1165
#kuhn_length = 0.015  # um
#chrv_size_effective_um = chrv_size_kuhn*kuhn_length
#location_ura_effective_um = location_ura_bp*(chrv_size_effective_um/chrv_size_bp)

#nuc_radius_um = 1.3  # Average of het5 msd convex hull distribution
#sim_nuc_radius_um = 1
#sim_D = 0.02  # um^2/s


def model_mscd(t, linkages, label_loc=location_ura_effective_um, chr_size=chrv_size_effective_um,
               nuc_radius=sim_nuc_radius_um, b=kuhn_length, D=sim_D, num_modes=10000):
    r"""
    Calculate the MSCD for the model of linked chromosomes

    Parameters
    ----------
    t : float array
        Time in seconds
    linkages : float array
        List of the link positions between the homologous chromosomes
    label_loc : float
        Location of the fluorescent label along the chromosome (microns)
    chr_size : float
        Length of the chromosome (microns)
    nuc_radius : float
        Radius of the nucleus (microns)
    b : float
        Kuhn length (microns)
    D : float
        Diffusivity (microns ** 2 / second)
    num_modes : int
        Number of normal modes used in the calculation

    Returns
    -------
    mscd_model : float array (size len(t))
        Calculated MSCD (microns ** 2) for the model with defined linkages

    """

    linkages = np.array(linkages) / b
    chr_size = chr_size / b
    label_loc = label_loc / b

    # Evaluate the MSCD if there are no linkages
    if len(linkages) == 0:
        mscd_model = 2 * linear_mid_msd(t, b, chr_size, D, num_modes)
        for i in range(len(t)):
            if mscd_model[i] > nuc_radius ** 2:
                mscd_model[i] = nuc_radius ** 2
        return mscd_model

    # Evaluate the MSCD if there are linkages between the chromosomes
    i = np.searchsorted(linkages, label_loc)
    if i == 0:
        mscd_func = linear_mscd
        Ndel = linkages[0] - label_loc
        N = 2 * linkages[0]
    elif i == len(linkages):
        mscd_func = linear_mscd
        Ndel = label_loc - linkages[-1]
        N = 2 * (chr_size - linkages[-1])
    else:
        mscd_func = ring_mscd
        Ndel = linkages[i] - label_loc
        N = 2 * (linkages[i] - linkages[i - 1])

    mscd_model = mscd_func(t, D=D, Ndel=Ndel, N=N, b=b, num_modes=num_modes)
    mscd_model = np.minimum(mscd_model, (nuc_radius ** 2) * np.ones(len(mscd_model)))

    return mscd_model


def model_mscd_confine(t, linkages, label_loc=location_ura_effective_um, chr_size=chrv_size_effective_um,
               nuc_radius=sim_nuc_radius_um, k_conf = 1, b=kuhn_length, D=sim_D, num_modes=10000):
    r"""
    Calculate the MSCD for the model of linked chromosomes

    Parameters
    ----------
    t : float array
        Time in seconds
    linkages : float array
        List of the link positions between the homologous chromosomes
    label_loc : float
        Location of the fluorescent label along the chromosome (microns)
    chr_size : float
        Length of the chromosome (microns)
    nuc_radius : float
        Radius of the nucleus (microns)
    k_conf : float
        Strength of confining potential
    b : float
        Kuhn length (microns)
    D : float
        Diffusivity (microns ** 2 / second)
    num_modes : int
        Number of normal modes used in the calculation

    Returns
    -------
    mscd_model : float array (size len(t))
        Calculated MSCD (microns ** 2) for the model with defined linkages

    """

    linkages = np.array(linkages) / b
    chr_size = chr_size / b
    label_loc = label_loc / b

    # Evaluate the MSCD if there are no linkages
    if len(linkages) == 0:
        mscd_model = 2 * linear_msd_confine(t, D, label_loc, chr_size, b, k_conf, num_modes)
        return mscd_model

    # Evaluate the MSCD if there are linkages between the chromosomes
    i = np.searchsorted(linkages, label_loc)
    if i == 0:
        mscd_func = linear_mscd
        Ndel = linkages[0] - label_loc
        N = 2 * linkages[0]
        mscd_plateau = 4 * b ** 2 * Ndel
    elif i == len(linkages):
        mscd_func = linear_mscd
        Ndel = label_loc - linkages[-1]
        N = 2 * (chr_size - linkages[-1])
        mscd_plateau = 4 * b ** 2 * Ndel
    else:
        mscd_func = ring_mscd
        Ndel = linkages[i] - label_loc
        N = 2 * (linkages[i] - linkages[i - 1])
        mscd_plateau = 2 * b ** 2 / (1 / (2 * Ndel) + 1 / (N - 2 * Ndel))

    if mscd_plateau < (nuc_radius ** 2):
        mscd_model = mscd_func(t, D, Ndel, N, b, num_modes)
    else:
        mscd_model = 2 * linear_msd_confine(t, D, label_loc, chr_size, b, k_conf, num_modes)

    return mscd_model


def generate_example_cell(mu, chr_size=chrv_size_effective_um):
    r"""
    Generate the number and location of linkage between homologous chromosomes

    Parameters
    ----------
    mu : float
        Average number of linkages between chromosomes (Poisson distributed)
    chr_size : float
        Size of the chromosome (microns)

    Returns
    -------
    cell : float array (length selected from Poisson distribution)
        List of linkage locations between the homologous chromosomes

    """
    cell = np.sort(chr_size * np.random.random_sample(size=np.random.poisson(lam=mu)))

    return cell


def model_plateau(linkages, label_loc=location_ura_effective_um, chr_size=chrv_size_effective_um,
                  nuc_radius=sim_nuc_radius_um, b=kuhn_length):
    r"""
    Evaluate the plateau values in the MSCD

    Parameters
    ----------
    linkages : float array
        List of the link positions between the homologous chromosomes
    label_loc : float
        Location of the fluorescent label along the chromosome (microns)
    chr_size : float
        Length of the chromosome (microns)
    nuc_radius : float
        Radius of the nucleus (microns)
    b : float
        Kuhn length (microns)

    Returns
    -------
    mscd_plateau : float
        Plateau value of the MSCD in the long-time asymptotic limit (microns ** 2)

    """

    chr_size /= b
    label_loc /= b
    linkages = np.array(linkages) / b

    # Evaluate the MSCD if there are no linkages
    if len(linkages) == 0:
        mscd_plateau = nuc_radius ** 2
        return mscd_plateau

    # Evaluate the MSCD if there are linkages between the chromosomes
    i = np.searchsorted(linkages, label_loc)
    if i == 0:
        Ndel = linkages[0] - label_loc
        mscd_plateau = 4 * b ** 2 * Ndel
    elif i == len(linkages):
        Ndel = label_loc - linkages[-1]
        mscd_plateau = 4 * b ** 2 * Ndel
    else:
        Ndel = linkages[i] - label_loc
        N = 2*(linkages[i] - linkages[i - 1])
        mscd_plateau = 2 * b ** 2 / (1 / (2 * Ndel) + 1 / (N - 2 * Ndel))

    if mscd_plateau > (nuc_radius ** 2):
        mscd_plateau = nuc_radius ** 2

    return mscd_plateau
