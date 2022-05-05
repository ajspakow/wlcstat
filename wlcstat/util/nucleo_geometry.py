# -*- coding: utf-8 -*-
r"""Geometry of a chromatin chain: a chain of nucleosomes

This module allows calculating various configurational properties of chains of
nucleosomes that rely solely on their idealized geometries (i.e. don't require
accounting for thermal fluctuations).

Example
-------
We can load in coordinates manually extracted from the nucleosome structure
(doi:10.1038/38444) PDB file, and fit a helix to them to extract the effective
helix parameters that describe the path DNA takes when wrapping the
nucleosome

.. plot::
    :context: reset
    :nofigs:

    >>> import nuc_chain.data as ncd
    >>> import nuc_chain.geometry as ncg
    >>> spokes = ncd.spokes_manual_deepti_2018_05_31
    >>> p_fit, opt = ncg.fit_spokes_helix(spokes)

We can convert the extracted helix parameters into more reasonable units and
center a 147bp-wound nucleosome

.. plot::
    :context:
    :nofigs:

    >>> p_nm = ncg.helix_params_angstrom_to_nm(p_fit)
    >>> p_nm = ncg.helix_params_to_centered_params(p_nm)

We can make a figure demonstrating the tangent, normals, and binormals of this
fit helix (see :py:func:`fit_spokes_helix` documentation for example of comparing the
fit to data

.. plot::
    :context: close-figs

    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> bp_i = np.arange(147)
    >>> H = ncg.H(bp_i, **p_nm)
    >>> U = ncg.Hu(bp_i, **p_nm)
    >>> N = ncg.Hn(bp_i, **p_nm)
    >>> B = ncg.Hb(bp_i, **p_nm)
    >>> ax.quiver3D(H[0], H[1], H[2], U[0], U[1], U[2], color=sns.color_palette()[0])
    >>> ax.quiver3D(H[0], H[1], H[2], N[0], N[1], N[2], color=sns.color_palette()[1])
    >>> ax.quiver3D(H[0], H[1], H[2], B[0], B[1], B[2], color=sns.color_palette()[2])

We can use the extracted helix parameters to plot a nucleosome chain with given
linker lengths and "unwrapping" states

.. plot::
    :context: close-figs

    >>> linker_lengths = np.array([35, 35, 35, 35, 35,
            45, 55, 60, 150, 35, 35, 35, 35, 35, 35, 35])
    >>> chain = ncg.minimum_energy_no_sterics(links=linker_lengths)
    >>> ncg.visualize_chain(*chain)

Coordinate Systems
------------------

"Standard Coordinates"
======================

The main coordinate system that will be used when possible is the coordinates
of the nucleosome itself. First we orient the best fit helix to the z-axis.
That is we fit the centerline of the DNA in the crystal structure of the
nucleosome to a helix defined by the equation (where i denotes the ith base
pair)

.. math::

    \begin{bmatrix}
    r*\cos(-2*\pi i\frac{T}{b}) \\
    r*\sin(-2*\pi i\frac{T}{b}) \\
    c*i
    \end{bmatrix}

Notice that this corresponds to a *left-handed* helix whose centerline is aligned with the
z-axis; whose base is on the x-y plane; and which "begins" at x=r, y=0.
:math:`b` is the number of bases (146 in the structure we
use), :math:`T` is the number of twists the helix undergoes, and :math:`c` the
total height of the helix. So :math:`c/T` is the pitch of the helix.

..todo::

    The following needs to be checked for off-by-one style error.

    The nucleosome crystal structure, as we've fit it, consists of 72 base pairs
    attached to one histone tetramer, a center basepair aligned with the dyad axis,
    and 73 base pairs attached to the other histone tetramer, in that order. We
    fit leaving 10 base pairs off on either side. So that leaves the fit
    containing 62/63 base pairs on the histones, resp.

    This means that the center of the helix is going to be at
    x=y=0,z=(62/125)*c. (62 out of 125 inter-nucleotide height changes,
    fenceposted between the 126 nucleotides we fit).

We will keep track of everything in terms of the entry orientation of the DNA into
the nucleosome core particle, and project all nucleosome-related vectors into
this basis in order to apply them in our "lab" frame.

Also notice that this means of the 10 parameters in our arbitrary helix fit,
only three of them (r, c, T, [b fixed]) determine the actual nucleosome
structure.

"Entry Coordinates"
===================

For actual propagation from nucleosome to nucleosome---and for plotting---we
will instead work in the rotated coordinate system defined by the three unit
material normals of the DNA parameterization.

The entry coordinate triad of material normals will be represented as the
rotation matrix [N, B, U] (alternately called t1, t2, t3).
Where by [N, B, U], we mean the triplet of (normalized) vectors corresponding to the normal,
binormal and tangent vector to the curve, respectively.

We can then simply use the [N, B, U] triplet at the terminal point on the helix
in these coordinates to define the next orientation. Conveniently, since the
third material normal always corresponds to the tangent vector of the curve, we
can always apply extra twist (or remove twist) *a posteriori* by rotating the
first two material normals about the third one.

.. caution:: This set of vectors is
    not torsion-free in general, and is in fact not for the case of a helix. So
    care must be taken to subtract out intrinsic twist before adding in known twist
    densities.

Units
-----
The fitting will be in ångströms, but then we will convert (and save) the
fitted values to nm and use those values henceforth. Parameterizations will be
in units of bp where possible.

A "linker length" is the number of base pairs starting not including the last
base pair bound to the previous nucleosome or the first base pair bound to the
next nucleosome. So the amount of DNA between nucleosomes if linker length is N
is (N+1)*length_per_bp

Notes
-----
This module is in active research-level development, it will be in constant
flux and not API is guaranteed to continue working. Consult Bruno Beltran
<brunobeltan0@gmail.com> and/or Deepti Kannan <dkannan@stanford.edu> if
attempting to integrate with this code base.

Some functions flat out don't work yet.

Well-tested functions include

- H, Hn, Hb, Hu
- H_oriented, Hn_oriented, Hb_oriented, Hu_oriented
- H_rot, H_rot_oriented
- all the fitting-related functions (*spokes*), H_arbitrary, fit_helix, etc.

Functions known to *not* work include

- minimum_energy_no_sterics, minimum_energy_no_sterics_const_linker
- visualize_chain

Tests
-----

I will organize these better at some point.

We have checked that H_oriented and H_rot_oriented are in fact equivalent to H
and H_rot (respectively) if the correct initial orientation is chosen::

    >>> import nuc_chain.geometry as ncg
    >>> h1 = ncg.H(i=None, **ncg.hp)
    >>> h2 = ncg.H_oriented(i=None, **ncg.hp, entry_rot=ncg.H_rot(0, **ncg.hp), unwrap=0)
    >>> assert(np.sum(h1 - h2) < 10e-15)
    >>> n1 = ncg.Hn(i=None, **ncg.hp)
    >>> n2 = ncg.Hn_oriented(i=None, **ncg.hp, entry_rot=ncg.H_rot(0, **ncg.hp), unwrap=0)
    >>> assert(np.sum(n1 - n2) < 10e-15)
    >>> b1 = ncg.Hb(i=None, **ncg.hp)
    >>> b2 = ncg.Hb_oriented(i=None, **ncg.hp, entry_rot=ncg.H_rot(0, **ncg.hp), unwrap=0)
    >>> assert(np.sum(b1 - b2) < 10e-15)
    >>> u1 = ncg.Hu(i=None, **ncg.hp)
    >>> u2 = ncg.Hu_oriented(i=None, **ncg.hp, entry_rot=ncg.H_rot(0, **ncg.hp), unwrap=0)
    >>> assert(np.sum(u1 - u2) < 10e-15)

We also checked that the unwrapped parameter in the *_oriented functions only
has an affect when we don't explictly ask for a particular i.

    >>> rot1 = ncg.H_rot(146, **ncg.hp)
    >>> rot2 = ncg.H_rot_oriented(146, **ncg.hp, entry_rot=ncg.H_rot(0, **ncg.hp), unwrap=10)
    >>> assert(np.sum(rot1 - rot2) < 10e-15)

    >>> h1 = ncg.H_oriented(146, **ncg.hp, entry_rot=np.identity(3), unwrap=0)
    >>> h2 = ncg.H_oriented(146, **ncg.hp, entry_rot=np.identity(3), unwrap=100)
    >>> assert(np.sum(h1 - h2) < 10e-15)

"""
import os
import sys
import inspect
from pathlib import Path
from string import Template
# from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool

import scipy
from scipy import stats
import numpy as np
import pandas as pd

from . import *
from . import data as ncd
from . import rotations as ncr
from . import wignerD as wd
from .rotations import Rz, Ry
from . import math as nuc_math
from . import linkers as ncl
from .linkers import convert

helix_params_doc = """r : float
        radius of helix
    c : float
        height of helix
    T : float
        number of turns
    b : float
        parametric unit in 1/T units. i.e. number of basepairs"""
helix_doc = """i : (N,) array_like
        parametric coordinate along helix. i.e. base pair number
    """ + helix_params_doc
helix_params_doc_extra = """x0 : float
        x-coordinate of center of base of cylinder,
    y0 : float
        y-coordinate of center of base of cylinder,
    z0 : float
        z-coordinate of center of base of cylinder,
    psi0 : float
        initial phase
    phi : float
        azimuthal angle of helical axis
    theta : float
        polar angle of helical axis"""
entry_rot_doc = """entry_rot : (3,3) array_like
        rotation matrix associated with the orientation of the entry site into
        the nucleosome"""
orientation_doc = entry_rot_doc + """
    entry_pos : (3,) array_like
        vector in lab frame pointing to starting location of helix, i.e. the
        "entry site" of the DNA into the nucleosome
    w_in : (optional) float
        number of nucleotides bound on entry side of central dyadic nucleotide.
    w_out : (optional) float
        number of nucleotides bound on exit side of central dyadic nucleotide."""
normed_doc = """normed : bool
        Whether to return normalized vectors (default True)"""

def H_arbitrary(i, r, c, T, b, x0, y0, z0, psi0, phi, theta):
    """Parametric description of a simple helix.

    Parameters
    ----------
    ${helix_doc}
    ${helix_params_doc_extra}

    Returns
    -------
    H : (N, 3) np.ndarray
        Position of helix at coordinate i.

    Notes
    -----
    The equation here has a non-standard convention (negative phase) to make
    the helix have the correct handedness.
    """
    ct = np.cos(theta)
    st = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)
    R = np.array([[1 - st**2*(1-ct)*cp**2, -st**2*(1-ct)*sp*cp,    st**2*cp],
                  [-st**2*(1-ct)*sp*cp,    1 - st**2*(1-ct)*sp**2, st**2*sp],
                  [-st**2*cp,                -st**2*sp,            1 - st**2*(1-ct)]])
    phase = (2*np.pi*i*T)/b
    straight_helix = np.array([r*np.cos(psi0 - phase), r*np.sin(psi0 - phase), c*i/b])
    return (np.dot(R, straight_helix).T + np.array([x0, y0, z0]))

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

def Hb_oriented(*args, **kwargs):
    """Bi-normal to :py:func:`H_oriented` at entry_rot.

    Parameters
    ----------
    ${helix_params_doc}
    ${entry_rot_doc}

    Returns
    -------
    b : (3,N) np.ndarray
        bi-normal vector for each i.
    """
    return HX_at_pos_rot_(Hb, *args, **kwargs)

def HX_at_pos_rot_(HX, i, r, c, T, b, entry_rot=None, entry_pos=None,
                 Lw=default_Lw):
    """General re-orienter used to implement the H*_oriented family of
    functions."""
    if i is None:
        i = np.arange(Lw + 1)
    # default to no rotation, in a kinda roundabout way
    cent_entry_rot = H_rot(0, r, c, T, b)
    cent_entry_rot = cent_entry_rot/np.linalg.norm(cent_entry_rot)
    if entry_rot is None:
        entry_rot = cent_entry_rot
    entry_rot = entry_rot/np.linalg.norm(entry_rot)
    transform_rot = entry_rot@np.linalg.solve(cent_entry_rot, np.identity(3))
    # default to no translation
    cent_entry_pos = transform_rot@HX(0, r, c, T, b)
    if entry_pos is None:
        entry_pos = cent_entry_pos
    entry_pos = np.reshape(entry_pos, (3,1))
    transform_pos = entry_pos - cent_entry_pos
    if len(transform_pos.shape) == 1:
        transform_pos = transform_pos[:,None]
    X = transform_rot@HX(i, r, c, T, b)
    return X + transform_pos

H_arbitrary.__doc__ = Template(H_arbitrary.__doc__).safe_substitute(
        helix_doc=helix_doc, helix_params_doc_extra=helix_params_doc_extra)
for h_func in [H]:
    h_func.__doc__ = Template(h_func.__doc__).safe_substitute(helix_doc=helix_doc)
for h_func in [Hu, Hn, Hb]:
    h_func.__doc__ = Template(h_func.__doc__).safe_substitute(
            helix_doc=helix_doc, normed_doc=normed_doc)
for h_func in [Hk, Htau]:
    h_func.__doc__ = Template(h_func.__doc__).safe_substitute(
            helix_params_doc=helix_params_doc)
for h_func in [H_oriented]:
    h_func.__doc__ = Template(h_func.__doc__).safe_substitute(
            helix_doc=helix_doc, orientation_doc=orientation_doc)
for h_func in [Hn_oriented, Hb_oriented, Hu_oriented]:
    h_func.__doc__ = Template(h_func.__doc__).safe_substitute(
            helix_params_doc=helix_params_doc, entry_rot_doc=entry_rot_doc)

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


def fit_helix(i, X, p0, bounds, p_fixed=None):
    """Fit a simple helix to data by minimizing l2 energy.

    Written as a helper function for fit_spokes_helix.

    Parameters
    ----------
    i : (N,) array_like
        Parametric position along the helix of each data point. See :py:func:`H`.
    X : (N,3) array_like
        Positions corresponding to positions in `i`.
    p0 : Dict[str, float]
        Initial guess for parameters to pass to energy minimization.
    bounds : (2,M)
        Min/max pairs for each parameter being varied. M == len(p0)
    p_fixed : Dict[str, float]
        Parameter names with values to leave them fixed to.

    Returns
    -------
    params : Dict[str, float]
        Best fit parameters for helix function H
    opt : OptimizeResult
        Result of scipy.optimize.least_squares
    """
    if p_fixed is None:
        p_fixed = []
    params_to_fit = [p for p in helix_param_order if p not in p_fixed]
    p0 = np.array([p0[h] for h in params_to_fit])
    bounds = np.array([[bounds[p][0] for p in params_to_fit],
                       [bounds[p][1] for p in params_to_fit]])
    # fit is by minimization of energy function
    def energy(params):
        # back out what parameters were which by closure of params_to_fit
        p = {p: params[i] for i, p in enumerate(params_to_fit)}
        return np.linalg.norm(np.linalg.norm(
                H_arbitrary(i, **p, **p_fixed) - X, axis=1))
    # remember to scale jacobian for angles vs lengths with xscale
    scale = {'r': 3*np.sqrt(2), 'x0': 1, 'y0': 1, 'z0': 1, 'psi0': np.pi/2000,
             'T': 0.01, 'c': 1, 'theta': np.pi/2000, 'phi': np.pi/2000}
    scale = [scale[p] for p in params_to_fit]
    opt = scipy.optimize.least_squares(energy, x0=p0,
                x_scale=scale, jac='2-point', bounds=bounds,
                method='trf', max_nfev=None)
    # back out the same as in energy function
    p_fit = {p: opt.x[i] for i, p in enumerate(params_to_fit)}
    return p_fit, opt

def fit_spokes_helix(basepair_id=None, spokes=None, p0=None, bounds=None,
                     uw_in=default_uw_in, uw_out=default_uw_out):
    """Fit a simple helix to the terminal carbons on each base pair in the
    nucleosome structure (doi:10.1038/38444).

    Allows specifying how many base pairs on each side should count as being
    "unwrapped" in the crystal structure, so as to exclude them from the fit.

    Parameters
    ----------
    basepair_id : (292,) array_like, optional
        The "i"'s corresponding to the bases in spokes. Defaults to in-order,
        i.e. [1,1,2,2,..,146,146].
        Ignored if spokes is not specified.
    spokes : (292, 3) array_like, optional
        The carbon locations extracted from the crystal structure. The default
        helix fit parameters assume units of ångströms.
        Defaults to :py:mod:`nuc_chain.data.spokes_manual_deepti_2018_05_31`.
        If default used, basepair_id is ignored.
    p0 : Dict[str, float], optional
        The initial guesses for the parameters to the helix function H. Order
        is not imporant. b should not be included.
        Defaults to values chosen by Deepti, saved in :py:mod:`nuc_chain.data`.
    bounds : (2, N) array_like, optional
        Min and max values to search for the helix parameters. N = number of
        helix parameters minux 'b'.
        Defaults to values chosen by Deepti, saved in :py:mod:`nuc_chain.data`.
    uw_in : int, optional
        How many bases are unwrapped from nucleosome crystal structure on
        "entry" side (arbitrarily chosen?).
        Defaults to 10.
    uw_out : int, optional
        How many bases are unwrapped from nucleosome crystal structure on
        "exit" side (arbitrarily chosen?)
        Defaults to 10.

    Returns
    -------
    params : Dict[str, float]
        Best fit parameters for helix function H
    opt : OptimizeResult
        Result of scipy.optimize.least_squares

    Example
    -------
    Performing the fit, and checking it against data visually::
        >>> import matplotlib as mpl
        >>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits.mplot3d import import Axes3D
        >>> import nuc_chain.geometry as ncg
        >>> import nuc_chain.data as ncd
        >>> # do fit
        >>> p_fit, opt = ncg.fit_spokes_helix()
        >>> print(p_fit)
            {'T': 1.7000000000000002,
            'b': 146,
            'c': 44.297700263176672,
            'phi': 1.6771849694713676,
            'psi0': 1.408727503155174,
            'r': 41.192061481066709,
            'theta': -1.7202108711768689,
            'x0': 10.47376213526613,
            'y0': -44.531679077103604,
            'z0': -47.455752173714558}
        >>> # load data to compare to
        >>> spokes = ncd.spokes_manual_deepti_2018_05_31
        >>> spokes_fit = ncg.H_arbitrary(spokes['basepair id'], **p_fit)
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> ax.scatter(*(spokes_fit[:,i] for i in [0,1,2]))
        >>> ax.scatter(*(spokes[x] for x in ['x', 'y', 'z']))

    Notes
    -----
    In (Luger et al., Nature 1997), the structure of the DNA is used in its
    entirety to estimate a superhelical radius (through base pairs) of 41.8Å
    with a pitch (height of one turn) of 23.9Å with 1.65 superhelical turns.

    Our fit above estimates a radius of 41.9Å and a pitch of ~26Å, with 1.7
    superhelical turns. So our fit estimates more turns and larger pitch...but
    we agree on the radius!
    """
    #fix b = 146 bp (number of basepairs in pdb nucleosome crystal structure)
    b = bp_in_structure
    if basepair_id is None:
        basepair_id = np.arange(0, b)
    if spokes is None:
        spokes = ncd.spokes_manual_deepti_2018_05_31
        basepair_id = spokes['basepair id'].values
        spokes = spokes[['x','y','z']].values
    if p0 is None:
        r0 = np.array([10, -45, -52])
        p0 = {'r': 43.0,
              'r0': np.array([10., -45., -52.]),
              'x0': r0[0],
              'y0': r0[1],
              'z0': r0[2],
              'psi0': 0.75*np.pi/2.0,
              'T': 1.8,
              'c': 42,
              'phi': (1.1*np.pi/2.0),
              'theta': (-1.1*np.pi / 2.0)
        }
    if bounds is None:
        bounds = {'r': (35, 50), 'c': (35, 50), 'T': (1.7, 2.1),
                'x0': (0, 20), 'y0': (-60, 40), 'z0': (-100, -40),
                'psi0': (np.pi/4, np.pi/2), 'phi': (np.pi/2, np.pi),
                'theta': (-3*np.pi/4, -2*np.pi/5)}
    # only fit to non-unwrapped bases
    ok_bp = (basepair_id >= uw_in) & (basepair_id < b - uw_out)
    basepair_id = basepair_id[ok_bp]
    spokes = spokes[ok_bp, :]
    # make sure to not count the base pairs we subtracted off!
    p_fit, opt = fit_helix(basepair_id, spokes, p0, bounds, p_fixed={'b': b-uw_in-uw_out})
    p_fit['b'] = b - uw_in - uw_out
    return p_fit, opt

def spoke_vectors(spokes=None):
    """Get the vector joining each pair of spokes (basepairs)"""
    if spokes is None:
        spokes = ncd.spokes_manual_deepti_2018_05_31
    spokes_watson = spokes[spokes['base id'] == 0]
    spokes_crick = spokes[spokes['base id'] == 1]
    for df in [spokes_watson, spokes_crick]:
        df.set_index('basepair id', inplace=True)
    base_pair_vectors = spokes_watson - spokes_crick
    return base_pair_vectors[['x', 'y', 'z']]

def spoke_cotangent_coordinates(spokes=None):
    """Get projection of the vector joining each pair of spokes (basepairs)
    into the cotangent plane of the helix at that location."""
    def vectors_to_proj(vec):
        "series with index ['x','y','z'] and name==i -> cotangent coordinates"
        return pd.Series({
            'x': np.dot(vec.values, Hn(vec.name, **helix_params_best)[:,0]),
            'y': np.dot(vec.values, Hb(vec.name, **helix_params_best)[:,0])
        })
    return spoke_vectors(spokes).apply(vectors_to_proj, axis=1)

def cotangent_angles(spokes=None):
    cotangent_proj = spoke_cotangent_coordinates(spokes)
    theta = np.arctan2(cotangent_proj['x'], cotangent_proj['y'])
    return theta

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

    :math:`\Omega_C \cdot R_z(L_w\left[\tau_n - \tau_H])`

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
    theta = Lw*(tau_n - Htau(**helix_params))
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

def entry_exit_vector(Lw, helix_params=helix_params_best):
    """Vector pointing from entry to exit of nucleosome (including correct
    magnitude) assuming that the entry orientation is such that the [N,B,U]
    basis aligns with the x,y,z axis.

    Multiply this on the left with your actual entry orientation matrix to get
    the vector that points from the entry location to the exit location."""
    # minus 1's since e.g. 147th nucleotide is at index 146
    v = H(Lw, **helix_params) - H(0, **helix_params)
    return np.linalg.solve(H_rot(0, **helix_params), np.identity(3))@v

def dp_omega(Ll, *, w_in=default_w_in, w_out=default_w_out,
             tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
             lpb=dna_params['lpb'], helix_params=helix_params_best,
             unwrap=None):
    """The displacement in nm and rotation between nucleosome entry sites from a
    single linker (or for the case of constant linker and unwrapping lengths).
    Here, "entry site" means the first bound base pair, so this function takes
    you from the first bound base pair, around the nucleosome, down the
    unwrapped base pairs on the outgoing side, down the linker, and down the
    unwrapped base pairs on the incoming side of the next nucleosome.

    Returns
    -------
    dp : (3,) np.ndarray[float64]
        Displacement from entry site to entry site.
    OmegaNext : (3,3) np.ndarray[float64]
        Rotation matrix from entry site to entry site.
    """
    if unwrap is not None:
        w_in, w_out = convert.resolve_unwrap(unwrap)
    b = helix_params['b']
    Lw = w_in + w_out
    unwrap = b - Lw - 1
    Onext = OmegaNextEntry(Ll, tau_n=tau_n, tau_d=tau_d,
            w_in=w_in, w_out=w_out, helix_params=helix_params)
    p0 = np.zeros((3,))
    Omega0 = np.identity(3)
    Omega1 = Omega0 @ Onext
    exit_pos = (Omega0@entry_exit_vector(Lw)).T[0] + p0
    exit_u = Omega1[:,2]
    exit_u = exit_u/np.linalg.norm(exit_u)
    p1 = exit_pos + exit_u*(Ll + unwrap)*lpb
    # we can save multiplying by Omega0^{-1} since it's the Id(3)
    return p1 - p0, Omega1

def dp_omega_linker_only(Ll, *, w_in=default_w_in, w_out=default_w_out,
             tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
             lpb=dna_params['lpb'], helix_params=helix_params_best,
             unwrap=None):
    """Fixing nucleosome entry and exit to be in the same place, the
    displacement and rotation between nucleosome entry sites from a single
    linker for the case of constant linker and unwrapping lengths.

    this function is useful to compare to the theory (where one linker starts
    exactly where the other ends.

    Returns
    -------
    dp : (3,) np.ndarray[float64]
        Displacement from entry site to entry site.
    OmegaNext : (3,3) np.ndarray[float64]
        Rotation matrix from entry site to entry site.
    """
    if unwrap is not None:
        w_in, w_out = convert.resolve_unwrap(unwrap)
    b = helix_params['b']
    Lw = w_in + w_out
    unwrap = b - Lw - 1
    Onext = OmegaNextEntry(Ll, tau_n=tau_n, tau_d=tau_d,
            w_in=w_in, w_out=w_out, helix_params=helix_params)
    p0 = np.zeros((3,))
    Omega0 = np.identity(3)
    Omega1 = Omega0 @ Onext
    exit_u = Omega1[:,2]
    exit_u = exit_u/np.linalg.norm(exit_u)
    p1 = p0 + exit_u*(Ll + unwrap)*lpb
    # we can save multiplying by Omega0^{-1} since it's the Id(3)
    return p1 - p0, Omega1

def dp_omega_exit(Ll, *, w_in=default_w_in, w_out=default_w_out,
             tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
             lpb=dna_params['lpb'], helix_params=helix_params_best,
             unwrap=None):
    """The displacement in nm and rotation between nucleosome exit sites from a
    single linker (or for the case of constant linker and unwrapping lengths).
    Here, "exit site" means the last bound base pair, so this function takes
    you from the last bound base pair, down the unwrapped base pairs from the
    outgoing side of the previous nucleosome, down the linker, down the
    unwrapped base pairs on the incoming side of the nucleosome, and
    finally around that nucleosome.

    Returns
    -------
    dp : (3,) np.ndarray[float64]
        Displacement from entry site to entry site.
    OmegaNext : (3,3) np.ndarray[float64]
        Rotation matrix from exit site to exit site.
    """
    if unwrap is not None:
        w_in, w_out = convert.resolve_unwrap(unwrap)
    b = helix_params['b']
    Lw = w_in + w_out
    unwrap = b - Lw - 1
    Onext = OmegaNextExit(Ll, tau_n=tau_n, tau_d=tau_d,
            w_in=w_in, w_out=w_out, helix_params=helix_params)
    p0 = np.zeros((3,)) # previous exit site
    Omega0 = np.identity(3) # previous exit orientation
    # we start off in the plus z direction (np.identity orientation)
    entry_pos = p0 + lpb*np.array([0, 0, Ll + unwrap])
    p1 = entry_pos + (Omega0@entry_exit_vector(Lw)).T[0]
    # we can save multiplying by Omega0^{-1} since it's the Id(3)
    return p1 - p0, Onext


def dp_omega_nuc(Ll, *, w_in=default_w_in, w_out=default_w_out,
             tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'],
             lpb=dna_params['lpb'], helix_params=helix_params_best,
             unwrap=None):
    """The displacement and rotation between nucleosome centers from a single
    linker for the case of constant linker and unwrapping lengths.

    ..todo::
        Test this function

    Returns
    -------
    dp : (3,) np.ndarray[float64]
        Displacement from entry site to entry site.
    OmegaNext : (3,3) np.ndarray[float64]
        Rotation matrix from entry site to entry site.
    """
    if unwrap is not None:
        w_in, w_out = convert.resolve_unwrap(unwrap)
    N_rot, N_pos = entry_to_nuc_center(w_in, helix_params=helix_params)
    b = helix_params['b']
    Lw = w_in + w_out
    unwrap = b - Lw - 1
    Onext = OmegaNextEntry(Ll, tau_n=tau_n, tau_d=tau_d,
            w_in=w_in, w_out=w_out, helix_params=helix_params)
    p0 = np.zeros((3,))
    Omega0 = np.identity(3)
    Omega1 = Omega0 @ Onext
    exit_pos = (Omega0@entry_exit_vector(Lw)).T[0] + p0
    exit_u = Omega1[:,2]
    exit_u = exit_u/np.linalg.norm(exit_u)
    p1 = exit_pos + exit_u*(Ll + unwrap)*lpb
    # now use these entry pos/rots to get the nucleosome center coordinates
    n0 = p0 + Omega0@N_pos
    n1 = p1 + Omega1@N_pos
    N0_rot_inv = np.linalg.solve(Omega0@N_rot, np.identity(3))
    N1_rot = Omega1@N_rot
    return n1 - n0, N0_rot_inv @ N1_rot

def entry_to_nuc_center(w_in, helix_params=helix_params_best):
    """Returns rotation and offset taking entry site to nucleosome center.

    Nucleosome center orientation defined to be t1,t2,t3 pointing away toward
    the dyad nucleotide (in plane perpendicular to helical axis), pointing
    along the helical axis, and their cross product t3 = t1 x t2."""
    # use the easy formula for the center of H, then move into entry_rot
    # coordinates
    N_pos = np.array([0, 0, H(w_in, **helix_params)[2]]) - H(0, **helix_params)[:,0]
    # now project these onto the plane perpendiculuar to the helical axis
    hu = Hu(w_in, **helix_params)[:,0]
    hn = Hn(w_in, **helix_params)[:,0]
    z = np.array([0,0,1])
    Nu = hu - np.dot(hu, z)*z
    Nu = Nu/np.linalg.norm(Nu)
    Nn = hn - np.dot(hn, z)*z
    Nn = Nn/np.linalg.norm(Nn)
    Nb = np.cross(Nu, Nn)
    N_rot = np.stack([Nn, Nb, Nu], axis=1)
    H_rot_inv = np.linalg.solve(H_rot(0, **helix_params), np.identity(3))
    return H_rot_inv@N_rot, H_rot_inv@N_pos

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
    w_ins, w_outs = convert.resolve_wrapping_params(unwraps, w_ins, w_outs, num_nucleosomes)
    # initialize to valid orientation matrix
    entry_rots = np.tile(np.identity(3), (num_nucleosomes, 1, 1))
    # and to start at the origin
    entry_pos = np.zeros((num_nucleosomes, 3))
    for i in range(num_linkers):
        # w_ins[i+1] because the ith linker is between i, and i+1 nucs
        Onext = OmegaNextEntry(links[i], tau_n=tau_n, tau_d=tau_d,
                w_in=w_ins[i+1], w_out=w_outs[i], helix_params=helix_params)
        if random_phi is not None:
            Onext = ncr.Rz(random_phi*np.random.rand()) @ Onext
        entry_rots[i+1] = entry_rots[i] @ Onext
        exit_u = entry_rots[i+1,:,2]
        exit_u = exit_u/np.linalg.norm(exit_u)
        mu_out = (b - 1)/2 - w_outs[i]
        mu_in = (b - 1)/2 - w_ins[i+1]
        #no translation included in this treatment for comparison to theory
        entry_pos[i+1] = entry_pos[i] + exit_u*(mu_out + links[i] + mu_in)*lpb
    return entry_rots, entry_pos

def minimum_energy_no_sterics(links, *, w_ins=default_w_in,
        w_outs=default_w_out, tau_n=dna_params['tau_n'],
        tau_d=dna_params['tau_d'], lpb=dna_params['lpb'],
        helix_params=helix_params_best, unwraps=None):
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
    w_ins, w_outs = convert.resolve_wrapping_params(unwraps, w_ins, w_outs, num_nucleosomes)
    # initialize to valid orientation matrix
    entry_rots = np.tile(np.identity(3), (num_nucleosomes, 1, 1))
    # and to start at the origin
    entry_pos = np.zeros((num_nucleosomes, 3))
    for i in range(num_linkers):
        # w_ins[i+1] because the ith linker is between i, and i+1 nucs
        Onext = OmegaNextEntry(links[i], tau_n=tau_n, tau_d=tau_d,
                w_in=w_ins[i+1], w_out=w_outs[i], helix_params=helix_params)
        entry_rots[i+1] = entry_rots[i] @ Onext
        exit_pos = (entry_rots[i]@entry_exit_vector(w_ins[i]+w_outs[i])).T[0] \
                + entry_pos[i]
        exit_u = entry_rots[i+1,:,2]
        exit_u = exit_u/np.linalg.norm(exit_u)
        mu_out = (b - 1)/2 - w_outs[i]
        mu_in = (b - 1)/2 - w_ins[i+1]
        entry_pos[i+1] = exit_pos + exit_u*(mu_out + links[i] + mu_in)*lpb
    return entry_rots, entry_pos

def wignerD_weights(links, weights, nlam, w_ins = default_w_in, w_outs = default_w_out,
    tau_n=dna_params['tau_n'], tau_d=dna_params['tau_d'], unwraps=None):
    """Computes weighted average of wigner Ds corresponding to the OmegaNext rotation operations
    due to a set of linker lengths and their associated probabilities.

    Returns
    -------
    Dlist : list of wigner_d_vals, weighted average of rotations
    """
    N = links.size
    if weights.size != N:
        raise ValueError("links and weights must be the same dimension")
    Dlists = np.ones_like(links) #list of Dlists, each element returned by ncr.wignerD_from_R
    for i in range(N):
        Onext = OmegaNextEntry(links[i], tau_n=tau_n, tau_d=tau_d,
                w_in=w_ins[i], w_out=w_outs[i], helix_params=helix_params_best)
        Dlist = ncr.wignerD_from_R(Onext, nlam)
        for j in range(nlam):
            #flatten each [mp-l, m-l] matrix in Dlist into a 1d array (rowwise)
            Dlist[j] = Dlist[j].flatten()
        #now Dlist is just a list of numbers
        Dlists[i] = np.concatenate(Dlist)
    #test Dlist compression: passed!
    #do weighted average: returns a Dlist
    return np.average(Dlists, weights=weights, axis=0)

def wignerD_probability(Dlist, R, nlam):
    """Returns the probability of rotation R based on probability distribution given by wigner D.
    TODO: test!

    Parameters
    ----------
    Dlist : array-like
        wignerD probability distribution function (output of wignerD_weights(...))
    R : (3,3)
        Rotation matrix
    nlam: int
        number of l values used to create Dlist

    Returns
    -------
    prob : float
        probability of rotation R, hopefully a number between 0 and 1(?)

    Notes
    -----
    :math:`P(R) = \sum_{l,m,j}{\mathcal{D}}^{mj}_{l}(R)(\sum_{i}p_i{\mathcal{D}}^{mj}_{l}(\Omega_i)^*`

    """
    DfromR = ncr.wignerD_from_R(R, nlam)
    #compress DfromR into a list of numbers (to match output of wignerD_weights())
    for l in range(nlam):
        DfromR[l] = DfromR[l].flatten()
    DfromR = np.concatenate(DfromR)
    #See formula in Deepti's Wigner D latex file
    return np.dot(DfromR, np.conjugate(Dlist))


def tabulate_e2e_vectors(*, tau_n=dna_params['tau_n'], unwrap=None):
    """Return a lookup table of entry->exit vectors with the right magnitude in
    nm. Multiply on the left with the entry orientation matrix to obtain the
    entry to exit displacement vector.

    One vector for each possible level of unwrapping.

    Returns
    -------
    entry_to_exit_vectors: pd.DataFrame
        columns are x, y, z components of vector
    """
    if unwrap is None:
        unwrap = np.arange(0, 147)
    Lws = np.array([(bp_in_nuc - u - 1) for u in unwrap])
    #list of vectors, one for each unwrapping level
    e2evecs = np.zeros((Lws.size, 3))
    for i, lw in enumerate(Lws):
        e2evecs[i,:] = entry_exit_vector(lw).T[0]
    df = pd.DataFrame(e2evecs, columns=['x', 'y', 'z'])
    return df

def tabulate_e2e_rotations(*, tau_n=dna_params['tau_n'], unwrap=None):
    """Return a lookup table of entry->exit rotation matrices.

    One value for each possible level of unwrapping.

    Returns
    -------
    OmegaE2E matrices : pd.DataFrame
        single column of matrices
    """
    if unwrap is None:
        unwrap = np.arange(0, 147)
    Lws = np.array([(bp_in_nuc - u - 1) for u in unwrap])
    #list of rotations, one for each unwrapping level
    rots = []
    for lw in Lws:
        rots.append(OmegaE2E(lw, tau_n=tau_n))
    df = pd.DataFrame(np.array(rots).reshape(3*Lws.size, 3), columns=['Rx', 'Ry', 'Rz'])
    return df

def tabulate_exit_angles(*, tau_n=dna_params['tau_n'], unwrap=None):
    """Return a lookup table of entry->exit (zyz) euler angles.

    One value for each possible level of unwrapping.

    Returns
    -------
    angles : pd.DataFrame
        columns are unwrap, alpha, beta, gamma.
    """
    if unwrap is None:
        unwrap = np.arange(-50, 147)
    alpha = np.zeros_like(unwrap).astype(float)
    beta = np.zeros_like(unwrap).astype(float)
    gamma = np.zeros_like(unwrap).astype(float)
    for i,u in enumerate(unwrap):
        w_in, w_out = convert.resolve_unwrap(u)
        R = OmegaE2E(w_in+w_out, tau_n=tau_n)
        alpha[i], beta[i], gamma[i] = ncr.zyz_from_matrix(R)
        R_reconstruct = Rz(alpha[i])@Ry(beta[i])@Rz(gamma[i])
        if np.sum(R_reconstruct - R) > 10e-10:
            raise ValueError('Failed to reconstruct R from zyz.')
    df = pd.DataFrame(np.array([unwrap, alpha, beta, gamma]).T,
            columns=['unwrap', 'alpha', 'beta', 'gamma'])
    return df

def tabulate_rise(dp_f=dp_omega):
    links = np.arange(10, 200)
    unwraps = np.arange(147)
    num_links = len(links)
    num_wraps = len(unwraps)
    rise = np.zeros((num_links,num_wraps))
    angle = np.zeros((num_links,num_wraps))
    radius = np.zeros((num_links,num_wraps))
    link_ix = np.zeros((num_links,num_wraps))
    unwrap_ix = np.zeros((num_links,num_wraps))
    for i, link in enumerate(links):
        for j, u in enumerate(unwraps):
            link_ix[i,j] = link
            unwrap_ix[i,j] = u
            dP, Onext = dp_f(link, unwrap=u)
            rise[i,j], angle[i,j], radius[i,j] \
                    = ncr.R_dP_to_screw_coordinates(Onext, dP)
    # df = pd.DataFrame(np.array([link_ix.flatten(), unwrap_ix.flatten(),
    #                             rise.flatten(), angle.flatten(),
    #                             radius.flatten()]).T,
    #                   columns=['linker length', 'unwrap', 'rise', 'angle',
    #                            'radius'])
    return link_ix, unwrap_ix, rise, angle, radius
    # return df

def hindered_rotating_kuhn_multiplier_expon(mu, w_in=default_w_in, w_out=default_w_out):
    links = np.arange(int(mu*100))
    phis = np.array([ncr.phi_theta_alpha_from_R(ncg.OmegaNextExit(
            link, w_in, w_out))[0] for link in links])
    cos_phi = np.mean(np.cos(phis))
    theta = ncr.phi_theta_alpha_from_R(ncg.OmegaNextExit(mu, w_in, w_out))[1]
    cos_theta = np.cos(theta)
    return (1 + cos_theta)/(1 - cos_theta)*(1 + cos_phi)/(1 - cos_phi)

def free_rotating_r2(links, theta):
    i = np.arange(len(links)) + 1
    return np.sum(links[:,None]*links[None,:] \
                  *np.power(np.cos(theta), np.abs(i[:,None] - i[None,:])))

def phi_from_link(Ll, w_in, w_out):
    if (Ll, w_in, w_out) in phi_from_link.cache:
        return phi_from_link.cache[(Ll, w_in, w_out)]
    return ncr.phi_theta_alpha_from_R(OmegaNextExit(Ll, w_in=w_in, w_out=w_out))[0]
phi_from_link.cache = {}

def tabulate_r2_heterogenous_chains_by_variance(num_chains, chain_length, sigmas, mu=35, pool_size=None, random_phi=None, **kwargs):
    n_sig = len(sigmas)
    links = np.zeros((n_sig, num_chains, chain_length-1))
    for i in range(n_sig):
        links[i,:,:] = ncl.fake_linkers_increasing_variance(mu, sigmas[i], size=(num_chains,chain_length-1), type='box')
    chains = np.zeros((n_sig, num_chains, chain_length, 3))
    rmax = np.zeros((n_sig, num_chains, chain_length))
    r2 = rmax.copy()
    variance = rmax.copy()
    chain_id = rmax.copy()
    def given_ij(ij):
        i, j = ij
        chain = minimum_energy_no_sterics_linker_only(links[i,j,:].flatten(), random_phi=random_phi, **kwargs)
        chains[i,j,:,:] = chain[1]
        r2[i,j] = nuc_math.r2(chains[i,j,:,:], axis=0)
        #Rmax in nm
        rmax[i,j] = np.concatenate(([0], convert.Rmax_from_links_unwraps(links[i,j], **kwargs))) * dna_params['lpb']
        variance[i] = sigmas[i]
        chain_id[i,j] = j
    if pool_size is None:
        for i in range(n_sig):
            for j in range(num_chains):
                given_ij((i,j))
    else:
        with Pool(processes=pool_size) as p:
            p.map(given_ij, [(i,j) for i in range(n_sig) for j in range(num_chains)])
    df = pd.DataFrame(np.stack([
        r2.flatten(), rmax.flatten(), variance.flatten(), chain_id.flatten()
        ], axis=1), columns=['r2', 'rmax', 'variance', 'chain_id'])
    return df

def tabulate_r2_heterogenous_chains_exponential(num_chains, chain_length, mu=35, pool_size=None, **kwargs):
    links = ncl.independent_linker_lengths(mu, size=(num_chains,chain_length-1))
    chains = np.zeros((num_chains, chain_length, 3))
    rmax = np.zeros((num_chains, chain_length))
    r2 = rmax.copy()
    chain_id = rmax.copy()
    def given_chain_i(i):
        chain = minimum_energy_no_sterics_linker_only(links[i,:].flatten(), **kwargs)
        chains[i,:,:] = chain[1]
        r2[i] = nuc_math.r2(chains[i,:,:], axis=0)
        #Rmax in nm
        rmax[i] = np.concatenate(([0], convert.Rmax_from_links_unwraps(links[i], **kwargs))) * dna_params['lpb']
        chain_id[i] = i
    if pool_size is None:
        for i in range(num_chains):
                given_chain_i(i)
    else:
        raise NotImplementedError('No pools here plz.')
    #     with Pool(processes=pool_size) as p:
    #         p.map(given_ij, [(i,j) for i in range(n_sig) for j in range(num_chains)])
    df = pd.DataFrame(
        np.stack([ r2.flatten(), rmax.flatten(), chain_id.flatten()], axis=1),
        columns=['r2', 'rmax', 'chain_id'])
    df['mu'] = mu
    return df

