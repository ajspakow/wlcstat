.. _wlcgreen_notes:
.. .. automodule:: wlcstat


Wormlike Chain Green's Function
===============================

We consider the end-to-end distribution function :math:`G(\vec{R},\vec{u}|\vec{u}_{0};L)` (or Green's function),
which gives the probability that a chain that begins at the origin with initial
orientation :math:`\vec{u}_{0}` will have end position :math:`\vec{R}` and final orientation :math:`\vec{u}`.
The Green's function
:math:`G(\vec{R},\vec{u}|\vec{u}_{0};L)` is written as

.. math::
    G(\vec{R},\vec{u}|\vec{u}_{0};L) = \int_{\vec{u}(s=0)=\vec{u}_{0}}^{\vec{u}(s=L)=\vec{u}}
    \hspace{-0.3in} \mathcal{D}[\vec{u}(s)]
    \exp \left\{ -\beta \mathcal{H}_{0} [\vec{u}(s)] \right\}
    \delta \left( \vec{R}-\int_{0}^{L}ds \vec{u} \right),

.. \label{eq:ORgreenfunc}

where :math:`\delta` is a Dirac delta function that restricts the path integration
to those paths that satisfy the fixed end separation :math:`\vec{R}`.
The path integral formulation of the Green function is used to find the
governing "Schrödinger" equation [Yamakawa1997]_

.. math::
    \left(
    \frac{\partial}{\partial L} - \frac{1}{2l_{p}} \tilde{\nabla}^{2}_{D} + \vec{u} \cdot \vec{\nabla}_{R}
    \right)
    G(\vec{R},\vec{u}|\vec{u}_{0};L)=
    \delta \left( L \right) \delta ( \vec{R} ) \delta \left( \vec{u} - \vec{u}_{0} \right),

.. \label{eq:schrodinger-r}

where :math:`\tilde{\nabla}^{2}_{D}` is the :math:`D`-dimensional angular Laplacian,
and :math:`\vec{\nabla}_{R}` is the :math:`D`-dimensional gradient operator.

Upon :math:`D`-dimensional Fourier transforming
from the variable :math:`\vec{R}` to the wavevector :math:`\vec{k}`,
our problem becomes that of a single wormlike chain
with Hamiltonian

.. math::
    \beta \mathcal{H} = \beta \mathcal{H}_{0} + \beta \mathcal{H}_{ext} =
    \frac{l_{p}}{2} \int_{0}^{L} ds
    \left( \frac{\partial \vec{u}}{\partial s} \right)^{\! \! 2}
    -i\vec{k} \cdot \int_{0}^{L}ds \vec{u}.

.. \label{eq:TotalHamiltonian}

The external Hamiltonian :math:`\beta \mathcal{H}_{ext}`
naturally emerges
in the Fourier transform of the end-to-end distribution function :math:`G(\vec{R},\vec{u}|\vec{u}_{0};L)`
due to the end-position constraint.
The corresponding "Schrödinger" equation in Fourier space is

.. math::
    \left(
    \frac{\partial}{\partial L} - \frac{1}{2l_{p}} \tilde{\nabla}^{2}_{D} - i \vec{k}\cdot \vec{u}
    \right)
    G(\vec{k},\vec{u}|\vec{u}_{0};L)=
    \delta \left( L \right) \delta \left( \vec{u} - \vec{u}_{0} \right),

.. \label{eq:schrodinger-k}

which demonstrates the direct correspondence between our problem and that of a rigid rotor
in a dipole field.

Exact results for the Green's function are provided in Refs.
[Spakowitz2004]_, [Spakowitz2005]_, [Spakowitz2006]_,  and [Mehraeen2008]_,
which give solutions with and without end constraints and in arbitrary dimensions :math:`D`.
For example, the solution to the orientation-independent Green function in Fourier-Laplace space
(i.e. Fourier transformed from :math:`\vec{R}` to :math:`\vec{k}` and Laplace transformed
from :math:`N` to :math:`p`)

.. math::
    G(\vec{K};p) = \frac{1}{P_{0}+\frac{\left(a_{1} K\right)^2}{P_{1}+
    \frac{\left(a_{2} K\right)^2}{P_{2}+\frac{\left(a_{3} K\right)^2}{\cdots}}}},
    :label: gwlc

where :math:`\vec{K}=2l_{p} \vec{k}`,
:math:`K` is the magnitude of :math:`\vec{K}`,
:math:`P_{\lambda}=p+ \lambda (\lambda+D-2)`, and
:math:`a_{\lambda} = \left[\frac{\lambda (\lambda+D-3)}{(2\lambda+D-2)(2\lambda+D-4)}\right]^{1/2}`.
A detailed derivation of this :math:`D`-dimensional solution are found in Ref. [Mehraeen2008]_.
We note that this exact solution produces identical results as those given in
Ref. [Spakowitz2004]_
for 2- and 3-dimensional solutions.

Methods to Evaluate the Wormlike Chain Green's Function
-------------------------------------------------------

The Green's function (e.g. :eq:`gwlc`) is written in Fourier-Laplace space.
Thus, evaluation of the Green's function requires inversion to real space
:math:`\vec{R}` and :math:`N`.
The python module 'wlcgreen' contains python code to numerically invert the
Fourier and Laplace transforms to render real-space solutions.
The procedure followed in this module is as follows:

- Laplace-space inversion (from :math:`p` to :math:`N`): Numerical Laplace inversion using residue theorem based on evaluation of the singularities in :math:`G`

- Fourier-space inversion (from :math:`\vec{k}` to :math:`\vec{R}`): Numerical Fourier inversion by numerical integration over :math:`\vec{k}`

These steps are performed sequentially, such that Laplace inversion is first followed by Fourier inversion.
The Spakowitz lab has utilized other numerical approaches for these steps.
For example, numerical Laplace inversion can also be performed by numerical
complex integration of the Bromwich integral, but this approach is not integrated into the 'wlcgreen' module.

The orientation-independent Green's function :math:`G(\vec{K};p)`
in :eq:`gwlc` serves as a useful example for our numerical procedures.
Below, we discuss the two steps for the real-space inversion for this function.

Step 1: Laplace-space inversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Laplace-space inversion requires us to determine the poles :math:`\mathcal{E}_{\alpha})`
(all simple and generally complex).
Mathematically, this is written as

.. math::
    j (K; p = \mathcal{E}_{\alpha}) = \frac{1}{G(K; p = \mathcal{E}_{\alpha})} = 0

for :math:`\alpha = 0, 1, \ldots` gives the infinite set of poles within the infinite continued
fraction.  In practice, the real part of :math:`\mathcal{E}_{\alpha}` becomes increasingly negative
with increasing :math:`\alpha`, and a cutoff :math:`\alpha = \alpha_{\mathrm{max}}` is invoked for
any given :math:`N`. We note that :math:`\mathcal{E}_{\alpha} \rightarrow -\alpha (\alpha+D-2)`
(unperturbed eigenvalues of the hyperspherical harmonics) as
:math:`K \rightarrow 0` and becomes more negative with increasing :math:`K` (see Ref [Mehraeen2008]_).
We determine :math:`\alpha_{\mathrm{max}}` based on the :math:`K=0` magnitude of :math:`G` for a given :math:`N`.

The Laplace inversion is then determined by evaluating the residues in the complex plane.
The is written as

.. math::
    G (K; N) = \sum_{\alpha = 0}^{\alpha_{\mathrm{max}}}
    \mathop{\mathrm{Res}}_{p = \mathcal{E}_{\alpha}} \left[
    G(K; p) \exp \left( p N \right) \right]
    = \sum_{\alpha = 0}^{\alpha_{\mathrm{max}}} \lim_{p \rightarrow \mathcal{E}_{\alpha}}
    \left[
    ( p - \mathcal{E}_{\alpha}) G(K; p) \exp \left( p N \right) \right] =
    \sum_{\alpha = 0}^{\alpha_{\mathrm{max}}} \frac{1}{\partial_{p} j (K; p = \mathcal{E}_{\alpha})}
    \exp \left( \mathcal{E}_{\alpha} N \right)

The final form arises from the fact that all of the poles :math:`\mathcal{E}_{\alpha}` are simple.
From this development, the Laplace inversion requires evaluation of two quantities:
:math:`\mathcal{E}_{\alpha}` (for :math:`\alpha = 0, 1, \ldots, \alpha_{\mathrm{max}}`)
and :math:`\partial_{p} j (K; p = \mathcal{E}_{\alpha}) = \frac{\partial j (K; p = \mathcal{E}_{\alpha})}{\partial p}`.
We leverage several numerical and asymptotic methods to find the poles :math:`\mathcal{E}_{\alpha}`
(see [Mehraeen2008]_ for a detailed description of these methods).

This basic approach is extended to Green's function evaluation with fixed end orientation
(i.e. the orientation-dependent Green's function).
This necessitates determination of poles and residues that depends
on an additional index :math:`\mu`, which defines the :math:`z`-component of the
angular momentum (noting the Quantum Mechanical analogy of our problem).
Inclusion of :math:`\mu` requires definition of the cutoff :math:`\mu = \mu_{\mathrm{max}}`
that restricts the otherwise infinite summations within the inversion.

We define a data structure that facilitates evaluation of the Laplace inversion
but also enables the 'wlcgreen' module to dynamically grow the data structure,
so residues and poles do not require re-evaluation if already evaluated at a given
:math:`K`.
We define the 'green_laplace' library as follows

.. code:: python

    green_laplace = { k_val1:[mu_zero_only, lam_zero_only, poles, residues], \
                      k_val2:[mu_zero_only, lam_zero_only, poles, residues], ... }

For a given key in the library 'k_val', the elements are as follows

- 'mu_zero_only' is a boolean (True or False) that determines whether the poles and residues are evaluated for non-zero values of :math:`\mu`

- 'lam_zero_only' is a boolean (True or False) that determines whether the residues are evaluated for non-zero :math:`\lambda`

- 'poles' is an array of poles [size :math:`(\alpha_{\mathrm{max}} + 1) \times (\mu_{\mathrm{max}} + 1)`], which accounts for the value of 'mu_zero_only'

- 'residues' is an array of residues [size :math:`(\alpha_{\mathrm{max}} + 1) \times (\mu_{\mathrm{max}} + 1)`] (to be discussed further)

Step 2: Fourier-space inversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Functions contained with the 'wlcgreen' module
----------------------------------------------

.. automodule:: wlcstat.wlcgreen
    :members:
