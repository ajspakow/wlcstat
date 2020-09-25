.. _wlctheory:

Wormlike Chain Theory
=====================

The wormlike chain model describes the polymer chain as a continuous
thread through space that resists bending deformation [Kratky1949]_.
Subjecting the chain to thermal fluctuations results in a wiggling wormlike
polymer chain.

The polymer conformation is defined by a space curve :math:`\vec{r}(s)` where
:math:`s` runs from zero to the contour length :math:`L`.  The wormlike chain model has
a fixed length (like the freely jointed chain model), which is
enforced by setting

.. math::
    \left( \frac{\partial \vec{r}(s)}{\partial s} \right)^{\! \! 2} =
    \left( \vec{u}(s) \right)^{2}
    = 1

at all values of :math:`s` along the chain, where we define the tangent vector
:math:`\vec{u}(s) = \frac{\partial \vec{r}(s)}{\partial s}`.
The chain opposes bending deformation, and the bending energy is 
quadratic in the curvature of the chain.  This gives the bending energy

.. math::
    E_{bend} = \frac{k_{B}T l_{p}}{2} \int_{0}^{L} ds \left(
    \frac{\partial^{2} \vec{r}(s)}{\partial s^{2}} \right)^{\! \! 2} =
    \frac{k_{B}T l_{p}}{2} \int_{0}^{L} ds \left(
    \frac{\partial \vec{u}(s)}{\partial s} \right)^{\! 2}

.. figure:: figures/wlc-model.pdf
    :width: 600
    :align: center
    :alt: Schematic representation of the wormlike chain model and definition of geometric parameters.

    Schematic representation of the wormlike chain model and definition of geometric parameters.

Wormlike Chain model Green function---Path integral formulation and the Diffusion equation
------------------------------------------------------------------------------------------

The governing differential equation for the Wormlike Chain model [Kratky1949]_ is
given through the application of the path integral in the 
necessary dimensions.
The differential Hamiltonian :math:`\mathcal{H}` of the Wormlike
chain model is given by

.. math::
    \beta \mathcal{H} = \frac{l_{p}}{2} \left( \frac{\partial \vec{u}}{\partial s} \right)^{2} +\beta V(\vec{r}, \vec{u}),

where the tangent :math:`\vec{u} = \frac{\partial \vec{r}}{\partial s}`
(with constraint :math:`\left| \frac{\partial \vec{r}}{\partial s} \right|=1` for all :math:`s`).
Define the Green function :math:`G(\vec{r}, \vec{u}|\vec{r}_{0},\vec{u}_{0}; L)` to be the 
probability of finding the chain at location :math:`\vec{r}` with
orientation :math:`\vec{u}` at length :math:`L` given that the chain begins at
location :math:`\vec{r}_{0}` with orientation :math:`\vec{u}_{0}`.
The Markov property of the chain is written as

.. math::
    G (\vec{r}_{C}, \vec{u}_{C}, s_{C}|\vec{r}_{A},\vec{u}_{A}, s_{A}) = \int d \vec{r}_{B} d \vec{u}_{B} G (\vec{r}_{C}, \vec{u}_{C}, s_{C}|\vec{r}_{B},\vec{u}_{B}, s_{B}) \times G (\vec{r}_{B}, \vec{u}_{B}, s_{B}|\vec{r}_{A},\vec{u}_{A}, s_{A}).

Set :math:`s_{A}=0`, :math:`s_{B}=L`, and :math:`s_{C}=L+\delta` with 
:math:`\vec{r}_{A}=\vec{r}_{0}`,
:math:`\vec{r}_{B}=\vec{r}_{1}`, and
:math:`\vec{r}_{C}=\vec{r}` to find the governing diffusion equation
for :math:`G(\vec{r},\vec{u}|\vec{r}_{0},\vec{u}_{0};L)`.


In the limit :math:`\delta \rightarrow 0`, the Green function
:math:`G(\vec{r},\vec{u}, N+\delta | \vec{r},\vec{u},N)` will be dominated by a single path whose action is

.. math::
    S \approx \frac{l_{p}}{2} \delta \left( \frac{ \vec{u} - \vec{u}_{1}}{\delta} \right)^{2} + \delta \beta V \!\! \left( \frac{\vec{r}+\vec{r}_{1}}{2} \right).

The Green function for this path is thus

.. math::
    G(\vec{r}, \vec{u}, L + \delta| \vec{r}_{1}, \vec{u}_{1}, L) =
    \frac{\delta (\vec{r} - \vec{r}_{1} - \vec{u} \delta )}{A} \exp \!
    \left( \frac{l_{p}}{\delta} \vec{u} \cdot \vec{u}_{1} \right) \exp (-\delta \beta V)

where the normalization constant :math:`A` is to be 
determined in the limit of :math:`\delta` approaching zero, and
the delta function is added to the Green function to conserve the path length.
The Green function is now given by

.. math::
    G(\vec{r}, \vec{u}, L + \delta| \vec{r}_{0}, \vec{u}_{0}, 0)   =
    \int \frac{\delta (\vec{r} - \vec{r}_{1} - \vec{u} \delta )}{A}
    \exp \left( \frac{l_{p}}{\delta} \vec{u} \cdot \vec{u}_{1}
    \right)\exp (-\delta \beta V)  \times G(\vec{r}_{1}, \vec{u}_{1},L
    | \vec{r}_{0}, \vec{u}_{0},0) d \vec{r}_{1} d \vec{u}_{1}.

Assume only small deviations from :math:`\vec{r}` and :math:`\vec{u}`
will contribute to the integral.  Expand the Green function 
:math:`G(\vec{r}_{1}, \vec{u}_{1}, L| \vec{r}_{0},\vec{u}_{0},0)` about :math:`\vec{r}`
and :math:`\vec{u}` to quadratic order in the deviation to find the
leading order change in the propagated Green function.

The resulting expansion is given by

.. math::
    G(\vec{r}_{1}, \vec{u}_{1}, L| \vec{r}_{0},\vec{u}_{0},0)
    & = &
    G(\vec{r}, \vec{u}, L| \vec{r}_{0},\vec{u}_{0},0)-
    \frac{\partial G}{\partial \vec{r}} \cdot \Delta \vec{r}
    - \frac{\partial G}{\partial \vec{u}} \cdot \Delta \vec{u} \nonumber \\
    &   &
    \hspace{-1in}
    + \frac{1}{2} \frac{\partial^{2} G}{\partial \vec{r}^{2}}
    : \Delta \vec{r} \Delta \vec{r}
    + \frac{1}{2} \frac{\partial^{2} G}{\partial \vec{u}^{2}}
    : \Delta \vec{u} \Delta \vec{u}
    + \frac{\partial^{2} G}{\partial \vec{r} \partial \vec{u}}
    : \Delta \vec{r} \Delta \vec{u}

where :math:`\Delta \vec{r} = \vec{r} - \vec{r}_{1}` and
:math:`\Delta \vec{u} = \vec{u} - \vec{u}_{1}`.
The double dot operation between two tensors acts as :math:`\textbf{A}:\textbf{B} = A_{ij} B_{ji}`.
The normalization constant is found by requiring the leading term
to approach unity as the time change :math:`\delta` approaches zero.
The following integral determines :math:`A`:

.. math::
    \frac{1}{A} \int \delta(\vec{r} - \vec{r}_{1} - \vec{u} \delta )
    \exp \left( \frac{l_{p}}{\delta}
    \vec{u} \cdot \vec{u}_{1} \right)
    d \vec{u}_{1} d \vec{r}_{1} & = & \nonumber \\
    \hspace{-0.5in}
    \frac{1}{A} \int \delta(\vec{r} - \vec{r}_{1} - \vec{u} \delta )
    \left[ \int_{0}^{2 \pi} \int_{0}^{\pi}
    \exp \left( \frac{l_{p}}{\delta} \cos \theta \right) \sin \theta
    d \theta d \phi \right] d \vec{r}_{1} & = & \nonumber \\
    \frac{1}{A} \frac{4 \pi \delta}{l_{p}} \sinh
    \left( \frac{l_{p}}{\delta} \right) & = & 1

thus :math:`A = \frac{4 \pi \delta}{l_{p}} \sinh \left( \frac{l_{p}}{\delta} \right)`.
The leading order
response of the Green function is found by evaluating the averaged
quantities within the integrand.  Define the average of a given
quantity to be

.. math::
    \left< ... \right> = \int d \vec{u}_{1} \frac{1}{A}
    \exp (h \vec{u} \cdot \vec{u}_{1}) (...),

where :math:`h = \frac{l_{p}}{\delta}`.

We must find the
following quantities: :math:`\left< \Delta \vec{r} \right>`,
:math:`\left< \Delta \vec{u} \right>`, and
:math:`\left< \Delta \vec{u} \Delta \vec{u} \right>`.  The other averages
( :math:`\left< \Delta \vec{r} \Delta \vec{r} \right>` and
:math:`\left< \Delta \vec{r} \Delta \vec{u} \right>` ) yield
terms of quadratic order in :math:`\delta`, thus they do not contribute
to the lowest order response.
For convenience, define the vector
quantity :math:`\vec{h} = h \vec{u}` for use in this derivation.
The average quantity
:math:`\langle \Delta \vec{u} \rangle` is given by

.. math::
    \langle \Delta \vec{u} \rangle & = & \vec{u} - \langle \vec{u}_{1} \rangle \nonumber \\
    & = & \vec{u} - \frac{1}{A} \int d \vec{u}_{1} \exp
    ( \vec{h} \cdot \vec{u}_{1} ) \vec{u}_{1} \nonumber \\
    & = & \vec{u} - \frac{1}{A} %\vec{\nabla}_{h}
    \frac{\partial A}{\partial \vec{h}} \nonumber \\
    & = & \vec{u} - \frac{4 \pi \vec{u}}{A}
    \left( \frac{ \cosh (h) }{h} - \frac{\sinh (h)}{h^{2}} \right) \nonumber \\
    & = & \vec{u} (1 - \coth (h) + h^{-1} ).

Thus, in the limit of :math:`h \rightarrow \infty`, the quantity
:math:`\langle \Delta \vec{u} \rangle = \frac{\vec{u} \delta}{l_{p}}`.
The average quantity :math:`\langle \Delta \vec{u} \Delta \vec{u} \rangle`
is given by

.. math::
    \langle \Delta \vec{u} \Delta \vec{u} \rangle & = & \vec{u} \vec{u}
    - \vec{u} \langle \vec{u}_{1} \rangle - \langle \vec{u}_{1} \rangle \vec{u}
    + \langle \vec{u}_{1} \vec{u}_{1} \rangle \nonumber \\
    & = & \frac{2 \delta}{l_{p}} \vec{u} \vec{u}
    - \vec{u} \vec{u} + \langle \vec{u}_{1} \vec{u}_{1} \rangle \nonumber \\
    & = & \frac{2 \delta}{l_{p}} \vec{u} \vec{u}
    - \vec{u} \vec{u} + \frac{1}{A}
    \frac{\partial^{2} A}{\partial \vec{h} \partial \vec{h}}.

The following quantity is found to complete the derivation:

.. math::
    \frac{\partial^{2} A}{\partial \vec{h} \partial \vec{h}}
    = 4 \pi \mathbf{I}
    \left( \frac{\cosh (h)}{h^{2}} - \frac{\sinh (h)}{h^{3}} \right) +
    4 \pi \vec{h} \vec{h}
    \left( \frac{ \sinh (h)}{h^{3}} - \frac{3 \cosh (h)}{h^{4}} + \frac{3 \sinh (h)}{h^{5}} \right).

This is evaluated and the limit as :math:`h` approaches infinity yields
the solution of :math:`\langle \Delta \vec{u} \Delta \vec{u} \rangle
= \frac{\delta}{l_{p}} ( \mathbf{I} - \vec{u} \vec{u} )`.
Altogether the Green function obeys

.. math::
    G + \frac{\partial G}{\partial s} \delta = G - \beta V G \delta
    - \frac{\partial G}{\partial \vec{r}} \cdot \vec{u} \delta
    - \frac{\partial G}{\partial \vec{u}} \cdot  \frac{\vec{u}}{l_{p}} \delta
    + \frac{1}{2} \vec{\nabla}_{u} \vec{\nabla}_{u} G :
    \frac{\delta}{l_{p}} ( \mathbf{I} - \vec{u} \vec{u} )

This differential equation is further simplified by noting that
:math:`\vec{u} \cdot \vec{\nabla}_{u} = 0` and
:math:`\vec{u} \vec{u} : \vec{\nabla}_{u} \vec{\nabla}_{u} = 0`.
The final differential equation for the Green function is given by

.. math::
    \frac{\partial G(\vec{r}, \vec{u}, L| \vec{r}_{0}, \vec{u}_{0}, 0)}{\partial L} & = &
    - \beta V G(\vec{r}, \vec{u}, L| \vec{r}_{0}, \vec{u}_{0}, 0) \nonumber \\
    &  &
    - \vec{u} \cdot \vec{\nabla}_{r} G(\vec{r}, \vec{u}, L| \vec{r}_{0}, \vec{u}_{0}, 0) + \frac{1}{2 l_{p}}
    \vec{\nabla}_{u}^{2} G(\vec{r}, \vec{u}, L | \vec{r}_{0}, \vec{u}_{0}, 0),

with the initial condition :math:`G(\vec{r}, \vec{u}, 0 | \vec{r}_{0}, \vec{u}_{0}, 0)=
\delta(\vec{r}-\vec{r}_{0})\delta(\vec{u}-\vec{u}_{0})`.

Orientation statistics
----------------------

Consider the case :math:`V=0`.
Define the Fourier transform (and inverse transform) of a function :math:`f(\vec{r})` from :math:`\vec{r}` to :math:`\vec{k}` to be

.. math::
    \tilde{f}(\vec{k}) =
    \int d \vec{r} f(\vec{r}) \exp \left( i \vec{k} \cdot \vec{r} \right) %= \tilde{f}(\vec{k})
    \hspace{0.1in}
    \mathrm{and}
    \hspace{0.1in}
    f(\vec{r}) =
    \frac{1}{(2\pi)^{3}} \int d \vec{k} \tilde{f} ( \vec{k})\exp \left( -i \vec{k} \cdot \vec{r} \right)

The governing differential equation for :math:`\tilde{G}(\vec{k},\vec{u}|\vec{u}_{0};L)` for :math:`V=0` is

.. math::
    \frac{\partial \tilde{G}(\vec{k},\vec{u}|\vec{u}_{0};L)}{\partial L} =
    \left(i\vec{k} \cdot \vec{u} +
    \frac{1}{2 l_{p}} \vec{\nabla}_{u}^{2}
    \right)
    \tilde{G}(\vec{k},\vec{u}|\vec{u}_{0};L)

with initial condition
:math:`\tilde{G}(\vec{k},\vec{u}|\vec{u}_{0};0)=\exp (i\vec{k} \cdot \vec{r}_{0}) \delta(\vec{u}-\vec{u}_{0})`

The integral of the Green function over the position :math:`\vec{r}` is equivalent to

.. math::
    G(\vec{u}|\vec{u}_{0};L) = \int d \vec{r} G(\vec{r}, \vec{u}, L | \vec{r}_{0}, \vec{u}_{0}, 0)
    = \tilde{G}(\vec{k}=\vec{0},\vec{u}|\vec{u}_{0};L),

which gives the orientation-only chain statistics, i.e. the probability that a chain begin with
orientation :math:`\vec{u}_{0}` and ends with orientation :math:`\vec{u}` regardless of the end positions.
The orientation Green function :math:`G(\vec{u}|\vec{u}_{0};L)` satisfies

.. math::
    \frac{\partial G(\vec{u}|\vec{u}_{0};L)}{\partial L} =
    \frac{1}{2 l_{p}} \vec{\nabla}_{u}^{2}
    G(\vec{u}|\vec{u}_{0};L)
    :label: gudiffeqn
.. \label{eq:gudiffeqn}

with initial condition :math:`G(\vec{u}|\vec{u}_{0};L) = \delta(\vec{u}-\vec{u}_{0})`.

Here, we note that the operator :math:`\vec{\nabla}_{u}^{2}` has eigenfunctions :math:`Y_{l}^{m}` that satisfy

.. math::
    \vec{\nabla}_{u}^{2} Y_{l}^{m} = - l(l+1) Y_{l}^{m}.

The eigenfunctions :math:`Y_{l}^{m}` are the spherical
harmonics [Arfken1999]_ that form a complete basis set for the 3-dimensional angular Laplacian :math:`\vec{\nabla}_{u}^{2}`.
This basis set is extended to arbitrary dimensions as the hyperspherical harmonics,
and we will make use of this extension in coming chapters.
The range of the indices :math:`l` and :math:`m` are :math:`l \in [0, \infty]` and :math:`m \in [-l, l]`.
The spherical harmonics satisfy

.. math::
    \int d \vec{u} Y_{l}^{m} (\vec{u})
    Y_{l'}^{m'*}(\vec{u}) = \delta_{ll'} \delta_{mm'}

The solution for :math:`G(\vec{u}|\vec{u}_{0};L)` is constructed as an expansion in the spherical harmonics, such that

.. math::
    G(\vec{u}|\vec{u}_{0};L) = \sum_{l=0}^{\infty} \sum_{m=-l}^{l}
    Y_{l}^{m} (\vec{u})
    Y_{l}^{m*}(\vec{u}_{0}) \exp [-l(l+1)N]
    :label: guwlc
.. \label{eq:guwlc}

where :math:`N=L/(2l_{p})`.
This form satisfies the governing Schrödinger equation (Eq. :eq:`gudiffeqn`), and the initial condition
is captured by noting that

.. math::
    \sum_{l=0}^{\infty} \sum_{m=-l}^{l}
    Y_{l}^{m} (\vec{u})
    Y_{l}^{m*}(\vec{u}_{0}) = \delta \left( \vec{u} - \vec{u}_{0} \right).

This development enables the evaluations of average quantities (as discussed in the Average Quantities section).
