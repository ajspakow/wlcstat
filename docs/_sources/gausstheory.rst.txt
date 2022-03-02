.. _gausstheory:

Chain statistics as a path integral
===================================

Consider a single polymer chain, defined by a space-curve representation.
The probability distribution (or Green function) :math:`G(b|a)` of the chain extending from :math:`a` to :math:`b` is

.. math::
    G(b|a) = \sum_{all \; paths} \phi [x(s)]

where the state :math:`a` denotes location :math:`x_{a}` at
arclength :math:`s_{a}` and state :math:`b` denotes location :math:`x_{b}`
at arclength :math:`s_{b}`.
The summation over \emph{all paths} indicates a sum of all possible chain conformations (or paths) that
begin at state :math:`a` and end at state :math:`b`.
This concept of a path integral is extensively developed by Feynman and Hibbs~\cite{feynman}
to describe the statistical treatment of quantum mechanical particles.
For polymer statistics, the weighting of a path (or chain conformation) is governed by statistical
mechanics.
The statistical weight :math:`\phi[x(s)]` of a path :math:`x(s)` is given by the
integration of the Boltzmann weighting of the path.
Define the Action
:math:`S` for the particular path to be

.. math::
    S[x(s)] = \int_{s_{a}}^{s_{b}} \mathcal{H}(x, x', x'', ...;s) ds

Thus, the statistical weight of a path is given by

.. math::
    \phi[x(s)] = \frac{1}{A} \exp (-\beta S[x(s)]),

where :math:`A` is a normalization constant to be discussed later.

The Hamiltonian :math:`\mathcal{H}` includes differential
elastic deformation energy :math:`e_{elas}(x',x'',...;s)`
and a potential energy per unit length :math:`V(x,s)`.
The Green function :math:`G(b|a)` is written in %the form
path-integral form

.. math::
    G(b|a) = \int_{a}^{b} \exp (-\beta S[x(s)]) \mathcal{D} x(s)

Currently the derivation of the path model is in
one dimension; however, the ideas are easily extended to higher
dimensions.

.. figure:: figures/path-int.pdf
    :width: 600
    :align: center
    :alt: Schematic of path summation in the space-arclength plane

    Schematic of path summation in the space-arclength plane



Gaussian Chain model Green function---explicit path integration
---------------------------------------------------------------


Consider a polymer chain defined by the Gaussian Chain model
in 1-D with no external field (:math:`V=0`).
The action is given by

.. math::
    S[x(n)]=
    \frac{k_{B}T}{2 b^{2}}
    \int_{n_{a}}^{n_{b}} \! dn \, \left( \frac{\partial x(n)}{\partial n} \right)^{\! 2}

where :math:`n=s/b` is a dimensionless arclength parameter.
This gives a statistical weight

.. math::
    \phi[x(n)] = \frac{1}{A}
    \exp \! \left[
    -\frac{1}{2 b^{2}}
    \int_{n_{a}}^{n_{b}} \! dn \, \left( \frac{\partial x(n)}{\partial n} \right)^{\! 2}
    \right]


For this treatment, convert the path integral into a discrete-chain
representation.  This procedure provides a practical method to
perform the summation over all possible chain conformations.
Break the chain into :math:`M` segments of length :math:`\delta = (n_{b}-n_{a})/M`.
The position of the :math:`m`th point in the chain is :math:`x_{m}`,
and the endpoints are :math:`x_{0}=x_{a}` and :math:`x_{M}=x_{b}`.
The limit :math:`M \rightarrow \infty` recovers our continuous chain representation

.. figure:: figures/discrete-path.pdf
    :width: 600
    :align: center
    :alt: Discrete-chain representation converts the chain into :math:`M` segments (in this case :math:`M=6`) of length :math:`\delta=(n_{b}-n_{a})/M`.

    Discrete-chain representation converts the chain into :math:`M` segments (in this case :math:`M=6`) of length :math:`\delta=(n_{b}-n_{a})/M`.

The statistical weight is now approximated as

.. math::
    \phi[x(n)] \approx \frac{1}{A}
    \exp \! \left[
    -\frac{1}{2 b^{2} \delta}
    %\int_{n_{a}}^{n_{b}} \! dn \, \left( \frac{\partial x(n)}{\partial n} \right)^{\! 2}
    \sum_{m=1}^{M} \left( x_{m}-x_{m-1} \right)^{2}
    \right]

and the Green function is approximated as

.. math::
    G(b|a) = \frac{1}{A}
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{1}
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{2}
    \ldots
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{M-1}
    %\left( \prod_{m'=1}^{M-1} d x_{m'} \right)
    \exp \! \left[
    -\frac{1}{2 b^{2} \delta}
    %\int_{n_{a}}^{n_{b}} \! dn \, \left( \frac{\partial x(n)}{\partial n} \right)^{\! 2}
    \sum_{m=1}^{M}
    \left( x_{m}-x_{m-1} \right)^{2}
    \right].

The :math:`x_{1}` integral is given by

.. math::
    \normalsize
    &  &
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{1}
    \exp \!  \left\{
    -\frac{1}{2 b^{2} \delta}
    \left[
    \left( x_{2}-x_{1} \right)^{2}+
    \left( x_{1}-x_{0} \right)^{2}
    \right]
    \right\}=
    \nonumber \\
    \normalsize
    &  &
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{1}
    \exp \! \left\{
    -\frac{1}{2 b^{2} \delta}
    \left[
    2 \left( x_{1}-x_{1}^{\star} \right)^{2}+
    \frac{1}{2} \left( x_{2}-x_{0} \right)^{2}
    \right]
    \right\}= \nonumber \\
    &  &
    \sqrt{\pi b^{2} \delta} \exp \! \left[
    - \frac{1}{2 b^{2} \delta} \frac{1}{2}(x_{2}-x_{0})^{2}
    \right],

where :math:`x_{1}^{\star}=\frac{x_{0}+x_{2}}{2}`.

The subsequent :math:`x_{2}` integral is given by

.. math::
    \normalsize
    &  &
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{2}
    \exp \! \left\{
    -\frac{1}{2 b^{2} \delta}
    \left[
    \left( x_{3}-x_{2} \right)^{2}+
    \frac{1}{2} \left( x_{2}-x_{0} \right)^{2}
    \right]
    \right\}=
    \nonumber \\
    \normalsize
    &  &
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{2}
    \exp \! \left\{
    -\frac{1}{2 b^{2} \delta}
    \left[
    \frac{3}{2} \left( x_{2}-x_{2}^{\star} \right)^{2}+
    \frac{1}{3} \left( x_{3}-x_{0} \right)^{2}
    \right]
    \right\}= \nonumber \\
    &  &
    \sqrt{\frac{4\pi b^{2} \delta}{3}} \exp \! \left[
    - \frac{1}{2 b^{2} \delta} \frac{1}{3}(x_{3}-x_{0})^{2}
    \right]

where :math:`x_{2}^{\star}=\frac{x_{0}+2x_{3}}{3}`.
The :math:`x_{m}` integral after performing integrals :math:`x_{1}`, :math:`x_{2}`,\ldots, :math:`x_{m-1}`
is

.. math::
    &  &
    \int_{-\infty}^{\infty} \! \! \! \! \! d x_{m}
    \exp \! \left\{
    -\frac{1}{2 b^{2} \delta}
    \left[
    \left( x_{m+1}-x_{m} \right)^{2}+
    \frac{1}{m} \left( x_{m}-x_{0} \right)^{2}
    \right]
    \right\}=
    \nonumber \\
    &  &
    \sqrt{\frac{2\pi b^{2} \delta}{1+\frac{1}{m}}} \exp \! \left[
    - \frac{1}{2 b^{2} \delta} \frac{1}{m+1}(x_{m+1}-x_{0})^{2}
    \right]

This gives the Green function

.. math::
    G(b|a) & = & \frac{1}{A}
    \prod_{m=1}^{M-1} \sqrt{\frac{2 \pi b^{2} \delta}{1+\frac{1}{m}}}
    \exp \! \left[
    - \frac{1}{2 b^{2} \delta} \frac{1}{M} ( x_{M}-x_{0} )^{2}
    \right]
    \nonumber \\
    & = & \frac{1}{A}
    \prod_{m=1}^{M-1} \sqrt{\frac{2 \pi b^{2} \delta}{1+\frac{1}{m}}}
    \exp \! \left[
    - \frac{( x_{b}-x_{a} )^{2}}{2 (n_{b}-n_{a}) b^{2}}
    \right]
    \nonumber \\
    & = &
    \frac{1}{\sqrt{2 \pi (n_{b}-n_{a}) b^{2} }}
    \exp \! \left[
    - \frac{( x_{b}-x_{a} )^{2}}{2 (n_{b}-n_{a}) b^{2}}
    \right]

where :math:`A`
is set to ensure :math:`\int dx_{b} G(b|a)=1`

Gaussian Chain model Green function---Schrödinger equation
--------------------------------------------------------------

Generally, it is more convenient to develop the path integral into
a diffusion equation (or Schr\"{o}dinger equation) to solve for the chain statistics.
Consider the Gaussian Chain model in 1-D with an external potential,
giving an action

.. math::
    S[x(n)]=
    \int_{0}^{N} \! dn \,
    \left[
    \frac{k_{B}T}{2 b^{2}}
    \left( \frac{\partial x(n)}{\partial n} \right)^{\! 2}
    + V[x(n)]
    \right]

where we set :math:`n_{a}=0` and :math:`n_{b}=N` (without loss of generality),
and assume :math:`V=V(x)` only.
The chain obeys Markovian statistics; thus, the Green
function for a chain that starts at :math:`A` (\emph{i.e.} position :math:`x_{A}` and arclength :math:`n_{A}`)
to state :math:`C` (\emph{i.e.} :math:`x_{C}, n_{C}`)
through intermediate state :math:`B` (\emph{i.e.} :math:`x_{B}, n_{B}`) is

.. math::
    G(x_{C},n_{C} |x_{A},n_{A}) =
    \int_{-\infty}^{\infty} \! \! \! dx_{B}
    G(x_{C}, n_{C} | x_{B}, n_{B})
    G(x_{B},n_{B} |x_{A},n_{A}).


.. figure:: figures/markovchain.pdf
    :width: 600
    :align: center
    :alt: Schematic of a chain that obeys Markovian statistics

    Schematic of a chain that obeys Markovian statistics

Take the state :math:`C` to differ an infinitesimal length from state :math:`B` by setting
:math:`n_{C}=N+\delta` and :math:`n_{B}=N` (:math:`n_{A}=0`).  For simplicity, set
:math:`x_{A}=x_{0}`,
:math:`x_{B}=y`,
:math:`x_{C}=x`, giving

.. math::
    G(x,N+\delta |x_{0},0) =
    \int_{-\infty}^{\infty} \! \! \! dy
    G(x, N+\delta | y,N)
    G(y,N |x_{0},0).

In the limit :math:`\delta \rightarrow 0`, the Green function
:math:`G(x, N+\delta | y,N)` will be dominated by a single path whose action is

.. math::
    S \approx \frac{k_{B} T}{2 b^{2}} \delta
    \left( \frac{x-y}{\delta} \right)^{2} + \delta \, V \! \! \left( \frac{x+y}{2} \right).

Assume the contribution to the integral is mainly near :math:`y = x + \eta`
where :math:`\eta \ll 1`.  The Green function is thus given by

.. math::
    G(x, N + \delta|x_{0},0) = \int_{-\infty}^{\infty} \! \! d \eta
    \frac{1}{A}
    \exp \left( - \frac{1}{2 b^{2}} \frac{ \eta^{2} }{ \delta } \right)
    \exp [- \beta \delta V( x + \eta/2)] G(x + \eta, N | x_{0},0).

The integrand decays to zero for :math:`\eta > \sqrt{2 b^{2} \delta }`
thus the integrand must be expanded to linear order in :math:`\delta`
and quadratic order in :math:`\eta` in order to capture the lowest order
behavior of the Green function.
The expansion of the Green function
yields the solvable integral

.. math::
    &  &
    G(x, N|x_{0},0) + \delta \frac{\partial G(x,N|x_{0},0)}{\partial N}  =
    \int_{-\infty}^{\infty} \! \! d \eta
    \frac{1}{A}
    \exp \left( - \frac{1}{2 b^{2}} \frac{ \eta^{2} }{ \delta } \right)
    \nonumber \\
    &  &
    \hspace{0.4in}
    \times [1 - \beta \delta V(x)]
    \left( G (x,N|x_{0},0) + \eta \frac{\partial G(x,N|x_{0},0)}{\partial x}
    + \frac{\eta^{2}}{2} \frac{\partial^{2} G(x,N|x_{0},0}{\partial x^{2}} \right).

The solution for the normalization constant :math:`A` is determined by
solution of the leading term on the right hand side:

.. math::
    \frac{1}{A} \int_{-\infty}^{\infty} \exp
    \left(- \frac{1}{2 b^{2}} \frac{ \eta^{2} }{ \delta } \right) d \eta = 1

thus :math:`A = b \sqrt{ 2 \pi \delta}`.
The two additional
integrals are given by

.. math::
    \frac{1}{A} \int_{-\infty}^{\infty} \eta
    \exp \left(- \frac{1}{2 b^{2}} \frac{ \eta^{2} }{ \delta } \right) d \eta
    & = & 0 \\
    \frac{1}{A} \int_{-\infty}^{\infty} \eta^{2}
    \exp \left(- \frac{1}{2 b^{2}} \frac{ \eta^{2} }{ \delta } \right) d \eta
    & = & b^{2} \delta.


The Green function is given by

.. math::
    &  &
    G(x,N|x_{0},0) + \delta \frac{\partial G(x,N|x_{0},0)}{\partial N} = \nonumber \\
    &  &
    G(x,N|x_{0},0) -
    \delta \beta V G(x,N|x_{0},0) + \delta \frac{b^{2}}{2}
    \frac{ \partial^{2} G(x,N|x_{0},0)}{\partial x^{2}}.

Thus, the Green function for the 1-D Gaussian Chain obeys

.. math::
    \frac{\partial G(x,N|x_{0},0)}{\partial N} = -\beta V G(x,N|x_{0},0) + \frac{b^{2}}{2}
    \frac{ \partial^{2} G(x,N|x_{0},0)}{\partial x^{2}},

with initial condition :math:`G(x,0|x_{0},0)=\delta(x-x_{0})`.
For :math:`V(x)=0`, the solution for the Green function is

.. math::
    G(x,N|x_{0},0) =
    \frac{1}{\sqrt{2 \pi N b^{2} }}
    \exp \! \left[
    - \frac{( x-x_{0} )^{2}}{2 N b^{2}}
    \right]

which agrees with our previous derivation (Sec.~\ref{chap:gaussian-explicit}).

The governing differential equation for the Gaussian Chain in 3D is written as an
extension to our 1D equation as

.. math::
    \frac{\partial G(\vec{R},N|\vec{R}_{0},0)}{\partial N} =
    -\beta V G(\vec{R},N|\vec{R}_{0},0) + \frac{b^{2}}{6}
    \vec{\nabla}^{2} G(\vec{R},N|\vec{R}_{0},0)

with initial condition :math:`G(\vec{R},0|\vec{R}_{0},0)=\delta(\vec{R}-\vec{R}_{0})`.
The solution of this differential equation subject to a given
potential :math:`V` gives the equilibrium probability distribution
for the Gaussian Chain model.
For :math:`V(\vec{R})=0`, the solution for the Green function
for a Gaussian Chain in 3-D is

.. math::
    G(\vec{R},N|\vec{R}_{0},0) =
    \left( \frac{2\pi N b^{2}}{3} \right)^{\! \! \! -3/2}
    \exp \! \left[ - \frac{3 (\vec{R} - \vec{R}_{0} )^{2}}{2Nb^{2}} \right]

