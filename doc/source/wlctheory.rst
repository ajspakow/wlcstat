.. _wlctheory:

Wormlike Chain Theory
======================

The wormlike chain model describes the polymer chain as a continuous
thread through space that resists bending deformation~\cite{kratky-wlc}.
Subjecting the chain to thermal fluctuations results in a wiggling wormlike
polymer chain.

The polymer conformation is defined by a space curve :math:`\vec{r}(s)` where
:math:`s` runs from zero to the contour length :math:`L`.  The wormlike chain model has
a fixed length (like the freely jointed chain model), which is
enforced by setting

.. math::
   \left( \frac{\partial \vec{r}(s)}{\partial s} \right)^{\! \! 2} = 1

at all values of :math:`s` along the chain.
The chain opposes bending deformation, and the bending energy is 
quadratic in the curvature of the chain.  This gives the bending energy

.. math::
   E_{bend} = \frac{k_{B}T l_{p}}{2} \int_{0}^{L} ds \left(
   \frac{\partial^{2} \vec{r}(s)}{\partial s^{2}} \right)^{\! \! 2}

Wormlike Chain model Green function---Path integral formulation and the Diffusion equation
-----------------------------------------------------------------------------------------------

The governing differential equation for the Wormlike Chain model~\cite{kratky-wlc} is
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
orientation :math:`\vec{u}` at length :math:` L ` given that the chain begins at
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

