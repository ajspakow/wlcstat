.. _active_brown:
.. .. automodule:: active_brown

Active-Brownian Dynamics
========================

Consider a particle undergoing
random motion over a potential
:math:`V`. For our development, we
focus on 1-dimensional transport,
but the final results are extended
to 3 dimensions.
The particle is subjected to
thermal (\emph{i.e.} Brownian) force
:math:`f_{B}` and active forces :math:`f_{A}` that
represent transient perturbations from the
surrounding enzyme activity.
The temporal evolution of the
position :math:`x(t)` is governed by a
modified Langevin equation

.. math::
    \xi \frac{d x}{d t}
    = f_{V} + f_{B} + f_{A},
    :label: langevin

where
:math:`f_{V} = - \frac{\partial V}{\partial x}`
is the potential force on the particle at
position :math:`x`.

The Brownian force :math:`f_{B}` is governed by the fluctuation dissipation theorem,
which states that :math:`f_{B}` is a
Gaussian-distributed random force
that satisfies

.. math::
    \langle f_{B}(t) \rangle & = & 0, \\
    \langle f_{B}(t) f_{B}(t') \rangle
    & = & \kappa_B ( |t - t'|)
    =
    2 k_{B}T \xi \delta (t - t').

This assumes that the environment is
purely Newtonian, leading to the
instantaneous decorrelation of Brownian
forces (i.e. Gaussian white noise).
Diffusive transport in a viscoelastic
fluid leads to temporal memory in the
Brownian force, reflecting temporal
correlation in the frictional stress.

In our work, we assume the active force :math:`f_{A}` are also Gaussian-distributed
with an arbitrary temporal correlation, such that

.. math::
    \langle f_{A}(t) \rangle & = & 0, \\
    \langle f_{A}(t) f_{A}(t') \rangle
    & = & \kappa_{A}(|t-t'|),

where
:math:`\kappa_{A}(t)` represents
the temporal correlation between active
forces.

We develop a path integral formulation of the Active-Brownian particle that results in an expression for the
joint probability :math:`\mathcal{P}[x(t)|x_0;t;f_A^{(0)}]`.  This function governs the probability that if a particle begins
at :math:`x_{0}` experiencing an active force :math:`f_{A}^{(0)}` at time :math:`t = 0`, the particle will be
located at position :math:`x` at time :math:`t`.
Carrying out the integral over the active forces is performed by noting a Gaussian form of active forces

.. math::
    \mathcal{P}[f_A(t)]\propto \exp \left\{
    -\frac{1}{2}\int_0^t dt_{1}\int_0^t dt_2
    f_A(t_1)\kappa_A^{-1}(|t_1-t_2|)f_A(t_2)
    \right\}

This is used in the path integral formulation, and after functional integration over Brownian and Active forces, we
arrive at the expression

.. math::
    &  &
    \mathcal{P}[x|x_0;t,f_A^{(0)}, t_{0}]
    \nonumber
    \\
    &  &
    =\int_{x_0}^{x} \mathcal{D}[x(t)]\int \mathcal{D}[w(t)] \int d\eta \exp
    \left\{
    -k_{B}T\xi\int_0^t [w(t_{1})]^2dt_{1}+i\xi\int_0^t \dot{w}(t_{1})x(t_{1})dt_{1}
    \right.
    \nonumber
    \\
    &  &
    \left.
    -i\xi w(t)x(t)+i\xi w(0)x(0)+i\int_0^t w(t_{1})f_V[x(t_{1})]dt_{1}
    -\frac{1}{2}
    \int_0^t
    \! \! dt_1 \int_{0}^{t}
    \! \!
    dt_2 w(t_1)\kappa_A(|t_1-t_2|)
    w(t_2)
    \right.
    \nonumber
    \\
    &  &
    \left.
    -i\eta f_A^{(0)}
    -\eta\int_0^t dt_1 \kappa_A(|t_1-t_0|)w(t_1)-\frac{1}{2}\eta^2\kappa_A(0)
    \right\}

The theoretical development so far is general to any form of the spatially varying potential
:math:`V(x)` and does not assume any special form for the conservative forces between particles.

Dynamic behavior of an Active-Brownian polymer
----------------------------------------------

We apply our results of an
active-Brownian particle in a harmonic
potential to the problem of an
active-Brownian polymer chain.
Since the normal-mode expansion of the
Rouse polymer adopts the form of
a particle in a harmonic potential,
we can frame our results in the previous
section to find the behavior of the
pth normal mode.

We defined the normal-mode correlation
function for the pth mode to be

.. math::
    C_{p}(\tau) & = &
    \langle
    \vec{X}_{p}(\tau) \cdot
    \vec{X}_{p'} (0)
    \rangle
    \\
    & = &
    \frac{3 k_{B}T}{k_{p}}
    \left\{
    \exp \left(-\frac{ p^{2} \tau}{N^2} \right)
    + \frac{\Gamma}{1- p^4/(K_{A}^{2} N^{4})}
    \left[
    \exp \left(-\frac{p^{2} \tau}{N^2} \right) -
    \frac{ p^{2}}{K_{A} N^{2}} \exp \left( - K_{A} \tau \right)
    \right]
    \right\} \delta_{pp'}

The dimensionless time :math:`\tau = t/t_{b}`
is scaled by the diffusive time :math:`t_{b} = b^{2} \xi/(k_{B}T)`,
and we define the dimensionless Active force :math:`F_{A}^{2} = \Gamma K_{A} = f_{A}^{2} b^{2}/(k_{B}T)^{2}`
and dimensionless active rate constant
:math:`K_{A} = t_{b} k_{A}`.

We find the center-of-mass
mean-square displacement to be

.. math::
    \mathrm{MSD}_{\mathrm{com}} =
    \langle
    \left(
    \vec{X}_{0}(t) - \vec{X}_{0}(0)
    \right)^{2} =\frac{6Nk_BT}{k_b}\Big[
    \Big(\frac{1+\Gamma}{N^2}\Big)\tau+\frac{\Gamma}{N^2K_A}\Big(e^{-K_A \tau}-1\Big)
    \Big]
    :label: msd-com

We note that Eq. :eq:`msd-com`
has the short-time asymptotic behavior
:math:`\mathrm{MSD}_{\mathrm{com}} \rightarrow
6 [k_{B}T/(\xi N)] t` as :math:`t \rightarrow 0`,
which coincides with the mean-square
displacement for a Brownian-only polymer.
The long-time asymptotic behavior
:math:`\mathrm{MSD}_{\mathrm{com}} \rightarrow
6 [k_{B}T/(\xi N)] (1 + F_{A}^{2}/2) t
` as :math:`t \rightarrow \infty`
reveals the long-time effective
temperature due to Brownian and active
fluctuations, given by
:math:`T_{AB} = T (1 + F_{A}^{2}/2)`.

The mean-square displacement of a segment of the polymer chain is define as

.. math::
    \text{MSD}(\tau) = \langle ( \vec{r}(n,t) - \vec{r}(n,0))^{2} \rangle
    :label: msd-active

for the chain position n. In this work, we focus on the midpoint motion at :math:`n=N/2`, but this is easily extended
to other polymer positions. We insert our normal-mode representation into Eq. :eq:`msd-active`,
resulting in the expression

.. math::
    \text{MSD}(\tau)
     = \text{MSD}_{\text{com}}(\tau) + 4 \sum_{p=1}^{\infty} \triangle C_{2p}(\tau)


where :math:`\triangle C_{p}(\tau) = C_{p}(0) - C_{p}(\tau)`.
In addition, we define the mean-square change in distance (MSCD) for a polymer chain.
This quantity is defined as

.. math::
    \mathrm{MSCD} = \langle \left( \Delta
    \vec{R}(t) - \Delta \vec{R}(0) \right)^{2}
    \rangle

where
:math:`\Delta \vec{R}(t) = \vec{r}(N/2 + \Delta, t) - \vec{r}(N/2 - \Delta,t)`.
Thus, MSCD is a measure of the change
in separation distance between two
points on the polymer that are :math:`2 \Delta`
segments from each other.
This results in the expression

.. math::
    \text{MSCD}(\tau) =
    16 \sum_{p=0}^{\infty}
    \triangle C_{2p+1}(\tau)
    \sin^{2} \!
    \left[
    \frac{\pi(2p+1)\Delta}{2N}
    \right],

where :math:`\triangle C_{p}(\tau) = C_{p}(0) - C_{p}(\tau)`


Functions contained with the 'active_brown' module
--------------------------------------------------

.. automodule:: wlcstat.active_brown
    :members:
