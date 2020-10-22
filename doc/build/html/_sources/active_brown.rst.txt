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
    & = & \kappa^{(B)} ( |t - t'|)
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
Our approach is amenable to address
such physics, and we develop this in
the appendix to this paper.

In our work, we assume the active force :math:`f_{A}` are also Gaussian-distributed
with an arbitrary temporal correlation, such that

.. math::
    \langle f_{A}(t) \rangle & = & 0, \\
    \langle f_{A}(t) f_{A}(t') \rangle
    & = & \kappa^{(A)}(|t-t'|),

where
:math:`\kappa^{(A)}(t)` represents
the temporal correlation between active
forces.

This theoretical model results in
a Smoluchowski equation (or
Schr\"{o}dinger equation) that is given by

.. math::
    \frac{\partial G(\vec{r}|\vec{r}_{0}; t)}
    {\partial t} & = &
    \frac{k_{B}T}{\xi} \nabla^{2}
    G(\vec{r}|\vec{r}_{0}; t)
    - \frac{1}{\xi} \vec{\nabla}
    \cdot
    \left[
    \vec{f}_{V} G(\vec{r}|\vec{r}_{0}; t)
    \right]
    \\
    &  & +
    \frac{1}{2 \xi^{2}}
    \int_{0}^{t} dt_{1} \int d \vec{r}_{1}
    K_{A}^{-1}(t-t_{1})
    \left[
    \vec{\nabla}_{r}
    G(\vec{r}|\vec{r}_{1}; t- t_{1})
    \right]
    \cdot
    \left[
    \vec{\nabla}_{r_{1}}
    G(\vec{r}_{1}|\vec{r}_{0}; t_{1})
    \right],
    :label: diff

which is composed of the thermal and
active contributions.

The active-Brownian Smoluchowski equation
(Eq. :eq:`diff`)
defines the conditions for the steady-state
probability distribution :math:`p_{ss}(\vec{r})`.
Upon Laplace transform of Eq. :eq:`diff`
and applying Final Value Theorem,
we arrive at the condition

.. math::
    0 =
    H_{0} p_{ss}(\vec{r}) +
    \frac{1}{2 \xi^{2}} \int d \vec{r}_{1}
    \left[
    \vec{\nabla}_{r}
    \tilde{G}_{A}(\vec{r}|\vec{r}_{1})
    \right]
    \cdot
    \left[
    \vec{\nabla}_{r_{1}} p_{ss} (\vec{r}_{1})
    \right],
    :label: pss

where :math:`H_{0} p_{ss}(\vec{r})` is the
standard Smoluchowski operator

.. math::
    H_{0} p_{ss}(\vec{r}) =
    \frac{k_{B}T}{\xi} \nabla^{2}
    p_{ss}(\vec{r})
    - \frac{1}{\xi} \vec{\nabla}
    \cdot
    \left[
    \vec{f}_{V} p_{ss}(\vec{r})
    \right]

and :math:`\tilde{G}_{A}` is defined as

.. math::
    \tilde{G}_{A}(\vec{r}|\vec{r}_{1}) =
    \int_{0}^{\infty} \! \! dt K_{A}^{-1}(t)
    G(\vec{r}|\vec{r}_{1};t).

The steady-state condition :math:`p_{ss}`
requires evaluation of Eq. :eq:`diff`,
which implies the non-equilibrium
nature of the active-Brownian particle.

Our current development is a generalized
representation of an active-Brownian particle.
However, the active-force correlation
function :math:`K_{A}^{-1}` defines a time-scale
of communication of active forces

.. math::
    t_{A} = \frac{\int_{0}^{\infty}
    d t K_{A}^{-1}(t) t }{\int_{0}^{\infty}
    d t K_{A}^{-1}(t)}

As this time scale approaches zero,
the active-force correlation can
be written as
:math:`K_{A}^{-1}(t) = 2 \xi k_{B} T_{A}
\delta (t)`,
where we define an active temperature
:math:`T_{A}`.
This limit implies the effect of active forces
can be interpreted as an effective
active-Brownian temperature
:math:`T_{AB} = T + T_{A}`.

However, the time scale of relaxation of the
particle in the absence of active forces
defines whether the
active forces perturb the behavior
beyond this effective temperature.
This is particularly true for a system
with a spectrum of relaxation time
scales, where some relaxation processes
are faster than the active-force correlation
time :math:`t_{A}`.
Such effects are prevalent for
the dynamics of a polymer chain
subjected to both active and Brownian
fluctuations.

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

The mapping from the harmonic particle to
the pth normal mode involves replacing
the following variables:
:math:`\xi \rightarrow N \xi`,
:math:`k \rightarrow k_{p} = \frac{3 \pi^{2} k_{B}T}{b^{2} N} p^{2}`,
:math:`f_{A}^{2} \rightarrow N f_{A}^{2}`.
We then define the Rouse time
:math:`t_{R} = N^{2} b^{2} \xi /(3 \pi^{2} k_{B}T)`
as the longest internal relaxation time of the
polymer, \emph{i.e.} :math:`t_{R} = N\xi/k_{1}`
for the :math:`p=1` normal mode.

We defined the normal-mode correlation
function for the pth mode to be

.. math::
    C_{p}(\tau) =
    \langle
    \vec{X}_{p}(\tau) \cdot
    \vec{X}_{p'} (0)
    \rangle
    =
    \frac{3 k_{B}T}{k_{p}} \left[
    1 + \frac{F_{A}^{2}}{2} \frac{K_{A}}{K_{A}+p^{2}}
    \right] \exp \left( - p^{2} \tau \right)
    \delta_{pp'}

where :math:`F_{A}^{2} =
f_{A}^{2}/(k_{A} \xi k_{B}T)`
and
:math:`K_{A} = t_{R} k_{A} =
N^{2} b^{2} \xi k_{A} /(3 \pi^{2} k_{B}T)`.
The dimensionless time :math:`\tau = t/t_{R}`
is scaled by the Rouse time of the
polymer.

We find the center-of-mass
mean-square displacement to be

.. math::
    \mathrm{MSD}_{\mathrm{com}} =
    \langle
    \left(
    \vec{X}_{0}(t) - \vec{X}_{0}(0)
    \right)^{2}
    \rangle
    =
    \frac{k_{B} T}{k_{1}}
    \left\{
    6 \tau + 3 F_{A}^{2}
    \left[
    \tau - \frac{1}{K_{A}} +
    \frac{1}{K_{A}}
    \exp \left(
    - K_{A} \tau
    \right)
    \right]
    \right\}.
    :label: msd-com

We note that Eq. :eq:`msd-com`
has the short-time asymptotic behavior
:math:`\mathrm{MSD}_{\mathrm{com}} \rightarrow
6 [k_{B}T/(\xi N)] t
` as :math:`t \rightarrow 0`,
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
    \mathrm{MSD} = \langle \left( \vec{r}(n, t) - \vec{r}(n, 0) \right)^{2} \rangle
    :label: msd-active

for the :math:`n`th segment of the chain with total length :math:`N`.
In this work, we focus on the midpoint
motion at :math:`n=N/2`, but this is
easily extended to other polymer positions.
We insert our normal-mode representation into Eq. :eq:`msd-active`,
resulting in the expression

.. math::
    \mathrm{MSD} & = &
    \sum_{p=0}^{\infty} \sum_{p'=0}^{\infty}
    \langle
    \left( \vec{X}_{p}(t) - \vec{X}_{p}(0) \right) \cdot
    \left( \vec{X}_{p'}(t) - \vec{X}_{p'}(0) \right) \rangle
    \phi_{p}(n=N/2) \phi_{p'}(n=N/2)
    \nonumber
    \\
    & = &
        \mathrm{MSD}_{\mathrm{com}} +
    \sum_{p = 1}^{\infty} 12 \frac{k_{B}T}{k_{2p}}
    \left[ 1 + \frac{F_{A}^{2}}{2} \frac{K_{A}}{K_{A}+(2p)^{2}} \right]
    \left[ 1 - \exp \! \left( - \frac{k_{2p}}{N \xi} t \right) \right].

In addition, we define the
the mean-square change in distance (MSCD) for a polymer chain.
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
    \mathrm{MSCD}
        =  \sum_{p = 0}^{\infty} 48 \frac{k_{B}T}{k_{2p+1}}
    \left[ 1 + \frac{F_{A}^{2}}{2} \frac{K_{A}}{K_{A}+(2p+1)^{2}} \right]
    \left[ 1 - \exp \! \left( - \frac{k_{2p+1}}{N \xi} t \right) \right]
    \sin^{2} \left[ \frac{\pi (2p+1) \Delta}{N} \right]


Functions contained with the 'active_brown' module
--------------------------------------------------

.. automodule:: wlcstat.active_brown
    :members:


Example usage of 'msd_active'
-----------------------------

We show the solution for
the MSD for chain length :math:`N=100` for
an active-Brownian polymer.
The top and middle plots show the behavior for :math:`K_{A} = 10` with varying :math:`F_{A}^{2}`,
and the bottom plot shows the behavior for :math:`F_{A}^{2} = 100` with varying :math:`K_{A}`.
The long-time asymptotic behavior for the active-Brownian polymer (dashed) is given by the
long-time behavior of the center-of-mass MSD, which is given by
:math:`\mathrm{MSD}_{\mathrm{com}} \rightarrow
6 [k_{B}T/(\xi N)] (1 + F_{A}^{2}/2) t`.


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.active_brown import *

    num_t = 100
    t_val_0 = 1e-6
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    Ndel=25
    ka = 10
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    fa2_vector = np.logspace(-1,3,5)
    for i in range(len(fa2_vector)):
        fa = np.sqrt(fa2_vector[i])
        msd = msd_active(t_val, 1, ka, fa, N, 1, 40000)
        plt.loglog(t_val, msd, label = '$F_{A}^{2}$ = ' + str(fa2_vector[i]))
        msd_inf = 6 / N * t_val * (1 + 0.5 * fa ** 2)
        plt.loglog(t_val, msd_inf,':')
    plt.legend()
    plt.xlabel(r'$time, t$')
    plt.ylabel(r'$MSD$')
    plt.tight_layout()
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.active_brown import *

    num_t = 100
    t_val_0 = 1e-6
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    Ndel=25
    ka = 10
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    fa2_vector = np.logspace(-1,3,5)
    for i in range(len(fa2_vector)):
        fa = np.sqrt(fa2_vector[i])
        msd = msd_active(t_val, 1, ka, fa, N, 1, 40000)
        plt.loglog(t_val, msd / t_val, label = '$F_{A}^{2}$ = ' + str(fa2_vector[i]))
        msd_inf = 6 / N * t_val * (1 + 0.5 * fa ** 2)
        plt.loglog(t_val, msd_inf / t_val,':')
    plt.legend()
    plt.xlabel(r'$time, t$')
    plt.ylabel(r'$MSD/t$')
    plt.tight_layout()
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.active_brown import *

    num_t = 100
    t_val_0 = 1e-6
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    fa = 10
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    ka_vector = np.logspace(-2,2,5)
    for i in range(len(ka_vector)):
        ka = ka_vector[i]
        msd = msd_active(t_val, 1, ka, fa, N, 1, 40000)
        plt.loglog(t_val, msd, label = '$K_{A}$ = ' + str(ka))
    msd_inf = 6 / N * t_val * (1 + 0.5 * fa ** 2)
    plt.loglog(t_val, msd_inf,':')
    plt.legend()
    plt.xlabel(r'$time, t$')
    plt.ylabel(r'$MSD$')
    plt.tight_layout()
    plt.show()



Example usage of 'mscd_active'
------------------------------

We show the solution for
the MSCD for chain length :math:`N=100` and :math:`\Delta = 25` for
an active-Brownian polymer.
The top plot shows the behavior for :math:`K_{A} = 10` with varying :math:`F_{A}^{2}`,
and the bottom plot shows the behavior for :math:`F_{A}^{2} = 100` with varying :math:`K_{A}`.
The long-time asymptotic behavior for the Brownian-only polymer is found by noting

.. math::
    \mathrm{MSCD} = \langle \left( \Delta \vec{R}(t) - \Delta \vec{R}(0) \right)^{2} \rangle =
    \langle \Delta \vec{R}(t)^{2} \rangle
    - 2 \langle \Delta \vec{R}(t) \cdot \Delta \vec{R}(0) \rangle
    + \langle \Delta \vec{R}(0) ^{2} \rangle \rightarrow
    2 \langle \Delta \vec{R}^{2} \rangle

which leads to the asymptotic solutions

.. math::
    \mathrm{MSCD} \rightarrow 4 \Delta


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.active_brown import *

    num_t = 100
    t_val_0 = 1e-4
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    Ndel=25
    ka = 10
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    fa2_vector = np.logspace(-1,3,5)
    for i in range(len(fa2_vector)):
        fa = np.sqrt(fa2_vector[i])
        mscd = mscd_active(t_val, 1, ka, fa, Ndel, N, 1, 20000)
        plt.loglog(t_val, mscd, label = '$F_{A}^{2}$ = ' + str(fa2_vector[i]))
    mscd_inf = 2 * 2 * Ndel
    plt.legend()
    plt.loglog(t_val, mscd_inf + 0*t_val,'--')
    plt.xlabel(r'$time, t$')
    plt.ylabel(r'$MSCD$')
    plt.tight_layout()
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.active_brown import *

    num_t = 100
    t_val_0 = 1e-4
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    Ndel=25
    fa = 10
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    ka_vector = np.logspace(-2,2,5)
    for i in range(len(ka_vector)):
        ka = ka_vector[i]
        mscd = mscd_active(t_val, 1, ka, fa, Ndel, N, 1, 20000)
        plt.loglog(t_val, mscd, label = '$K_{A}$ = ' + str(ka))
    mscd_inf = 2 * 2 * Ndel
    plt.legend()
    plt.loglog(t_val, mscd_inf + 0*t_val,'--')
    plt.xlabel(r'$time, t$')
    plt.ylabel(r'$MSCD$')
    plt.tight_layout()
    plt.show()