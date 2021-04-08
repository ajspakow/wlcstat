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

We develop a path integral formulation of the Active-Brownian particle that results in an expression for the
joint probability :math:`\mathcal{P}[x(t)|x_0;t;f_A^{(0)}]`.  This function governs the probability that if a particle begins
at :math:`x_{0}` experiencing an active force :math:`f_{A}^{(0)}` at time :math:`t = 0`, the particle will be
located at position :math:`x` at time :math:`t`.
Carrying out the integral over the active forces is performed by noting a Gaussian form of active forces

.. math::
    \mathcal{P}[f_A(t)]\propto \exp{\Big(-\frac{1}{2}\int_0^t dt_{1}\int_0^t dt_2
    f_A(t_1)\kappa_A^{-1}(|t_1-t_2|)f_A(t_2)\Big)}

This is used in the path integral formulation, and after functional integration over Brownian and Active forces, we
arrive at the expression

.. math::
    &  &
    \mathcal{P}[x(t)|x_0;t;f_A^{(0)}] =  \\
    &  &
    \int_{x_0}^{x(t)} \mathcal{D}[x(t)]\int \mathcal{D}[w(t)] \int d\eta \exp{\Big(-k_BT\xi\int_0^t [w(t^\prime)]^2dt^\prime+i\xi\int_0^t \dot{w}(t^\prime)x(t^\prime)dt^\prime\Big)}
    \\
    &  &
    \times\exp{\Big(-i\xi w(t)x(t)+i\xi w(0)x(0)+i\int_0^t w(t^\prime)f_V[x(t^\prime)]dt^\prime-i\eta f_A^{(0)}\Big)}
    \\
    &  &
    \times\exp{\Big(-\frac{1}{2}\int_0^t dt_1\int_0^t dt_2 w(t_1)\kappa_A(|t_1-t_2|)w(t_2)-\eta\int_0^t dt_1 \kappa_A(|t_1-t_0|)w(t_1)-\frac{1}{2}\eta^2\kappa_A(0)\Big)}

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
    C_{p}(\tau) =
    \langle
    \vec{X}_{p}(\tau) \cdot
    \vec{X}_{p'} (0)
    \rangle
    =
    \frac{3 k_{B}T}{k_{p}}
    \left\{
    \exp \left(-\frac{3 \pi^{2} p^{2} \tau}{N^2} \right)
    + \frac{\Gamma K_{A}^{2} }{K_A^2-(3 \pi^{2})^{2} p^4/N^{4}}
    \left[
    \exp \left(-\frac{3 \pi^{2} p^{2} \tau}{N^2} \right) -
    \frac{3 \pi^{2} p^{2}}{N^{2} K_{A}} \exp \left( - K_{A} \tau \right)
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