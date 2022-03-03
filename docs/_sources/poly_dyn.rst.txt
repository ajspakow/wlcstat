.. _poly_dyn:
.. .. automodule:: poly_dyn

Polymer Dynamics
================






Dynamics of a linear and ring flexible Gaussian chains (Rouse model)
--------------------------------------------------------------------

We consider a flexible polymer chain
that is subject
Brownian forces [Doi1988]_.
Our goals is to determine specific dynamic properties of a linear and ring
Gaussian chain. For example, we consider the dynamic motion of two
chain segments relative to each other.  We define this relative motion
as the mean-square change in displacment (MSCD).
In this development, we present a detailed derivation for the linear chain,
and we provide the results for the ring polymer.

.. figure:: figures/rouse-mscd.pdf
    :width: 600
    :align: center
    :alt: Schematic representation of the linear and ring Rouse model with definition of segments for MSCD calculation.

    Schematic representation of the linear and ring Rouse model with definition of segments for MSCD calculation.
    The total chain length is :math:`N = 2 N_{s}`, and we consider two points on the chain located at
    :math:`n_{1} = N_{s} - \Delta = N/2 - \Delta` and :math:`n_{2} = N_{s} + \Delta = N/2 + \Delta`.


We start by defining a discrete
polymer chain with bead positions
:math:`\vec{R}_{m}`, where :math:`m` runs from 0
to :math:`n_{\mathrm{b}}`.
Each bead is connected to their
neighboring beads by Hookean springs,
resulting in a potential force
on the :math:`m`-th bead
:math:`\vec{F}_{V_{m}}` that is given by


.. math::
    & & \vec{F}_{V_{m}}=
    \left\{ \begin{array}{c}
		\frac{3 k_{B}T}{g b^{2}} \left( \vec{R}_{m+1} - 2 \vec{R}_{m} + \vec{R}_{m-1} \right),
        \hspace{0.3in}
		m = 1,\ldots, n_{\mathrm{b}}-1 \\
		\frac{3 k_{B}T}{g b^{2}} \left( \vec{R}_{1} - \vec{R}_{0} \right),
        \hspace{2in}
		m = 0 \\
		-\frac{3 k_{B}T}{g b^{2}}   \left( \vec{R}_{n_{\mathrm{b}}} - \vec{R}_{n_{\mathrm{b}}-1} \right),
        \hspace{1.65in}
		m = n_{\mathrm{b}},
		\end{array}
	\right.

.. \label{eq:restoringforce}

where :math:`b` is the Kuhn statistical segment length of the polymer [Doi1988]_,
and :math:`g` is the number of Kuhn lengths
per bead.  We set :math:`g = N/n_{\mathrm{b}}`, giving
a total chain length of :math:`N` Kuhn lengths.

The Langevin equation of motion for the :math:`m`-th bead is given by

.. math::
    g \xi \frac{d\vec{R}_{m}(t)}{dt}=\vec{F}_{V_{m}} + \vec{F}_{B_{m}}

where the drag coefficient :math:`\xi` is defined as the viscous drag per Kuhn length.
We now take the limit of
:math:`n_{\mathrm{b}} \rightarrow \infty`
for fixed chain length :math:`N`, resulting in
a continuous-chain representation
of the chain :math:`\vec{r}(n,t)` where :math:`n`
is a path-length variable that runs
from :math:`0` to :math:`N`.
The Langevin equation for the chain is
given by

.. math::
    \xi \frac{\partial \vec{r}(n, t)}{\partial t}=
    \frac{3 k_{B}T}{b^{2}}
    \frac{\partial^{2} \vec{r}(n, t)}{\partial n^{2}}
    + \vec{f}_{B}(n,t),
    :label: eom-cont

which is subject to the end conditions

.. math::
    \frac{\partial \vec{r}(n=0,t)}{\partial n} =
    \frac{\partial \vec{r}(n=N,t)}{\partial n} = 0.

The Brownian forces satisfy the fluctuation dissipation theorem

.. math::
    \langle \vec{f}_{B}(n,t)\vec{f}_{B}(n',t')\rangle =
    2k_{B}T \xi \delta (n - n') \delta(t-t') \textbf{I}.


It is convenient to define a set of normal coordinates that
effectively decouple the interactions implicit within the equation of motion (Eq. :eq:`eom-cont`).
We define the normal modes

.. math::
    & & \phi_{p}(n)=
    \left\{ \begin{array}{c}
		\sqrt{2} \cos \left( \frac{\pi p n}{N} \right), \ \ p = 1,2,\ldots \\
		1, \ \ \ \ \ \ \ \ \ \  p = 0.
		\end{array}
	\right.

These modes represent a complete basis set that satisfy the boundary conditions for :math:`\vec{r}(n,t)`.
Orthogonality is demonstrated by noting

.. math::
    \int_{0}^{N} \! \! dn \phi_{p}(n) \phi_{p'}(n) = N \delta_{p,p'}.

The amplitude of the :math:`p`-th mode :math:`\vec{X}_{p}(t)` is given by

.. math::
    \vec{X}_{p}(t) = \frac{1}{N} \int_{0}^{N} \! \! dn \, \vec{r}(n,t) \phi_{p}(n),

and the inversion back to chain coordinates is written as

.. math::
    \vec{r}(n,t) = \sum_{p=0}^{\infty} \vec{X}_{p}(t) \phi_{p}(n).

Upon performing a transform to normal
coordinates, we find the governing
equation of motion

.. math::
    \xi N \frac{d \vec{X}_{p}}{d t} =
    - k_{p} \vec{X}_{p}
    + \vec{f}_{B_{p}}
    :label: eom-normal

where :math:`k_{p} = \frac{3 \pi^{2} k_{B}T}{b^{2} N} p^{2}`.
The :math:`p`-mode Brownian force :math:`\vec{f}_{B_{p}}`
is given by

.. math::
    \vec{f}_{B_{p}} = \int_{0}^{N} \! \! dn \, \vec{f}_{B}(n,t) \phi_{p}(n)

and satisfies the fluctuation dissipation theorem

.. math::
    \langle
    \vec{f}_{B_{p}} (t) \vec{f}_{B_{p'}}(t')
    \rangle =
    2 k_{B}T \xi N \delta_{pp'} \delta (t - t') \mathbf{I}.

A similar derivation for a ring polymer results in a treatment that
is identical to the linear chain, but the normal modes for the
ring polymer are continuous across the ends [i.e. :math:`\vec{r}(n=0,t) = \vec{r}(n=N,t)`].
Specifically, the complete normal-mode set is separated into even and odd functions,
respectively
defined by

.. math::
    & & \phi_{p}^{(e)}(n)=
    \left\{ \begin{array}{c}
		\sqrt{2} \cos \left( \frac{2 p \pi n}{N} \right), \ \ p = 1,2,\ldots \\
		1, \ \ \ \ \ \ \ \ \ \  p = 0.
		\end{array}
	\right.

and

.. math::
    \phi_{p}^{(o)}(n)=
		\sqrt{2} \sin \left( \frac{2 \pi p  n}{N} \right), \ \ p = 1,2,\ldots

The even and odd normal modes satisfy the equation of motion defined in Eq. :eq:`eom-normal`
with
:math:`k_{p} = \frac{12 \pi^{2} k_{B}T}{b^{2} N} p^{2}`

Mean-square displacement (MSD) for linear polymers
**************************************************

The mean-square displacement of a segment of the polymer chain is define a

.. math::
    \mathrm{MSD} = \langle \left( \vec{r}(n, t) - \vec{r}(n, 0) \right)^{2} \rangle
    :label: msd

for the nth segment of the chain with total length :math:`N`.
We insert our normal-mode representation into Eq. :eq:`msd`,
resulting in the expression

.. math::
    \mathrm{MSD} =
    \sum_{p=0}^{\infty} \sum_{p'=0}^{\infty}
    \langle
    \left( \vec{X}_{p}(t) - \vec{X}_{p}(0) \right) \cdot
    \left( \vec{X}_{p'}(t) - \vec{X}_{p'}(0) \right) \rangle
    \phi_{p}(n) \phi_{p'}(n)

The equation of motion (Eq. :eq:`eom-normal`) can be used to determine the correlation function
:math:`\langle \vec{X}_{p}(t) \cdot \vec{X}_{p'}(0) \rangle` (detailed discussion is found in Ref. [Doi1988]_).
This results in the expression

.. math::
    \langle \vec{X}_{p}(t) \cdot \vec{X}_{p'}(0) \rangle = 3 \frac{k_{B}T}{k_{p}}
    \exp \! \left( - \frac{k_{p}}{N \xi} t \right) \delta_{pp'}

for :math:`p \ge 1` and

.. math::
    \langle
    \left( \vec{X}_{0}(t) - \vec{X}_{0}(0) \right)^{2} \rangle = 6 \frac{k_{B}T}{N \xi} t

We focus on the :math:`\mathrm{MSD}` for the midpoint of a linear chain,
thus :math:`n=N/2`.
Inserting this into our definition of :math:`\mathrm{MSD}` results in the
expression for the MSD of the midpoint of a linear chain

.. math::
    \mathrm{MSD}
    & = &
    6 \frac{k_{B}T}{\xi N} t +
    \sum_{p \, \mathrm{even}} 12 \frac{k_{B}T}{k_{p}}
    \left[ 1 - \exp \! \left( - \frac{k_{p}}{N \xi} t \right) \right] \\
    & = &
    6 \frac{k_{B}T}{\xi N} t +
    \sum_{p = 1}^{\infty} 12 \frac{k_{B}T}{k_{2p}}
    \left[ 1 - \exp \! \left( - \frac{k_{2p}}{N \xi} t \right) \right]

Mean-squared change in distance (MSCD) for linear and ring polymers
*******************************************************************

We now consider the mean-square change in distance (MSCD) for a linear polymer
chain.
This quantity is defined as

.. math::
    \mathrm{MSCD} = \langle \left( \Delta \vec{R}(t) - \Delta \vec{R}(0) \right)^{2} \rangle
    :label: mscd

where :math:`\Delta \vec{R}(t) = \vec{r}(N/2 + \Delta, t) - \vec{r}(N/2 - \Delta,t)` where
the total chain length is :math:`N = 2 N_{s}`.
We insert our normal-mode representation into Eq. :eq:`mscd`
to find

.. math::
    \mathrm{MSCD} =
    \sum_{p=1}^{\infty} \sum_{p'=1}^{\infty}
    \langle
    \left( \vec{X}_{p}(t) - \vec{X}_{p}(0) \right) \cdot
    \left( \vec{X}_{p'}(t) - \vec{X}_{p'}(0) \right)
    \rangle
    \left[ \phi_{p}(N/2 + \Delta) - \phi_{p}(N/2 - \Delta) \right]
    \left[ \phi_{p'}(N/2 + \Delta) - \phi_{p'}(N/2 - \Delta) \right]

The equation of motion (Eq. :eq:`eom-normal`) can be used to determine the correlation function
:math:`\langle \vec{X}_{p}(t) \cdot \vec{X}_{p'}(0) \rangle` (detailed discussion is found in Ref. [Doi1988]_).
This results in the expression

.. math::
    \langle \vec{X}_{p}(t) \cdot \vec{X}_{p'}(0) \rangle = 3 \frac{k_{B}T}{k_{p}}
    \exp \! \left( - \frac{k_{p}}{N \xi} t \right) \delta_{pp'}

Inserting this into our definition of :math:`\mathrm{MSCD}` results in the
expression for the linear chain

.. math::
    \mathrm{MSCD}^{\mathrm{(linear)}}
    & = & \sum_{p \, \mathrm{odd}} 48 \frac{k_{B}T}{k_{p}}
    \left[ 1 - \exp \! \left( - \frac{k_{p}}{N \xi} t \right) \right]
    \sin^{2} \left( \frac{\pi p \Delta}{N} \right) \\
    & = & \sum_{p = 0}^{\infty} 48 \frac{k_{B}T}{k_{2p+1}}
    \left[ 1 - \exp \! \left( - \frac{k_{2p+1}}{N \xi} t \right) \right]
    \sin^{2} \left[ \frac{\pi (2p+1) \Delta}{N} \right]

where :math:`k_{p} = \frac{3 \pi^{2} k_{B}T}{b^{2} N} p^{2}` for the linear chain.

We follow a parallel derivation for the ring polymer.
We note that only the odd set of normal modes contribute to :math:`\mathrm{MSCD}`
for the ring polymer.
Taking similar steps as in the linear case,
we arrive at the expression for the ring polymer

.. math::
    \mathrm{MSCD}^{\mathrm{(ring)}}
    = \sum_{p = 1}^{\infty} 48 \frac{k_{B}T}{k_{p}}
    \left[ 1 - \exp \! \left( - \frac{k_{p}}{N \xi} t \right) \right]
    \sin^{2} \left( \frac{2 \pi p \Delta}{N} \right)

where :math:`k_{p} = \frac{12 \pi^{2} k_{B}T}{b^{2} N} p^{2}` for the ring polymer.

Functions contained with the 'poly_dyn' module
----------------------------------------------

.. automodule:: wlcstat.poly_dyn
    :members:


Example usage of 'linear_mscd' and 'ring_mscd'
----------------------------------------------

We show the solution for the MSCD for chain length :math:`N=100` and :math:`\Delta = 25` for both
linear and ring polymers. In this plot, we also show the MSD (using 'linear_mid_msd') multiplied by
2, which coincides with the short-time asymptotic behavior.
The long-time asymptotic behavior is found by noting

.. math::
    \mathrm{MSCD} = \langle \left( \Delta \vec{R}(t) - \Delta \vec{R}(0) \right)^{2} \rangle =
    \langle \Delta \vec{R}(t)^{2} \rangle
    - 2 \langle \Delta \vec{R}(t) \cdot \Delta \vec{R}(0) \rangle
    + \langle \Delta \vec{R}(0) ^{2} \rangle \rightarrow
    2 \langle \Delta \vec{R}^{2} \rangle

which leads to the asymptotic solutions

.. math::
    \mathrm{MSCD}^{\mathrm{(linear)}} \rightarrow 4 \Delta

and

.. math::
    \mathrm{MSCD}^{\mathrm{(ring)}} \rightarrow \frac{2}{\frac{1}{2 \Delta} + \frac{1}{N-2 \Delta}}.

These asymptotic solutions are also included in the plot.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.poly_dyn import *

    num_t = 100
    t_val_0 = 1e-4
    t_val_f = 1e6
    t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
    N=100
    Ndel=25
    mscdl = linear_mscd(t_val, 1, Ndel, N, 1, 20000)
    mscdr = ring_mscd(t_val, 1, Ndel, N, 1, 20000)
    msd = linear_mid_msd(t_val, 1, N, 1, 10000)
    mscdl_inf = 2 * 2 * Ndel
    mscdr_inf = 2 * 1 / (1 / (2 * Ndel) + 1  / (N - 2 * Ndel))
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    plt.loglog(t_val, mscdl)
    plt.loglog(t_val, mscdr)
    plt.loglog(t_val, 2 * msd)
    plt.loglog(t_val, 2 * 6 * t_val /N,'--')
    plt.loglog(t_val, mscdl_inf + 0*t_val,'--')
    plt.loglog(t_val, mscdr_inf + 0*t_val,'--')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$MSCD$')
    plt.tight_layout()
    plt.show()


Application of 'linear_mscd' and 'ring_mscd' to homologue pairing in meiosis
----------------------------------------------------------------------------

We apply 'linear_mscd' and 'ring_mscd' to the coordinated dynamics of loci on two homologous chromosomes
(chromosome V in S. cerevisiae). The average number of linkages between the chromosomes is :math:`\mu = 4`.
The length of chromosome V is approximately 577 kb, and the position of the fluorescent locus is at 177 kb.
The top plot below shows the locations of the random linkages for the five example "cells".
The effective tethers holding the loci together are highlighted in white.
The nearest linkage is highlighted by a thicker line. Notice that when two tethers exist,
they effectively form a large chromatin loop on which the two loci live.
The bottom plot shows analytical MSCDs for five different theoretical
"cells", where a Poisson-distributed number of linkages (:math:`\mu = 4`)
have been distributed uniformly along chr.~V, which we assume to be composed of
approximately 116,000 Kuhn lengths of length.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.poly_dyn import *

    mu = 4
    cells = [generate_example_cell(mu) for i in range(5)]
    ax, all_closest_links = draw_cells(cells)
    plt.show()
    t = np.logspace(1,4,50).astype(float)
    plt.figure()
    for i, linkages in enumerate(cells):
        plt.loglog(t, model_mscd(t,linkages), label = 'Cell ' + str(i+1))
    plt.legend()
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    plt.xlabel(r'time (sec)')
    plt.ylabel(r'$MSCD$ ($\mu m^{2}$)')
    plt.show()