.. _wlc_lcpoly:
.. .. automodule:: wlc_lcpoly

Polymer Field Theory: Liquid Crystalline Solutions
==================================================

We define a system Hamiltonian [Spakowitz2003]_ of a collection of polymer chains and solvent molecules

.. math::
    \beta \mathcal{H} =
    \sum_{i=1}^{n_{p}} \frac{\epsilon}{2} \int_{0}^{L}
    \left( \frac{\partial \vec{u}_{i}}{\partial s} \right)^{\! \! 2} d s +
    \chi \int d \vec{r} \, \hat{\phi}_{s}(\vec{r} \,) \hat{\phi}_{p}(\vec{r} \,) -
    \frac{a}{2} \int d \vec{r} \, \hat{\mathbf{S}}(\vec{r} \,):\hat{\mathbf{S}}(\vec{r} \,)
    :label: sysH

:math:`\hat{\phi}_{s}` and :math:`\hat{\phi}_{p}` are the local volume fractions of the solvent and polymer molecules, and :math:`\hat{\mathbf{S}}`
is the tensorial nematic order parameter

.. math::
    \hat{\phi}_{s}(\vec{r} \,) & = & v \sum_{j=1}^{n_{s}}
    \delta (\vec{r}-\vec{r}_{j})
    \\
    \hat{\phi}_{p}(\vec{r} \,) & = & A \sum_{i=1}^{n_{p}}
    \int_{0}^{L} d s \delta (\vec{r}-\vec{r}_{i}(s))
    \\
    \hat{\mathbf{S}} (\vec{r} \,) & = & A \sum_{i=1}^{n_{p}}
    \int_{0}^{L} d s \delta (\vec{r}-\vec{r}_{i}(s)) %\nonumber \\
    \left( \vec{u}_{i}(s) \vec{u}_{i}(s) -
    \frac{1}{3} \mathbf{I} \right)

The Flory-Huggins parameter :math:`\chi` gives the strength of polymer-solvent
interactions.
The Maier-Saupe parameter :math:`a` gives the strength of aligning interactions
that lead to liquid crystal phase phenomena.

The thermodynamic behavior is found by finding the grand canonical partition
function

.. math::
    \Xi & = & \exp [-\beta G(T,V,\mu)] \nonumber \\
    & = &\sum_{n_{s},n_{p}=0}^{\infty} \frac{1}{n_{s}! n_{p}!}
    \frac{e^{\beta \mu L A n_{p}}}{v^{n_{s}} (L A)^{n_{p}}}
    \int \prod_{j=1}^{n_{s}} d \vec{r}_{j} \int \prod_{i=1}^{n_{p}}
    \mathcal{D}[\vec{r}_{i}] \nonumber \\
    &   & \hspace{0.1in} \times \prod_{\vec{r}}
    \delta (\hat{\phi}_{s}+\hat{\phi}_{p}-1)
    \prod_{s}\delta (\vec{u}_{i} - \partial_{s} \vec{r}_{i})
    e^{-\beta \mathcal{H}[\hat{\phi}_{s},\hat{\phi}_{p},\hat{\mathbf{S}}]}
    :label: xidefine

where :math:`\prod_{\vec{r}}
\delta (\hat{\phi}_{s}+\hat{\phi}_{p}-1)` accounts for the incompressibility constraint
:math:`\hat{\phi}_{s}+\hat{\phi}_{p} = 1` at all locations,
and :math:`\prod_{s}\delta (\vec{u}_{i} - \partial_{s} \vec{r}_{i})` denotes the
fixed chain length constraint
:math:`\vec{u}_{i}=\partial_{s} \vec{r}_{i}`.
This cannot be solved exactly due to the many-body interactions that are implicit
within the Hamiltonian.

To make progress, we
use field-theoretical techniques~\cite{kn:fredreview,kn:matsen}
to transform the many-chain problem into a single-chain problem in fluctuating
effective potential fields.
The field-theoretic representation emerges from a series of
identity transformations that leave the theory unaffected.

First, we introduce collective variables
:math:`\phi_{s}`,
:math:`\phi_{p}`, and
:math:`\mathbf{S}` into the grand canonical partition function, giving

.. math::
    \Xi & = &
    \int \! \! \mathcal{D} \phi_{s}
    \int \! \! \mathcal{D} \phi_{p}
    \int \! \! \mathcal{D} \mathbf{S} \! \! \!
    \sum_{n_{s},n_{p}=0}^{\infty} \frac{1}{n_{s}! n_{p}!}
    \frac{e^{\beta \mu L A n_{p}}}{v^{n_{s}} (L A)^{n_{p}}}
    \int \prod_{j=1}^{n_{s}} d \vec{r}_{j} \nonumber \\
    &   & \hspace{0.1in} \times
    \int \prod_{i=1}^{n_{p}}
    \mathcal{D}[\vec{r}_{i}]
    \prod_{\vec{r}}
    \delta (\phi_{s} - \hat{\phi}_{s})
    \prod_{\vec{r}}
    \delta (\phi_{p} - \hat{\phi}_{p})
    \prod_{\vec{r}}
    \delta (\mathbf{S} - \hat{\mathbf{S}})
    \nonumber \\
    &   & \hspace{0.1in} \times \prod_{\vec{r}}
    \delta (\phi_{s}+\phi_{p}-1)
    \prod_{s}\delta (\vec{u}_{i} - \partial_{s} \vec{r}_{i})
    e^{-\beta \mathcal{H}[\phi_{s},\phi_{p},\mathbf{S}]}

Second, rewrite the delta functions using Fourier representations

.. math::
    \prod_{\vec{r}}
    \delta (\phi_{s} - \hat{\phi}_{s}) =
    \int \! \! \mathcal{D} W_{s}
    \exp \left[
    i \int \! \! d \vec{r} W_{s}(\vec{r})
    \left(
    \phi_{s}(\vec{r}) - \hat{\phi}_{s}(\vec{r})
    \right)
    \right] \\
    \prod_{\vec{r}}
    \delta (\phi_{p} - \hat{\phi}_{p}) =
    \int \! \! \mathcal{D} W_{p}
    \exp \left[
    i \int \! \! d \vec{r} W_{p}(\vec{r})
    \left(
    \phi_{p}(\vec{r}) - \hat{\phi}_{p}(\vec{r})
    \right)
    \right] \\
    \prod_{\vec{r}}
    \delta (\mathbf{S} - \hat{\mathbf{S}}) =
    \int \! \! \mathcal{D} \mathbf{\lambda}
    \exp \left[
    i \int \! \! d \vec{r} \mathbf{\lambda}(\vec{r}) :
    \left(
    \mathbf{S}(\vec{r}) - \hat{\mathbf{S}}(\vec{r})
    \right)
    \right]

The grand canonical partition function is now written as

.. math::
    \Xi & = & \int \mathcal{D} \phi_{p} \mathcal{D} W_{s}
    \mathcal{D} W_{p} \mathcal{D} \mathbf{S} \mathcal{D} \mathbf{\lambda}
    \exp \left\{ i \int d \vec{r} \,
    W_{s} (1-\phi_{p}) + \right. \nonumber \\
    &   &  i  \int d \vec{r} \, \left[ W_{p} \phi_{p} + \mathbf{\lambda}:\mathbf{S} \right] -
    \chi   \int d \vec{r} \,  \phi_{p} (1-\phi_{p}) + \nonumber \\
    &   & \left. \frac{a}{2} \int d \vec{r} \, \mathbf{S}:\mathbf{S}
    + \frac{1}{v} z_{s}(W_{s}) +
    \frac{e^{\beta \mu L A}}{L A} z_{p}(W_{p}, \mathbf{\lambda}) \right\}
    :label: grandcanon

Here, we have defined the
single-solvent-molecule partition function and single-polymer-chain
partition function, given by

.. math::
    z_{s}(W_{s}) = \int d \vec{r} \, e^{-i v W_{s}(\vec{r} \,)}
    :label: zs

and

.. math::
    &   & z_{p}(W_{p},\mathbf{\lambda}) = \int \mathcal{D}[\vec{r} \,(s)]
    \exp \left\{- \frac{\epsilon}{2} \int_{0}^{L}
    \left( \frac{\partial \vec{u}}{\partial s} \right)^{\! \! 2} d s -
    \right. \nonumber \\
    &   &  i A \int_{0}^{L} d s \left[ W_{p}(\vec{r} \,(s)) +
    \left. \mathbf{\lambda}(\vec{r} \,(s)) \! : \! \!
    \left(\vec{u} \vec{u} - \frac{1}{3} \mathbf{I}\right) \right] \right\}
    :label: zp

This form still cannot be solved exactly, but we can find a
saddle-point approximation that essentially gives the maximum
term in the functional integrals.
To do this, we write the free energy functional

.. math::
    &  &
    - \beta G =
    i \int d \vec{r} \,
    W_{s} (1-\phi_{p}) +
    i  \int d \vec{r} \, \left[ W_{p} \phi_{p} + \mathbf{\lambda}:\mathbf{S} \right]  \\
    &  &
    -  \chi   \int d \vec{r} \,  \phi_{p} (1-\phi_{p}) +
    \frac{a}{2} \int d \vec{r} \, \mathbf{S}:\mathbf{S}
    + \frac{1}{v} z_{s}(W_{s}) +
    \frac{e^{\beta \mu L A}}{L A} z_{p}(W_{p}, \mathbf{\lambda}) \nonumber

which will be evaluated at values of
:math:`\phi_{s}`,
:math:`\phi_{p}`,
:math:`\mathbf{S}`,
:math:`W_{s}`,
:math:`W_{p}`, and
:math:`\mathbf{\lambda}`
that maximize :math:`-\beta G` (or free energy minimum).

We find the first variation of the free energy to be

.. math::
    &  &
    - \delta \beta G =
    - i \int \! \! d \vec{r} W_{s} \delta \phi_{p}
    +i \int \! \! d \vec{r} (1- \phi_{p}) \delta W_{s}
    + i \int \! \! d \vec{r} W_{p} \delta \phi_{p}
    + i \int \! \! d \vec{r} \phi_{p} \delta W_{p} \nonumber \\
    &  &
    + i \int \! \! d \vec{r} \mathbf{\lambda} : \delta \mathbf{S}
    + i \int \! \! d \vec{r} \mathbf{S} : \delta \mathbf{\lambda}
    - \chi
    \int \! \! d \vec{r} (1 - 2 \phi_{p}) \delta \phi_{p}
    + a \int \! \! d \vec{r} \mathbf{S} : \delta \mathbf{S} \nonumber \\
    &  &
    +\frac{1}{v} \delta z_{s}(W_{s}) +
    \frac{e^{\beta \mu L A}}{L A} \delta z_{p}(W_{p}, \mathbf{\lambda})


where

.. math::
    &  &
    \delta z_{s}(W_{s}) = - i v
    \int \! \! d \vec{r} \exp
    \left[
    - i v W_{s}(\vec{r})
    \right] \delta W_{s} \\
    &  &
    \delta z_{p}(W_{p}, \mathbf{\lambda}) =
    \int \! \! d \vec{r}
    \int \mathcal{D}[\vec{r} \,(s)]
    \int_{0}^{L} \! \! ds
    \left[
    -iA \delta (\vec{r} - \vec{r}(s)) \delta W_{p}(\vec{r})
    - i A \delta (\vec{r} - \vec{r}(s))
    \left(\vec{u} \vec{u} - \frac{1}{3} \mathbf{I}\right) :
    \delta \mathbf{\lambda}
    \right] \nonumber \\
    &  &
    \exp \left\{- \frac{\epsilon}{2} \int_{0}^{L}
    \left( \frac{\partial \vec{u}}{\partial s} \right)^{\! \! 2} d s -
    i A \int_{0}^{L} d s \left[ W_{p}(\vec{r} \,(s)) +
    \mathbf{\lambda}(\vec{r} \,(s)) \! : \! \!
    \left(\vec{u} \vec{u} - \frac{1}{3} \mathbf{I}\right) \right] \right\}

Setting :math:`\delta \beta G` to zero gives the self-consistent field equations:

.. math::
    1-\phi_{p} & = & e^{-v w_{s}},
    \\
    \phi_{p} & = & - \frac{e^{\beta \mu L A}}{L A}
    \frac{\delta z_{p}}{\delta w_{p}},
    \\
    w_{p} - w_{s} & = & \chi (1 - 2 \phi_{p}),
    \\
    \mathbf{h} & = & - a \mathbf{S}, \label{saddles} \\
    \mathbf{S} & = & - \frac{ e^{\beta \mu L A}}{L A}
    \frac{\delta z_{p}}{\delta \mathbf{h}},

where :math:`w_{s} = i W_{s}`, :math:`w_{p} = i W_{p}`,
and :math:`\mathbf{h}=i \mathbf{\lambda}`.


Now consider the special case of having a homogeneous solution
:math:`\phi_{s}(\vec{r}) = \phi_{s}`,
:math:`\phi_{p}(\vec{r}) = \phi_{p}`, and nematic order parameter
:math:`\mathbf{S}(\vec{r}) = S_{0} \left( \hat{z}\hat{z} - \frac{1}{3} \mathbf{I} \right)`.
Define a
normalized scalar order parameter :math:`m` by

.. math::
    m = \frac{S_{0}}{\phi_{p}} = \frac{1}{L} \int_{0}^{L} d s
    \left< \left( \frac{3}{2} u_{z}^{2}(s) - \frac{1}{2} \right) \right>,
    :label: mdefine

where :math:`\left< \ldots \right>`
indicates an average with respect to the single-chain, self-consistent-field
Hamiltonian,

.. math::
    \beta \mathcal{H}_{0} = \frac{\epsilon}{2} \int_{0}^{L}
    \left( \frac{\partial \vec{u}}{\partial s} \right)^{2} d s -
    a \phi_{p} m A \int_{0}^{L} d s
    \left( u_{z}^{2}(s) - \frac{1}{3} \right).
    :label: hsaddle

The grand potential density
:math:`g = G/V` is equivalently expressed as the osmotic pressure defined through :math:`p=-[g(\phi_p)-g(0)]` and given by

.. math::
    \beta p & = & -\frac{1}{v} \log (1-\phi_{p}) - \chi \phi_{p}^{2}
    -\frac{1}{3} a \phi_{p}^{2} m^{2}  \nonumber \\
    &   & \hspace{0.5in} - \frac{\phi_{p}}{v}
    + \frac{\phi_{p}}{L A}
    :label: gsaddle

The chemical potential :math:`\mu` is given by

.. math::
    \beta \mu & = & \frac{1}{L A} \log \phi_{p} - \frac{1}{v} \log (1-\phi_{p})+ \nonumber \\
    &   & \hspace{0.5in} \chi (1-2 \phi_{p}) - \frac{1}{L A} \log q
    :label: musaddle

where the single-polymer molecule orientation partition function :math:`q`
is calculated from the self-consistent field Hamiltonian, Eq. :eq:`hsaddle`.
The Helmholtz free energy density (up to an additive constant) for a single, homogeneous phase is

.. math::
    \beta f & = &
    -\beta p + \beta \mu \phi_{p} \nonumber \\
    & = & \frac{\phi_{p}}{L A} \log \phi_{p} +
    \frac{(1-\phi_{p})}{v} \log (1-\phi_{p})
    +\chi \phi_{p} (1-\phi_{p}) \nonumber \\
    &   & + \frac{a}{3} \phi_{p}^{2} m^{2} +\frac{\phi_{p}}{v}
    - \frac{\phi_{p}}{L A} - \frac{\phi_{p}}{L A} \log q
    :label: fsaddle

This agrees with the Flory-Huggins theory when :math:`a=0` and :math:`q` is
a constant (with unimportant factors that are linear in :math:`\phi_{p}`).

Functions contained with the 'wlc_lcpoly' module
------------------------------------------------

.. automodule:: wlcstat.wlc_lcpoly
    :members:

Example usage of 'q_lcpoly'
---------------------------

We reproduce Fig. 6 from [Spakowitz2003]_ to demonstrate the use of 'q_lcpoly'.
We show free energy relative to the isotropic state :math:`\Delta f` versus
the order parameter :math:`m` for a thermotropic liquid-crystalline polymer system with
:math:`N = L/(2lp)=3.3333` and
:math:`2 l_{p} \kappa = 20.7582` (blue),
:math:`2 l_{p} \kappa = 21.0606` (orange), and
:math:`2 l_{p} \kappa = 23.6844` (green).

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.wlc_lcpoly import *

    lam_0 = 0
    lam_f = 25
    n_lam = 500
    lam_vec = np.linspace(lam_0, lam_f, n_lam)
    length_kuhn = 3.333
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    kappa_vec = np.array([20.7582, 21.0606, 23.6844])
    for i_kap in range(len(kappa_vec)):
        kappa = kappa_vec[i_kap]
        m_val = np.zeros(n_lam)
        q_val = np.zeros(n_lam)
        f_val = np.zeros(n_lam)
        for i_lam in range(n_lam):
            lam = lam_vec[i_lam]
            q_val[i_lam] = q_lcpoly(length_kuhn,lam)
            m_val[i_lam] = lam / kappa
            f_val[i_lam] = kappa * length_kuhn * m_val[i_lam] ** 2 / 3 - np.log(q_val[i_lam])
        plt.plot(m_val, f_val,'-')
    plt.ylabel(r'$\beta \Delta f$')
    plt.xlabel(r'$m$')
    plt.ylim((-0.6, 0.6))
    plt.xlim((0, 0.8))
    plt.tight_layout()
    plt.show()