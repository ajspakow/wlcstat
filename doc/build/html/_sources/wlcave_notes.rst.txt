.. _wlcave_notes:
.. .. automodule:: wlcstat


Average Quantities
==================

The orientation Green function can be used to evaluate any average quantity that can be
expressed in terms of the tangent orientations.
The Green function operates as a conditional probability, such that :math:`G(\vec{u}|\vec{u}_{0};L)` gives the 
probability that a chain of length :math:`L` with an initial tangent orientation :math:`\vec{u}_{0}` will have a final
tangent orientation :math:`\vec{u}`.
We consider an arbitrary orientation-only average quantity :math:`A^{n} =
\langle A_{n}(s_{n}) A_{n-1}(s_{n-1}) \ldots A_{1}(s_{1}) \rangle`,
where :math:`A_{m}(s_{m})` is a function that is expressible only in terms of the tangent orientation at arclength
position :math:`s_{m}`.
The arclength positions are ordered such that :math:`s_{1} < s_{2} < \ldots < s_{n-1} < s_{n}`.
The average quantity :math:`A^{n}` is evaluated using the following procedure:

- An orientation Green function :math:`G(\vec{u}_{1}|\vec{u}_{0};s_{1})` is inserted into :math:`A^{n}` to propagate from the initial orientation :math:`\vec{u}_{0}` at :math:`s=0` to the first orientation :math:`\vec{u}_{1}` at :math:`s_{1}`

- An orientation Green function :math:`G(\vec{u}_{m+1}|\vec{u}_{m};s_{m+1}-s_{m})` is inserted into :math:`A^{n}` between :math:`A_{m+1}` and :math:`A_{m}`, resulting in a product of :math:`n-1` Green functions that propagate between the internal orientations in the average

- An orientation Green function :math:`G(\vec{u}|\vec{u}_{n};L-s_{n})` is inserted into :math:`A^{n}` to propagate from the :math:`n` th orientation :math:`\vec{u}_{n}` at :math:`s_{n}` to the final orientation :math:`\vec{u}` at :math:`L`

- The average :math:`A^{n}` is normalized by a factor that represents the orientational integral over the full-chain Green function, i.e. the same quantity as the average without the product of :math:`A_{m}`

- The average is then evaluated by integrating over the orientations :math:`\vec{u}_{0}`, :math:`\vec{u}_{1}`, \ldots, :math:`\vec{u}_{n}`, and :math:`\vec{u}` 

This procedure is amenable to the evaluation of arbitrary average quantities that are expressible in terms of the
tangent orientations :math:`\vec{u}_{m}` along the chain at positions :math:`s_{m}`.

For example, the average :math:`\langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle` (where :math:`s_{1} < s_{2}`) is given by

.. math::
    \langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle & = &
    \frac{1}{\mathcal{N}}
    \int d \vec{u} d \vec{u}_{2} d \vec{u}_{1} d \vec{u}_{0}
    G(\vec{u}|\vec{u}_{2};L-s_{2})
    G(\vec{u}_{2}|\vec{u}_{1};s_{2}-s_{1})
    G(\vec{u}_{1}|\vec{u}_{0};s_{1})
    \vec{u}_{2} \cdot \vec{u}_{1}
    \nonumber \\
    & = &
    \frac{3}{\mathcal{N}}
    \int d \vec{u} d \vec{u}_{2} d \vec{u}_{1} d \vec{u}_{0}
    G(\vec{u}|\vec{u}_{2};L-s_{2})
    G(\vec{u}_{2}|\vec{u}_{1};s_{2}-s_{1})
    G(\vec{u}_{1}|\vec{u}_{0};s_{1})
    u_{2}^{(z)} u_{1}^{(z)},
    :label: u1u2
..   \label{eq:u1u2}

where rotational symmetry implies
:math:`\langle u^{(x)}(s_{2}) u^{(x)}(s_{1}) \rangle =\langle u^{(y)}(s_{2}) u^{(y)}(s_{1}) \rangle = \langle u^{(z)}(s_{2}) u^{(z)}(s_{1}) \rangle` and, consequently, :math:`\langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle = 3 \langle u^{(z)}(s_{2}) u^{(z)}(s_{1}) \rangle`.
The normalization factor :math:`\mathcal{N}` is given by

.. math::
   \mathcal{N} = \int d \vec{u} d \vec{u}_{0} G(\vec{u}|\vec{u}_{0};L) = 4 \pi.

Inserting :eq:`guwlc` into :eq:`u1u2` gives

.. math::
   \langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle & = &
   \frac{3}{4 \pi} \int d \vec{u} d \vec{u}_{2} d \vec{u}_{1} d
   \vec{u}_{0} \sum_{l_{2}=0}^{\infty} \sum_{m_{2}=-l_{2}}^{l_{2}}
   \sum_{l_{1}=0}^{\infty} \sum_{m_{1}=-l_{1}}^{l_{1}}
   \sum_{l_{0}=0}^{\infty} \sum_{m_{0}=-l_{0}}^{l_{0}} \times
   \nonumber \\ &  & Y_{l_{2}}^{m_{2}} (\vec{u})
   Y_{l_{2}}^{m_{2}*}(\vec{u}_{2})  Y_{l_{1}}^{m_{1}} (\vec{u}_{2})
   Y_{l_{1}}^{m_{1}*}(\vec{u}_{1})  Y_{l_{0}}^{m_{0}} (\vec{u}_{1})
   Y_{l_{0}}^{m_{0}*}(\vec{u}_{0})  \cos \theta_{2} \cos \theta_{1}
   \times \nonumber \\ &  & \exp \! \left[
   -l_{2}(l_{2}+1)\frac{L-s_{2}}{2l_{p}}
   -l_{1}(l_{1}+1)\frac{s_{2}-s_{1}}{2l_{p}}
   -l_{0}(l_{0}+1)\frac{s_{1}}{2l_{p}}  \right],

where :math:`u_{2}^{(z)}=\cos \theta_{2}`.
Using the properties of the spherical harmonics [Arfken1999]_,
we note that :math:`\cos \theta = 2 \sqrt{\pi/3} Y_{1}^{0}(\vec{u})`,
and the average is written as

.. math::
   \langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle & = & 
   \int d \vec{u} d \vec{u}_{2} d \vec{u}_{1} d \vec{u}_{0} 
    \sum_{l_{2}=0}^{\infty} \sum_{m_{2}=-l_{2}}^{l_{2}} 
    \sum_{l_{1}=0}^{\infty} \sum_{m_{1}=-l_{1}}^{l_{1}} 
    \sum_{l_{0}=0}^{\infty} \sum_{m_{0}=-l_{0}}^{l_{0}} \times
    \nonumber \\
    &  &
    Y_{l_{2}}^{m_{2}} (\vec{u})
    Y_{l_{2}}^{m_{2}*}(\vec{u}_{2}) 
    Y_{l_{1}}^{m_{1}} (\vec{u}_{2})
    Y_{l_{1}}^{m_{1}*}(\vec{u}_{1}) 
    Y_{l_{0}}^{m_{0}} (\vec{u}_{1})
    Y_{l_{0}}^{m_{0}*}(\vec{u}_{0}) 
    Y_{1}^{0*}(\vec{u}_{2})
    Y_{1}^{0}(\vec{u}_{1})
    \times
    \nonumber \\
    &  &
    \exp \! \left[ 
    -l_{2}(l_{2}+1)\frac{L-s_{2}}{2l_{p}} 
    -l_{1}(l_{1}+1)\frac{s_{2}-s_{1}}{2l_{p}} 
    -l_{0}(l_{0}+1)\frac{s_{1}}{2l_{p}} 
    \right]

Integration over :math:`\vec{u}` and :math:`\vec{u}_{0}` force :math:`l_{2}=m_{2}=l_{0}=m_{0}=0`,
and the subsequent integration over :math:`\vec{u}_{1}` and :math:`\vec{u}_{2}` result in 
:math:`l_{1}=1` and :math:`m_{1}=0`.
The final expression is given by

.. math::
     \langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle = \exp \! \left( 
    - \frac{s_{2} - s_{1}}{l_{p}}
    \right),

which demonstrates the role of the persistence length :math:`l_{p}` as a correlation length for 
the tangent orientation.

This average shows that the orientation statistics can be used to directly evaluate averages that are expressed in 
terms of the tangent vector :math:`\vec{u}`.
However, other average quantities can be evaluated using the orientation Green function if the quantity can be expressed in terms of the
tangent orientation.
For example, the mean-square end-to-end vector is written as

.. math::
    \langle
    \vec{R}^{2}
    \rangle & = & 
    \int_{0}^{L} \int_{0}^{L} ds_{2} ds_{1} \langle \vec{u}(s_{2}) \cdot \vec{u}(s_{1}) \rangle
    = 2 \int_{0}^{L} \int_{0}^{s_{2}} ds_{2} ds_{1} \exp \! \left( 
    - \frac{s_{2} - s_{1}}{l_{p}}
    \right)  \\
    & = &
    2l_{p}^{2} \left[
    2N - 1 + \exp \left( -2N \right)
    \right],

where :math:`N=L/(2 l_{p})`.

An alternative approach to calculating averages is to use the Fourier-transformed Green's
function as a generator of averages
(see Refs. [Spakowitz2004]_, [Spakowitz2005]_, [Spakowitz2006]_, and [Mehraeen2008]_.
Development of this approach is provided in our
discussion of the Green's function.

The 'wlcave.py' module provides scripts to calculate a number of
average quantities for the wormlike chain model. These include the
following:

- The mean-square end-to-end distance :math:`\langle R^{2} \rangle`

- The mean-square radius of gyration :math:`\langle \vec{R}_{G}^{2} \rangle`

- The 4th moment of the end-to-end distribution :math:`\langle R_{z}^{4} \rangle / (2 l_{p})^{4}`

Functions contained with the 'wlcave' module
--------------------------------------------

.. automodule:: wlcstat.wlcave
    :members:


Example usage of 'r2_ave'
-------------------------

This example gives the mean-square end-to-end distance :math:`\langle R^{2} \rangle / (2 l_{p})^{2}`
(i.e. length non-dimensionalized by :math:`2 l_{p}`) versus
chain length :math:`N = L/(2 l_{p})` for 3 dimensions.  The short-length asymptotic behavior
:math:`\langle R^{2} \rangle / (2 l_{p})^{2} \rightarrow N^{2}`
and long-length asymptotic behavior
:math:`\langle R^{2} \rangle / (2 l_{p})^{2} \rightarrow 2 N/(d-1)` are included to show the
limiting behaviors.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.wlcave import *

    num_pts = 100
    length_kuhn_0 = 1e-3
    length_kuhn_f = 1e3
    length_kuhn = np.logspace(np.log10(length_kuhn_0), np.log10(length_kuhn_f), num_pts)
    dimensions = 3
    r2 = r2_ave(length_kuhn, dimensions)
    r2short = length_kuhn ** 2
    r2long = length_kuhn * 2 / (dimensions - 1)
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    plt.loglog(length_kuhn, r2)
    plt.loglog(length_kuhn[0:60], r2short[0:60])        # Short length asymptotic solution
    plt.loglog(length_kuhn[40:100], r2long[40:100]) # Long length asymptotic solution
    plt.xlabel(r'Chain length $N = L/(2l_{p})$')
    plt.ylabel(r'End-to-end distance $\langle R^{2} \rangle$')
    plt.tight_layout()
    plt.show()

Example usage of 'rg2_ave'
--------------------------

This example gives the mean-square radius of gyration :math:`\langle \vec{R}_{G}^{2} \rangle / (2 l_{p})^{2}`
(i.e. length non-dimensionalized by :math:`2 l_{p}`) versus
chain length :math:`N = L/(2 l_{p})` for 3 dimensions.  The short-length asymptotic behavior
:math:`\langle \vec{R}_{G}^{2} \rangle / (2 l_{p})^{2} \rightarrow N^{2}/12`
and long-length asymptotic behavior
:math:`\langle \vec{R}_{G}^{2} \rangle / (2 l_{p})^{2} \rightarrow N/[3 (d-1)]` are included to show the
limiting behaviors.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.wlcave import *

    num_pts = 100
    length_kuhn_0 = 1e-3
    length_kuhn_f = 1e3
    length_kuhn = np.logspace(np.log10(length_kuhn_0), np.log10(length_kuhn_f), num_pts)
    dimensions = 3
    rg2 = rg2_ave(length_kuhn, dimensions)
    rg2short = length_kuhn ** 2 / 12
    rg2long = length_kuhn / (3 * (dimensions - 1))
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    plt.loglog(length_kuhn, rg2)
    plt.loglog(length_kuhn[0:60], rg2short[0:60])        # Short length asymptotic solution
    plt.loglog(length_kuhn[40:100], rg2long[40:100]) # Long length asymptotic solution
    plt.xlabel(r'Chain length $N = L/(2l_{p})$')
    plt.ylabel(r'Mean-square radius of gyration $\langle \vec{R}_{G}^{2} \rangle$')
    plt.tight_layout()
    plt.show()

Example usage of 'rz4_ave'
--------------------------

This example gives the 4th moment of the end-to-end distribution :math:`\langle R_{z}^{4} \rangle / (2 l_{p})^{4}`
(i.e. length non-dimensionalized by :math:`2 l_{p}`) versus
chain length :math:`N = L/(2 l_{p})` for 3 dimensions.  The short-length asymptotic behavior
:math:`\langle R_{z}^{4} \rangle / (2 l_{p})^{2} \rightarrow N^{4} 3 / [d (d + 2)]`
and long-length asymptotic behavior
:math:`\langle R_{z}^{4} \rangle / (2 l_{p})^{2} \rightarrow 12 N^2/[d^2 (d-1)^2]` are included to show the
limiting behaviors.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.wlcave import *

    num_pts = 100
    length_kuhn_0 = 1e-3
    length_kuhn_f = 1e3
    length_kuhn = np.logspace(np.log10(length_kuhn_0), np.log10(length_kuhn_f), num_pts)
    dimensions = 3
    rz4 = rz4_ave(length_kuhn, dimensions)
    rz4short = length_kuhn ** 4 * 3 / (dimensions * (dimensions + 2))
    rz4long = length_kuhn ** 2 * 12 / (dimensions * (dimensions - 1)) ** 2
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    plt.loglog(length_kuhn, rz4)
    plt.loglog(length_kuhn[0:60], rz4short[0:60])        # Short length asymptotic solution
    plt.loglog(length_kuhn[40:100], rz4long[40:100]) # Long length asymptotic solution
    plt.xlabel(r'Chain length $N = L/(2l_{p})$')
    plt.ylabel(r'4th moment $\langle R_{z}^{4} \rangle$')
    plt.tight_layout()
    plt.show()
