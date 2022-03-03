.. _wlcstruc_notes:
.. .. automodule:: wlcstat


Structure Factor
================

The multi-scale statistical behavior of a polymer chain is frequently characterized by the
single-chain structure factor.
The structure factor is directly measured in a
scattering experiment, and correlating scattering experiments with theoretical
models provides insight into the physical behavior in
polymeric fluids.
From a theoretical perspective, the single-chain structure factor acts as an
input for treating many-chain systems including
concentration fluctuations.
Although the structure factors for both a flexible Gaussian chain
and a rigid rod are easily expressible, the structure factor of a semiflexible
polymer is considerably more difficult to predict theoretically,
largely due to the challenges in acquiring analytical expressions for the
governing chain statistics.
Our goal in this section is to apply our exact results to predict the
structure factor for the wormlike chain model that is valid over the
entire range of scattering vectors and chain lengths.

Two-point Structure Factor
--------------------------

We define the structure factor to be

.. math::
    S(\vec{k}) = \frac{1}{L^{2}}
    \int_{0}^{L} ds_{1}
    \int_{0}^{L} ds_{2}
    \left< \exp \left[ i \vec{k} \cdot \left( \vec{r}(s_{1})-\vec{r}(s_{2}) \right) \right] \right>,

.. \label{eq:sk}

where :math:`\vec{k}` is the scattering vector, and
the angular brackets denote an average with respect to the
specific model that governs the chain statistics.
Our current manuscript contains exact results that permit the evaluation of
the structure factor for the wormlike chain model without excluded volume
interactions.
Our definition for the structure factor differs from the conventional definition
by an extra factor of :math:`1/L`, thus rendering a dimensionless
structure factor that tends to one at zero scattering vector.

For the wormlike chain model, the structure factor is given by

.. math::
    S(\vec{k}) = \frac{2}{N^{2}} \mathcal{L}^{-1} \left[ \frac{G(K;p)}{p^{2}} \right],

.. \label{eq:sklaplace}

where :math:`K = 2 l_{p} | \vec{k} |`, :math:`\mathcal{L}^{-1}` implies a Laplace
inversion from :math:`p` to :math:`N=L/(2l_{p})`, and :math:`G(K;p)` is the wormlike chain
Green function.
The structure factor is expressible in terms of the magnitude of the scattering vector
due to the rotational invariance of the governing statistics.
Using results given in :ref:`Wormlike Chain Green's Function`,
we write the structure factor as

.. math::
    S(\vec{k}) =
    \frac{2}{N^{2}} \left[
    \frac{N}{j_{0}^{(+)}(K;0)} -
    \frac{\partial_{p} j_{0}^{(+)}(K;0)}{j_{0}^{(+)}(K;0)^{2}}
    + \sum_{\alpha=0}^{\infty}
    \frac{\exp \left( \mathcal{E}_{\alpha} N \right)}{\mathcal{E}_{\alpha}^{2} \partial_{p} j_{0}^{(+)}(K;\mathcal{E}_{\alpha})}
    \right],
    :label: skwlc

.. \label{eq:skwlc}

where :math:`j_{0}^{(+)}(K;p)` and :math:`\partial_{p} j_{0}^{(+)}(K;p)` are
continued-fraction functions that are defined in
:ref:`Wormlike Chain Green's Function`.
Equation :eq:`skwlc` is evaluated using methods outlined in
:ref:`Wormlike Chain Green's Function`.
to find the eigenvalues :math:`\mathcal{E}_{l}` and to calculate
:math:`j_{0}^{(+)}` and :math:`\partial_{p} j_{0}^{(+)}`.
The structure factor given by Eq. :eq:`skwlc` is expressed in arbitrary number of
dimensions :math:`D`, which is generally useful for treating many-chain systems.

In Ref. [Spakowitz2004]_, we present realizations of the
structure factor for the wormlike chain model in 3 dimensions,
demonstrating that our exact results
for wormlike chain statistics in 3 dimensions capture the structure factor
over all chain lengths and scattering vectors.
Equation :eq:`skwlc`, along with the techniques provided in
:ref:`Wormlike Chain Green's Function`,
provide a convenient methodology for calculating the structure factor
as a comparison with scattering experiments.
Specifically, the infinite summation in Eq. :eq:`skwlc` is truncated at a finite cutoff,
and the partial summation is evaluated using methods to evaluate
:math:`\mathcal{E}_{l}` found in :ref:`Wormlike Chain Green's Function`.
For most calculations, only a couple of terms are necessary to achieve
accurate realizations of the structure factor.
For example, the structure factor in 3 dimensions for a chain
of length :math:`N=0.1` requires only 4 terms in the summation to achieve
accuracy within 1.76 percent difference from a more accurate calculation with
30 terms in the summation (\emph{i.e.} completely converged in its numerical behavior).
Since the terms in the summation decay exponentially with :math:`N`, larger :math:`N` values
require fewer terms.
To illustrate this, we note that the structure factor in 3 dimensions for a slightly
larger chain of length :math:`N=0.5` calculated with 4 terms in the summation
is within :math:`1.05 \times 10^{-4}` percent difference from the calculation
with 30 terms in the summation.
Therefore, most practical applications of Eq. :eq:`skwlc` require the inclusion of
only a couple of terms in the summation.

Three-point and Four-point Density Correlation Functions
--------------------------------------------------------

We now extend the definition of the structure factor to include 3-point and 4-point
correlation functions, as these serve as necessary input for
a density expansion of the free energy
under the random phase approximation (RPA) to a field-theoretic treatment of polymer solutions.
The structure factors form the basis of the cubic and quartic vertex functions in the free energy
expansion due to the connection between the structure factor and correlated density fluctuations
in the solution.
We discuss the evaluation of these functions below.

Three-point structure factor
****************************

We define the 3-point structure factor :math:`S^{(3)}(\vec{k}_{1}, \vec{k}_{2})` as

.. math::
    S^{(3)}(\vec{k}_{1}, \vec{k}_{2}; L) =
    \frac{1}{L^{3}}
    \int_{0}^{L} \! \! ds_{1}
    \int_{0}^{L} \! \! ds_{2}
    \int_{0}^{L} \! \! ds_{3}
    \left< \exp \left[
    i \vec{k}_{1} \cdot \vec{r}(s_{1}) +
    i \vec{k}_{2} \cdot \vec{r}(s_{2}) -
    i \left( \vec{k}_{1} + \vec{k}_{2} \right) \cdot \vec{r}(s_{3})
    \right] \right>.

We note that the 3-point density corelations depend on 3 k-vectors, which satisfy
:math:`\vec{k}_{1}+\vec{k}_{2}+\vec{k}_{3}=0` or :math:`\vec{k}_{3} = -\vec{k}_{1}-\vec{k}_{2}`.
For the wormlike chain model,
we leverage the Green's function to evaluate this average.
For this treatment, we specialize the discussion to 3 dimensions,
since the evaluation of the average in arbitrary dimensions requires
rotation operations that are concretely defined in 3 dimensions.
The 3-point structure factor is
written in terms of the Green's function to be

.. math::
    S^{(3)}(\vec{k}_{1}, \vec{k}_{2}; L) & = &
    \frac{2}{N^{3}}
    \int_{0}^{N} \! \! ds_{3}
    \int_{0}^{s_{3}} \! \! ds_{2}
    \int_{0}^{s_{2}} \! \! ds_{1}
    \int \! \! d \vec{u}_{3}
    \int \! \! d \vec{u}_{2}
    \int \! \! d \vec{u}_{1}
    G(\vec{K}_{1} + \vec{K}_{2}, \vec{u}_{3} | \vec{u}_{2}; s_{3} - s_{2})
    G(\vec{K}_{1}, \vec{u}_{2} | \vec{u}_{1}; s_{2} - s_{1})
    \\
    &  &
    + \frac{2}{N^{3}}
    \int_{0}^{N} \! \! ds_{3}
    \int_{0}^{s_{3}} \! \! ds_{2}
    \int_{0}^{s_{2}} \! \! ds_{1}
    \int \! \! d \vec{u}_{3}
    \int \! \! d \vec{u}_{2}
    \int \! \! d \vec{u}_{1}
    G(\vec{K}_{1} + \vec{K}_{2}, \vec{u}_{3} | \vec{u}_{2}; s_{3} - s_{2})
    G(\vec{K}_{2}, \vec{u}_{2} | \vec{u}_{1}; s_{2} - s_{1})
    \\
    &  &
    + \frac{2}{N^{3}}
    \int_{0}^{N} \! \! ds_{3}
    \int_{0}^{s_{3}} \! \! ds_{2}
    \int_{0}^{s_{2}} \! \! ds_{1}
    \int \! \! d \vec{u}_{3}
    \int \! \! d \vec{u}_{2}
    \int \! \! d \vec{u}_{1}
    G(\vec{K}_{1}, \vec{u}_{3} | \vec{u}_{2}; s_{3} - s_{2})
    G(-\vec{K}_{2}, \vec{u}_{2} | \vec{u}_{1}; s_{2} - s_{1})

where :math:`\vec{K}_{i} = 2 l_{p} \vec{k}_{i}`.
In this calculation, we need to evaluate the rotated spherical harmonics (see :ref:`Rotation of spherical harmonics`).

We develop two methods to evaluate these integrals, which are implemented at small-:math:`k` values and
large-:math:`k` values (where :math:`k = |\vec{k}|`.  We describe these methods below.

Functions contained with the 'wlcstruc' module
----------------------------------------------

.. automodule:: wlcstat.wlcstruc
    :members:

Example usage of 's2_wlc'
------------------------------------------

This plot shows the structure factor for the wormlike chain in 3 dimensions over a range of :math:`K`
for :math:`N = 0.1`, :math:`0.5`, and :math:`1.0`.
This plot is a reproduction of Fig. 2 in Ref. [Spakowitz2004]_.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.wlcstruc import *

    num_k = 100
    k_val_0 = 1e-2
    k_val_f = 20
    length_kuhn_vec = np.array([0.1, 0.5, 1])
    dimensions = 3

    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)

    for ind_length in range(0, len(length_kuhn_vec)):
        length_kuhn = float(length_kuhn_vec[ind_length])
        k_val = np.linspace(k_val_0, k_val_f / length_kuhn, num_k)
        structure_factor = s2_wlc(k_val, length_kuhn, dimensions)
        plt.plot(k_val * length_kuhn, np.real(k_val * structure_factor[:, 0] * length_kuhn), '-')

    plt.xlabel(r'$Lk$')
    plt.ylabel(r'Structure Factor, $k L S(K;N)$')
    plt.ylim((2, 3.5))
    plt.xlim((2, 20))
    plt.tight_layout()
    plt.show()

