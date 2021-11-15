.. _poly_confine:

Polymer Confinement
===================

We consider a flexible Gaussian chain in a spherical confinement of radius :math:`a`.
Within the confinement, the potential is set to zero, i.e. :math:`V = 0` for :math:`R < a`.
The Green function is zero outside the confinement, so we enforce the confinement by setting
:math:`G = 0` at :math:`R = a` and :math:`R_{0} = a`.
The Green function for the 3-D Gaussian Chain within the confinement obeys

.. math::
    \frac{\partial G(\vec{R},N|\vec{R}_{0},0)}{\partial N} =
    \frac{b^{2}}{6}
    \vec{\nabla}^{2} G(\vec{R},N|\vec{R}_{0},0)

with the boundary condition

.. math::
    G(\vec{R} = a \hat{e}_{R},N|\vec{R}_{0},0) & = & 0 \\
    G(\vec{R},N|\vec{R}_{0} = a \hat{e}_{R_{0}},0) & = & 0

where :math:`\hat{e}_{R}` and :math:`\hat{e}_{R_{0}}` are unit vectors in the
direction of :math:`\vec{R}` and :math:`\vec{R}_{0}`, respectively.

We use separation of variables to solve for the Green function within the confinement.
We define the position-dependent wavefunction :math:`\psi_{nlm} = R_{nl}(R) Y_{l}^{m}(\hat{e})`.
We non-dimensionalize the radial distance by :math:`a`, defining the dimensionless position vector
:math:`\vec{r} = \vec{R}/a`.
With this definition, we determine the Green function by solving the appropriate
eigenfunction problem, given by

.. math::
    \mathcal{H} \psi_{nlm} =
    \vec{\nabla}^{2} \psi_{nlm} =
    - \lambda_{nlm}^{2} \psi_{nlm}

Inserting our wavefunction into this equation gives the governing equation

.. math::
    \left[
    \frac{1}{r^{2}} \frac{d }{dr } r^{2} \frac{d}{d r} - \frac{l(l+1)}{r^{2}}
    \right] R_{nl}(r) = - \lambda_{nl}^{2} R_{nl}(r)

where we note that :math:`\lambda` only depends
on the indices :math:`n` and :math:`l`.
This adopts the form of the spherical Bessel differential equation,
leading to solutions of the form
:math:`R_{nl} = A j_{l}(\lambda_{nl} r) + B y_{l}(\lambda_{nl} r)`,
where :math:`j_{l}` and :math:`y_{l}` are the spherical Bessel functions
of the first and second kind, respectively.
Since :math:`r = 0` is part of the domain of interest, we exclude :math:`y_{l}`
from our solution, since these functions diverge as :math:`r \rightarrow 0`.
The boundary condition is satisfied by forcing :math:`R(r = 1) = 0`, leading
to the condition :math:`j_{l} ( \lambda_{nlm}) = 0`.

We then write the full solution for the Green function in the form

.. math::
    G(\vec{R}|\vec{R}_{0};N) =
    \sum_{nlm} \psi_{nlm} (\vec{R})
    \psi_{nlm}^{*} (\vec{R}_{0})
    \exp \left[
    - \frac{1}{6} \lambda_{nl}^{2}  \frac{b^{2} N}{a^{2}}
    \right]

where

.. math::
    \psi_{nlm} (\vec{R}) = A_{nlm} j_{l}(\lambda_{nl} R/a) Y_{l}^{m} (\hat{e})

The normalization constant :math:`A_{nlm}` dictates that

.. math::
    \int d \hat{e} \int_{0}^{a} dR R^{2} \psi_{nlm} (\vec{R}) \psi_{n'l'm'}^{*} (\vec{R})
    =
    \delta_{nn'}
    \delta_{ll'}
    \delta_{mm'}


Chain segmental density and average position for free ends
**********************************************************

From this Green function,
we now determine the segmental density and average position of a polymer chain within a
confinement.
Specifically, we address the behavior for both a chain with free ends and a chain
with both ends attached to the sphere surface (as
applied to organization within the yeast nucleus during prophase I of meiosis).
We define the segmental density for a
segment at position :math:`\Delta` for a
chain of total length :math:`N`.
For free ends, the segment density is given by

.. math::
    \rho_{\mathrm{free}}(\vec{R}_{\Delta})
    & = &
    \frac{1}{\mathcal{N}_{\mathrm{free}}} \int d \vec{R} d \vec{R}_{0}
    G(\vec{R}|\vec{R}_{\Delta}; N - \Delta)
    G(\vec{R}_{\Delta}|\vec{R}_{0}; \Delta) \\
    & = &
    \frac{4 a^{-2}}{\mathcal{N}_{\mathrm{free}}}
    \int_{0}^{a} dR R^{2}
    \int_{0}^{a} dR_{0} R_{0}^{2}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    \frac{\sin (n \pi R / a)}{R}
    \frac{\sin (n \pi R_{\Delta} / a)}{R_{\Delta}}
    \frac{\sin (n_{0} \pi R_{\Delta} / a)}{R_{\Delta}}
    \frac{\sin (n_{0} \pi R_{0} / a)}{R_{0}} \\
    &  &
    \hspace{1in}
    \times
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right] \\
    & = &
    \frac{4 a^{2}}{\mathcal{N}_{\mathrm{free}}}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    \frac{(-1)^{n}}{(n \pi)}
    \frac{(-1)^{n_{0}}}{(n_{0} \pi)}
    \frac{\sin (n \pi R_{\Delta} / a) \sin (n_{0} \pi R_{\Delta} / a)}{R_{\Delta}^{2}}
    \\
    &  &
    \hspace{1in}
    \times
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right]

where :math:`\mathcal{N}_{\mathrm{free}}` is a normalization constant that ensures
:math:`\int d \vec{R}_{\Delta} \rho_{\mathrm{free}}(\vec{R}_{\Delta}) = 1`.
With this,
we find the normalization constant to be

.. math::
    \mathcal{N}_{\mathrm{free}}
    = 8 \pi a^{3}
    \sum_{n = 1}^{\infty}
    \frac{1}{(n \pi)^{2}}
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} N}{a^{2}}
    \right]

From this, we find the average squared position of the segment within the sphere
to be

.. math::
    \langle
    R_{\Delta}^{2}
    \rangle_{\mathrm{free}} =
    \frac{16 \pi a^{5}}{\mathcal{N}_{\mathrm{free}}}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    \frac{(-1)^{n}}{(n \pi)}
    \frac{(-1)^{n_{0}}}{(n_{0} \pi)}
    I_{nn_{0}}
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right]

where the integral factor :math:`I_{nn_{0}}` is given by

.. math::
    I_{nn_{0}} =
    \frac{4 n n_{0}(-1)^{n} (-1)^{n_{0}}}{(\pi n^{2} - \pi n_{0}^{2})^{2}}

if :math:`n \ne n_{0}`, and

.. math::
    I_{nn} = \frac{1}{6} - \frac{1}{4 \pi^{2} n^{2}}

for the case :math:`n = n_{0}`.


Chain segmental density and average position for attached ends
**************************************************************

We now find the segmental density for the case where the chain ends are attached to the
sphere surface.
We assume the ends are free to move on the sphere surface, such that
:math:`\vec{R} = (a-\epsilon) \hat{e}` and :math:`\vec{R}_{0} = (a-\epsilon) \hat{e}_{0}`.
Since the Green function approaches zero on the sphere surface, we set the position a
small distance :math:`\epsilon` within the sphere.
We will take the limit :math:`\epsilon \rightarrow 0` after we find the average to
determine the limiting behavior for surface attachment.

With this development, we find the segmental density within the sphere to be

.. math::
    \rho_{\mathrm{surf}}(\vec{R}_{\Delta})
    & = &
    \frac{1}{\mathcal{N}_{\mathrm{surf}}} \int d \hat{e} d \hat{e}_{0}
    G(\vec{R}= (a-\epsilon) \hat{e}|\vec{R}_{\Delta}; N - \Delta)
    G(\vec{R}_{\Delta}|\vec{R}_{0}= (a-\epsilon) \hat{e}_{0}; \Delta) \\
    & = &
    \frac{4 a^{-2}}{\mathcal{N}_{\mathrm{surf}}}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    \frac{\sin [n \pi (a - \epsilon) / a]}{(a - \epsilon)}
    \frac{\sin (n \pi R_{\Delta} / a)}{R_{\Delta}}
    \frac{\sin (n_{0} \pi R_{\Delta} / a)}{R_{\Delta}}
    \frac{\sin [n_{0} \pi (a - \epsilon) / a]}{(a - \epsilon)} \\
    &  &
    \hspace{1in}
    \times
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right] \\
    & = &
    \frac{4\epsilon^{2} a^{-6}}{\mathcal{N}_{\mathrm{surf}}}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    (-1)^{n}(n \pi)
    (-1)^{n_{0}}(n_{0} \pi)
    \frac{\sin (n \pi R_{\Delta} / a) \sin (n_{0} \pi R_{\Delta} / a)}{R_{\Delta}^{2}}
    \\
    &  &
    \hspace{1in}
    \times
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right]

We find the normalization constant to be

.. math::
    \mathcal{N}_{\mathrm{surf}}
    = 8 \pi \epsilon^{2} a^{-5}
    \sum_{n = 1}^{\infty}
    (n \pi)^{2}
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} N}{a^{2}}
    \right]

From this, we find the average squared position of the segment within the sphere
to be

.. math::
    \langle
    R_{\Delta}^{2}
    \rangle_{\mathrm{surf}} =
    \frac{16 \pi \epsilon^{2} a^{-3}}{\mathcal{N}_{\mathrm{surf}}}
    \sum_{n = 1}^{\infty}
    \sum_{n_{0} = 1}^{\infty}
    (-1)^{n}(n \pi)
    (-1)^{n_{0}}(n_{0} \pi)
    I_{nn_{0}}
    \exp \left[
    - \frac{1}{6} (n \pi)^{2}  \frac{b^{2} (N-\Delta)}{a^{2}}
    \right]
    \exp \left[
    - \frac{1}{6} (n_{0} \pi)^{2}  \frac{b^{2} \Delta}{a^{2}}
    \right]

where the integral factor :math:`I_{nn_{0}}` is given by

.. math::
    I_{nn_{0}} =
    \frac{4 n n_{0}(-1)^{n} (-1)^{n_{0}}}{(\pi n^{2} - \pi n_{0}^{2})^{2}}

if :math:`n \ne n_{0}`, and

.. math::
    I_{nn} = \frac{1}{6} - \frac{1}{4 \pi^{2} n^{2}}

for the case :math:`n = n_{0}`.

Functions contained with the 'poly_confine' module
--------------------------------------------------

.. automodule:: wlcstat.poly_confine
    :members:

