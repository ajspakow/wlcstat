.. _chromo:
.. .. automodule:: chromo

Chromosomal DNA
===============

We present a model ([Beltran2019]_) combining thermal fluctuations in the DNA linkers with the geometric effects of experimentally-relevant linker length heterogeneity. We show that almost any linker length heterogeneity is sufficient to produce the disordered chromatin structures thatare now believed to dominate nuclear architecture. The intuition behind our structural claims extends to any polymer composed of aperiodic kinks, such as the dihedral "kinks" found at junctions of block copolymers.

.. figure:: figures/nucleosome.png
    :width: 600
    :align: center
    :alt: Nucleosomal DNA structure with variable length linkers

    The structure of a nucleosome with straight linkers extrapolated from the entry (:math:`\Omega_\text{entry}`) and exit (:math:`\Omega_\text{exit}`) orientations of the bound DNA.

We model each DNA linker as a twistable wormlike chain (TWLC), and the nucleosomes as the points where these linker strands connect. We label the orientation of the DNA entering the *i* th nucleosome by the matrix :math:`\Omega^{(i)}_\text{entry}`. The exit orientation :math:`\Omega^{(i)}_\text{exit}` must then satisfy :math:`\Omega^{(i)}_\text{exit} = \Omega^{(i)}_\text{entry} \cdot \Omega_\text{kink}`, where :math:`\Omega_\text{kink}` is the fixed "kink" rotation.

We represent a TWLC of length :math:`L` as a space curve :math:`\vec{R}(s)`, where :math:`0 < s < L`. The chain's orientation at each point along this curve, :math:`\Omega(s)`, is represented by an orthonormal triad :math:`\vec{t}_{i}`, where :math:`\vec{t}_{3} = \partial_s \vec{R}(s)`. We track the bend and twist of our polymer via the angular "velocity" vector :math:`\vec{\omega}(s)`, which operates as :math:`\partial_s \vec{t}_{i}(s) = \vec{\omega}(s) \times \vec{t}_{i}(s)`. The Green's function

.. math::
    G_{\text{TWLC}} (\vec{R}, \Omega | \Omega_0;L_1 ) =\!
    \int_{\Omega(0)=\Omega_0}^{\Omega(s)=\Omega}
    \hspace{-0.3in}
    \mathcal{D} [\Omega(s)]
    e^{-\beta \mathcal{E}}\delta(\vec{R} - \!\int_{0}^{L_1} \!\!\vec{t_3} ds),

of the first linker represents the probability that a   polymer of length :math:`L_1` that begins at the origin with initial orientation :math:`\Omega_0` ends at position :math:`\vec{R}` with end orientation :math:`\Omega`. For a TWLC with no kinks, the energy is quadratic in bending and twist deformations

.. math::
    \beta \mathcal{E} = \frac{l_p}{2}\int_{0}^{L_1} \!\ \! ds \,
    (\omega_1^2 + \omega_2 ^2) + \frac{l_t}{2}\int_{0}^{L_1} \! \! ds \,
    {\left(\omega_3 - \tau\right)}^2

The natural twist of DNA gives :math:`\tau = 2 \pi (\text{10.5} \, \text{bp})^{-1}`, and we set the persistence length :math:`l_p = 50 \, \text{nm}` and twist persistence length :math:`l_t = 10 \, \text{nm}` to match measurements of DNA elasticity.

Reference [Spakowitz2006]_ solves :math:`G_{\text{TWLC}}` analytically in Fourier space (:math:`\vec{R} \rightarrow \vec{k}`) by computing the coefficients of the Wigner D-function expansion

.. math::
    \hat{G} (\vec{k}, \Omega | \Omega_0;L_1 ) =
    \sum_{\substack{ l;\hspace{0.1em}l_0 m_0 j_0}} \!\!
    g_{l_0 m_{0} j_{0}}^{l m_0 j_0}
    \mathcal{D}_{l}^{m_0j_0}(\Omega) \mathcal{D}_{l_0}^{m_0 j_0 *}(\Omega_0).

To account for :math:`\Omega_\text{kink}`, we rotate the final orientation of the linker DNA, :math:`\Omega = \Omega_\text{entry}`, to :math:`\Omega_\text{exit}` using the formula

.. math::
    \mathcal{D}_l^{m_0j_0}(\Omega \cdot \Omega_\text{kink}) =
    \sum_j \sqrt{\frac{8\pi}{2l+1}} \mathcal{D}_l^{m_0j}(\Omega_\text{kink}) \mathcal{D}_l^{jj_0}(\Omega).

The resulting Green's function combines the effects of a DNA linker and a nucleosome, but is still a Wigner D-function expansion with modified coefficients :math:`B_{l_{0} m_{0} j_{0}}^{l m_0 j}`.

Functions contained with the 'chromo' module
--------------------------------------------

.. automodule:: wlcstat.chromo
    :members:

Example usage of 'r2_kinked_twlc'
---------------------------------

We plot the mean-square inter-nucleosome distance :math:`\langle R^{2} \rangle` versus total DNA length
(i.e. including both linker DNA and nucleosomal DNA) for a fixed linker length of 36 basepairs.
The blue points show predictions using standard physical properties of DNA
[:math:`l_{p} = 50 \, \text{nm}`, :math:`l_{t} = 100 \, \text{nm}`, :math:`\tau = 2 \pi / (10.5 \, \text{basepairs})`].
The orange and green curves show zero temperature predictions (i.e. no conformation fluctuations) using the
'r2_kinked_twlc' code in the limit of :math:`l_{p}, l_{t} \gg 1` (orange) and using a straight-linker algorithm (green),
demonstrating quantitative agreement between these approaches.

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from wlcstat.chromo import *

    links = 36 * np.ones(50)
    r2_t0chk, lengthDNA, kuhn, w_ins, w_outs = r2_kinked_twlc(links, lp=100000, lt=100000)
    r2, lengthDNA, kuhn, w_ins, w_outs = r2_kinked_twlc(links)
    rots, pos = minimum_energy_no_sterics_linker_only(links)
    r2_t0 = np.array([np.linalg.norm(pos[i, :] - pos[0, :])**2 for i in range(pos.shape[0])])
    plt.figure(figsize=(10,8))
    font = {'family' : 'serif',
        'weight':'normal',
        'size': 18}
    plt.rc('font', **font)
    length_total = np.cumsum(links + w_outs[0:-1] + w_ins[1:] + 1)
    plt.plot(length_total, r2,'o-','color','C0')
    plt.plot(length_total, r2_t0chk,'o-','color','C2')
    plt.plot(length_total, r2_t0[1:],'.-','color','C3')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Length DNA (basepairs)')
    plt.ylabel(r'Mean-square inter-nucleosome distance $\langle R^{2} \rangle$')
    plt.tight_layout()
    plt.show()

Example usage of 'gen_chromo_conf' and 'gen_chromo_pymol_file'
--------------------------------------------------------------

We demonstrate the usage of 'gen_chromo_conf' to generate the conformation of a chromosome array with 100 nucleosomes and random linker lengths between 35 and 50 basepairs. The output conformation is used to generate a PyMOL file using 'gen_chromo_pymol_file'. The resulting file 'r_chromo.pdb' is loaded in PyMOL to generate an image of the chromosomal DNA.

.. code:: python

    links = np.random.randint(35, 50, 100)
    r, rdna1, rdna2, rn, un = gen_chromo_conf(links)
    gen_chromo_pymol_file(r, rdna1, rdna2, rn, un, filename='r_chromo.pdb', ring=False)

.. figure:: figures/chromo_conf.png
    :width: 600
    :align: center
    :alt: Chromosome conformation

    Image of a 100-nucleosome array generated in PyMOL using output from 'gen_chromo_conf' and 'gen_chromo_pymol_file'.
