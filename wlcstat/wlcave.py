import numpy as np


def r2_ave(length_kuhn, dimensions=3):
    r"""
    r2_ave - Calculate the average end-to-end distance squared :math:`\langle R^{2} \rangle / (2 l_{p})^{2}`
    for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    r2 : float (array)
        The mean-square end-to-end distance for the wormlike chain model (non-dimensionalized by :math:`2 l_{p})`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """
    r2 = 2 * (length_kuhn / (dimensions - 1)
              - (1 - np.exp(-(dimensions - 1) * length_kuhn)) / (dimensions - 1) ** 2)

    return r2


def rg2_ave(length_kuhn, dimensions=3):
    r"""
    rg2_ave - Calculate the radius of gyration
    :math:`\langle \vec{R}_{G}^{2} \rangle / (2 l_{p})^{2}` for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    rg2 : float (array)
        The mean-square radius of gyration for the wormlike chain model (non-dimensionalized by :math:`2 l_{p})`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """
    rg2 = 2 * (length_kuhn / (6 * (dimensions - 1)) - 1 / (2 * (dimensions - 1) ** 2)
               + length_kuhn ** -1 / (dimensions - 1) ** 3
               - length_kuhn ** -2 * (1 - np.exp(-(dimensions - 1) * length_kuhn)) / (dimensions - 1) ** 4)

    return rg2


def rz4_ave(length_kuhn, dimensions=3):
    r"""
    rz4_ave - Calculate the 4th moment of the end-to-end distribution
    :math:`\langle R_{z}^{4} \rangle / (2 l_{p})^{4}` for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    rz4 : float (array)
        The mean-square end-to-end distance for the wormlike chain model (non-dimensionalized by :math:`2 l_{p})`

    Notes
    -----
    See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008). (Ref [Mehraeen2008]_)
    """

    # Calculate the coefficients from the ladder operations on the hyperspherical harmonics (dimensions)
    a1 = np.sqrt(1 / dimensions)
    a2 = np.sqrt(2 * (dimensions - 1) / ((dimensions + 2) * dimensions))

    # Calculate the contributions from the 2 contributing stone-fence diagrams
    diagram1 = (length_kuhn ** 2 / (2 * (dimensions - 1) ** 2) - 2 * length_kuhn / (dimensions - 1) ** 3
                + 3 / (dimensions - 1) ** 4
                - length_kuhn * np.exp(-(dimensions - 1) * length_kuhn) / (dimensions - 1) ** 3
                - 3 * np.exp(-(dimensions - 1) * length_kuhn) / (dimensions - 1) ** 4)
    diagram1 *= a1 ** 4

    diagram2 = (length_kuhn / (2 * dimensions * (dimensions - 1) ** 2)
                - (5 * dimensions - 1) / (4 * dimensions ** 2 * (dimensions - 1) ** 3)
                + np.exp(-(dimensions - 1) * length_kuhn) * (
                        length_kuhn / ((dimensions + 1) * (dimensions - 1) ** 2)
                        + (dimensions + 3) / ((dimensions - 1) ** 3 * (dimensions + 1) ** 2))
                + np.exp(-2 * dimensions * length_kuhn) / (4 * dimensions ** 2 * (dimensions + 1) ** 2))
    diagram2 *= a1 ** 2 * a2 ** 2

    rz4 = 24 * (diagram1 + diagram2)

    return rz4


def gen_conf_wlc(length_kuhn, num_beads, b=1, tangent_calc=False):
    delta = length_kuhn / (num_beads - 1)
    return
