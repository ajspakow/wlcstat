"""
.. module:: wlcave
   :platform: Mac, Unix
   :synopsis: Python module to calculate average quantities for the wormlike chain model

.. moduleauthor:: Andrew Spakowitz <ajspakow@stanford.edu>


"""

import numpy as np


def r2_ave(length_kuhn, dimensions=3):
    """
    r2_ave - Calculate the average end-to-end distance squared for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    r2 : float (array)
        The mean-square end-to-end distance for the wormlike chain model

    Notes: See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008).
    """
    r2 = 2 * (length_kuhn / (dimensions - 1)
              - (1 - np.exp(-(dimensions - 1) * length_kuhn)) / (dimensions - 1) ** 2)

    return r2


def rz4_ave(length_kuhn, dimensions=3):
    """
    rz4_ave - Calculate the < rz ** 4 > for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float (array)
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    rz4 : float (array)
        The mean-square end-to-end distance for the wormlike chain model

    Notes: See Mehraeen, et al, Phys. Rev. E, 77, 061803 (2008).
    """
    rz4 = 2 * (length_kuhn / (dimensions - 1)
              - (1 - np.exp(-(dimensions - 1) * length_kuhn)) / (dimensions - 1) ** 2)

    return rz4


def rg2_ave(length_kuhn):
    pass


def r_inv_ave():
    pass

