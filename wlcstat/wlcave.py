import numpy as np


def r2_ave(length_kuhn, dimensions=3):
    """
    r2_ave

    Calculate the average end-to-end distance squared for the wormlike chain model

    Parameters
    ----------
    length_kuhn : float
        The length of the chain in Kuhn lengths
    dimensions : int
        The number of dimensions (default to 3 dimensions)

    Returns
    -------
    r2 : float
        The mean-square end-to-end distance for the wormlike chain model
    """
    r2 = length_kuhn - 0.5 * (1 - np.exp(-2 * length_kuhn))

    return r2


def rg2_ave(length_kuhn):
    pass


def r_inv_ave():
    pass

