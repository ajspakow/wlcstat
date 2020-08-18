import numpy as np


def r2_ave(length_kuhn):
    """
    r2_ave

    Calculate the average end-to-end distance squared for the wormlike chain model

    """
    r2 = length_kuhn - 0.5 * (1 - np.exp(-2 * length_kuhn))

    return r2


def rg_ave(length_kuhn):
    pass


def r_inv_ave():
    pass

