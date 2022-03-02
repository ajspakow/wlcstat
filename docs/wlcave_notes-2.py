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