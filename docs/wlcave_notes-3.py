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