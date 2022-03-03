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