import matplotlib.pyplot as plt
import numpy as np
from wlcstat.wlcgreen import *

length_kuhn = 0.1
alpha_max = 25
num_k = 100
k_val_0 = 1e-1
k_val_f = 1e5
k_val = np.logspace(np.log10(k_val_0), np.log10(k_val_f), num_k)
mu=0
dimensions = 3
num_poles = min(18, alpha_max + 1 - mu)
poles = np.zeros((num_k, num_poles), dtype=type(1 + 1j))
for i_k_val in range(num_k):
    poles_k_val, resi_k_val = eval_poles_and_residues(k_val[i_k_val],mu,dimensions)
    for i_pole in range(num_poles):
        poles[i_k_val, i_pole] = poles_k_val[i_pole]

plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)

for i_pole in range(num_poles):
    plt.loglog(k_val, np.exp(np.real(poles[:, i_pole]) * length_kuhn))

plt.ylim((10 ** -14, 1))
plt.xlim((10 ** -1, 10 ** 5))
plt.xlabel(r'$K = (2l_{p}) k$')
plt.ylabel(r'$\exp \left[ \mathrm{Real} (\mathcal{E}_{\alpha}) N \right]$')
plt.tight_layout()
plt.show()