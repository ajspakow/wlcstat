import matplotlib.pyplot as plt
import numpy as np
from wlcstat.wlcgreen import *

num_k = 100
k_val_0 = 1e-1
k_val_f = 1e3
k_val = np.logspace(np.log10(k_val_0), np.log10(k_val_f), num_k)
mu=0
dimensions = 3
num_poles = min(12, 26-mu)
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
    plt.semilogx(k_val, np.real(poles[:, i_pole]))

plt.xlabel(r'$K = (2l_{p}) k$')
plt.ylabel(r'Real ($\mathcal{E}_{\alpha}$)')
plt.tight_layout()
plt.show()