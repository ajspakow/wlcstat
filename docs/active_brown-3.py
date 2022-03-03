import matplotlib.pyplot as plt
import numpy as np
from wlcstat.active_brown import *

num_t = 100
t_val_0 = 1e-6
t_val_f = 1e6
t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
N=100
fa = 10
plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
ka_vector = np.logspace(-2,2,5)
for i in range(len(ka_vector)):
    ka = ka_vector[i]
    msd = msd_active(t_val, 1, ka, fa, N, 1, 40000)
    plt.loglog(t_val, msd, label = '$K_{A}$ = ' + str(ka))
msd_inf = 6 / N * t_val * (1 + 0.5 * fa ** 2)
plt.loglog(t_val, msd_inf,':')
plt.legend()
plt.xlabel(r'$time, t$')
plt.ylabel(r'$MSD$')
plt.tight_layout()
plt.show()