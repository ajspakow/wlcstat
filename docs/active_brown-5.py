import matplotlib.pyplot as plt
import numpy as np
from wlcstat.active_brown import *

num_t = 100
t_val_0 = 1e-4
t_val_f = 1e6
t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
N=100
Ndel=25
fa = 10
plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
ka_vector = np.logspace(-2,2,5)
for i in range(len(ka_vector)):
    ka = ka_vector[i]
    mscd = mscd_active(t_val, 1, ka, fa, Ndel, N, 1, 20000)
    plt.loglog(t_val, mscd, label = '$K_{A}$ = ' + str(ka))
mscd_inf = 2 * 2 * Ndel
plt.legend()
plt.loglog(t_val, mscd_inf + 0*t_val,'--')
plt.xlabel(r'$time, t$')
plt.ylabel(r'$MSCD$')
plt.tight_layout()
plt.show()