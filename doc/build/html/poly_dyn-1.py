import matplotlib.pyplot as plt
import numpy as np
from wlcstat.poly_dyn import *

num_t = 100
t_val_0 = 1e-4
t_val_f = 1e6
t_val = np.logspace(np.log10(t_val_0), np.log10(t_val_f), num_t)
N=100
Ndel=25
mscdl = linear_mscd(t_val, 1, Ndel, N, 1, 20000)
mscdr = ring_mscd(t_val, 1, Ndel, N, 1, 20000)
msd = linear_mid_msd(t_val, 1, N, 1, 10000)
mscdl_inf = 2 * 2 * Ndel
mscdr_inf = 2 * 1 / (1 / (2 * Ndel) + 1  / (N - 2 * Ndel))
plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
plt.loglog(t_val, mscdl)
plt.loglog(t_val, mscdr)
plt.loglog(t_val, 2 * msd)
plt.loglog(t_val, 2 * 6 * t_val /N,'--')
plt.loglog(t_val, mscdl_inf + 0*t_val,'--')
plt.loglog(t_val, mscdr_inf + 0*t_val,'--')
plt.xlabel(r'$t$')
plt.ylabel(r'$MSCD$')
plt.tight_layout()
plt.show()