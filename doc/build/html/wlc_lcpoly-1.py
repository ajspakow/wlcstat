import matplotlib.pyplot as plt
import numpy as np
from wlcstat.wlc_lcpoly import *

lam_0 = 0
lam_f = 25
n_lam = 500
lam_vec = np.linspace(lam_0, lam_f, n_lam)
length_kuhn = 3.333
plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
kappa_vec = np.array([20.7582, 21.0606, 23.6844])
for i_kap in range(len(kappa_vec)):
    kappa = kappa_vec[i_kap]
    m_val = np.zeros(n_lam)
    q_val = np.zeros(n_lam)
    f_val = np.zeros(n_lam)
    for i_lam in range(n_lam):
        lam = lam_vec[i_lam]
        q_val[i_lam] = q_lcpoly(length_kuhn,lam)
        m_val[i_lam] = lam / kappa
        f_val[i_lam] = kappa * length_kuhn * m_val[i_lam] ** 2 / 3 - np.log(q_val[i_lam])
    plt.plot(m_val, f_val,'-')
plt.ylabel(r'$\beta \Delta f$')
plt.xlabel(r'$m$')
plt.ylim((-0.6, 0.6))
plt.xlim((0, 0.8))
plt.tight_layout()
plt.show()