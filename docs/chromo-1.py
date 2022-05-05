import matplotlib.pyplot as plt
import numpy as np
from wlcstat.chromo import *

links = 36 * np.ones(50)
r2_t0chk, lengthDNA, kuhn, w_ins, w_outs = r2_kinked_twlc(links, lp=100000, lt=100000)
r2, lengthDNA, kuhn, w_ins, w_outs = r2_kinked_twlc(links)
rots, pos = minimum_energy_no_sterics_linker_only(links)
r2_t0 = np.array([np.linalg.norm(pos[i, :] - pos[0, :])**2 for i in range(pos.shape[0])])
plt.figure(figsize=(10,8))
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
length_total = np.cumsum(links + w_outs[0:-1] + w_ins[1:] + 1)
plt.plot(length_total, r2,'o-','color','C0')
plt.plot(length_total, r2_t0chk,'o-','color','C2')
plt.plot(length_total, r2_t0[1:],'.-','color','C3')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Length DNA (basepairs)')
plt.ylabel(r'Mean-square inter-nucleosome distance $\langle R^{2} \rangle$')
plt.tight_layout()
plt.show()