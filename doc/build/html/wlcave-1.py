import matplotlib.pyplot as plt
import numpy as np
from wlcstat.wlcave import *

length_kuhn = np.logspace(-4, 4, 50)
r2 = r2_ave(length_kuhn)
plt.figure(figsize=(5,4))

font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
plt.loglog(length_kuhn, r2)
plt.xlabel(r'Chain length $N = L/(2l_{p})$')
plt.ylabel(r'End-to-end distance $\langle R^{2} \rangle$')
plt.tight_layout()
plt.show()