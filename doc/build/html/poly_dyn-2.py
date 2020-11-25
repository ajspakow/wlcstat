import matplotlib.pyplot as plt
import numpy as np
from wlcstat.poly_dyn import *

mu = 4
cells = [generate_example_cell(mu) for i in range(5)]
ax, all_closest_links = draw_cells(cells)
plt.show()
t = np.logspace(1,4,50).astype(float)
plt.figure()
for i, linkages in enumerate(cells):
    plt.loglog(t, model_mscd(t,linkages), label = 'Cell ' + str(i+1))
plt.legend()
font = {'family' : 'serif',
    'weight':'normal',
    'size': 18}
plt.rc('font', **font)
plt.xlabel(r'time (sec)')
plt.ylabel(r'$MSCD$ ($\mu m^{2}$)')
plt.show()