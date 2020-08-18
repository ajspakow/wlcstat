import matplotlib.pyplot as plt
import numpy as np
length_kuhn = np.logspace(-4, 4, 50)
r2 = r2_ave(length_kuhn)
plt.loglog(length_kuhn, r2)
plt.show()