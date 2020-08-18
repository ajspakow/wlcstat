import matplotlib.pyplot as plt
import numpy as np
length_kuhn = np.logspace(-4, 4, 50)
plt.loglog(length_kuhn, length_kuhn)
plt.show()