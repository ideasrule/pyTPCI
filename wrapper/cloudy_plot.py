import numpy as np
import matplotlib.pyplot as plt
import sys

unit_length = 2.9 * 6.378136e8

for filename in sys.argv[1:]:
    depth, rho, n, mu, Te = np.loadtxt(filename, usecols=range(5), unpack=True)
    cloudy_radii = 1 + (np.max(depth) - depth) / unit_length
    plt.plot(cloudy_radii, Te)

plt.show()
