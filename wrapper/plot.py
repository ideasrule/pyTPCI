import numpy as np
import matplotlib.pyplot as plt
import sys

for i, filename in enumerate(sys.argv[1:]):
    #if i % 2 != 1:
    #    continue
    
    r, rho, v, P = np.loadtxt(filename, usecols=(0,2,3,6), unpack=True)
    #r, rho, v = np.loadtxt(filename, usecols=(0,2,3), unpack=True)
    plt.figure(0)
    plt.loglog(r, rho)
    plt.title("rho")

    plt.figure(1)
    plt.semilogx(r, v)
    plt.title("v")

    plt.figure(2)
    plt.loglog(r, P)
    plt.title("P")

    plt.figure(3)
    plt.loglog(r, P*2.3*1.67e-24/(1.38e-16*rho*1e-10))
    plt.title("T?")

plt.show()
