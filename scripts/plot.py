import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# PLUTO units
UNIT_DENSITY = 1e-10
UNIT_VELOCITY = 1e5
UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY**2

# start [stop, step]
if len(sys.argv) == 2:
    start = str(sys.argv[1])
    files = ['data.'+start.zfill(4)+'.tab']
else:
    start, stop, step = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    file_nums = list(map(str, np.arange(start, stop+1, step)))
    files = []
    for i in range(len(file_nums)):
        fname = 'data.'+file_nums[i].zfill(4)+'.tab'
        if os.path.exists(fname): # skip files that don't exist
            files.append(fname)

# oldest is purple, proceeds through rainbow to newest (red)
colors = plt.cm.rainbow(np.linspace(0,1, len(files)))

for i, filename in enumerate(files):
    
    r, rho, v, P = np.loadtxt(filename, usecols=(0,2,3,6), unpack=True)
    plt.figure(0)
    plt.loglog(r, rho*UNIT_DENSITY, color=colors[i])
    plt.ylabel("Density (g/cm3)")
    plt.title("rho")

    plt.figure(1)
    plt.loglog(r, v, color=colors[i])
    plt.ylabel("Velocity (km/s)")
    plt.title("v")

    plt.figure(2)
    plt.loglog(r, P*UNIT_PRESSURE, color=colors[i])
    plt.ylabel("Pressure (g/(cm s^2))")
    plt.title("P")

    plt.figure(3)
    plt.loglog(r, P*2.3*1.67e-24/(1.38e-16*rho*1e-10), color=colors[i])
    plt.title("Estimated Temperature (K)")

plt.show()
