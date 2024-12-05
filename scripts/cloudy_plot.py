import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# replace with planet radius (Earth radii)
unit_length = 12.54 * 6.378136e8


# start, [stop, step]
if len(sys.argv) == 2:
    start = str(sys.argv[1])
    files = ['cl_data.'+start.zfill(4)+'.over.tab']
else:
    start, stop, step = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    file_nums = list(map(str, np.arange(start, stop+1, step)))
    files = []
    for i in range(len(file_nums)):
        fname = 'cl_data.'+file_nums[i].zfill(4)+'.over.tab'
        if os.path.exists(fname): # skip files that don't exist
            files.append(fname)

# oldest is purple, proceeds through rainbow to newest (red)
colors = plt.cm.rainbow(np.linspace(0, 1, len(files)))

for i, file in enumerate(files):
    depth, rho, n, mu, Te = np.loadtxt(file, usecols=range(5), unpack=True)
    cloudy_radii = 1 + (np.max(depth) - depth) / unit_length
    plt.plot(cloudy_radii, Te, color=colors[i])

plt.show()
