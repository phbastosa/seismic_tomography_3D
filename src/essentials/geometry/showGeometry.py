import numpy as np
import matplotlib.pyplot as plt

nodes = np.loadtxt("nodes_n121.bin", delimiter=",")
shots = np.loadtxt("shots_n567.bin", delimiter=",")

plt.figure(1)
plt.scatter(shots[:,0], shots[:,1])
plt.scatter(nodes[:,0], nodes[:,1])
plt.show()