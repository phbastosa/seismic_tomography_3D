import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

parameters = sys.argv[1]

geometry_folder = catch_parameter(parameters, "geometry_folder")[1:-1]

shots = np.loadtxt(geometry_folder + "xyz_shots.txt", delimiter = ',')
nodes = np.loadtxt(geometry_folder + "xyz_nodes.txt", delimiter = ',')

first_arrival = readBinaryArray(len(nodes), f"../outputs/first_arrivals/times_nr{len(nodes)}_shot_1.bin")

tcut = float(catch_parameter(parameters,"tcut"))

plt.figure(1, figsize = (15, 5))

plt.plot(first_arrival, "o", markersize = 2)

plt.title("First arrival for shot 1", fontsize = 18)
plt.xlabel("Trace index", fontsize = 15)
plt.ylabel("Time [s]", fontsize = 15)
plt.ylim([0, tcut])
plt.xlim([0,len(nodes)])

plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

