import numpy as np
import matplotlib.pyplot as plt

from functions import *

from scipy.ndimage import gaussian_filter

nx = 351
ny = 251
nz = 51

dx = 20
dy = 20
dz = 20

types = ["high_pod", "high_fim", "high_fsm"]
iterations = [5, 10]

vp_location = "../../../outputs/recovered_models/"
geometry_folder = "../../../inputs/geometry/"

shots = np.loadtxt(geometry_folder + "xyz_shots.txt", delimiter = ',')
nodes = np.loadtxt(geometry_folder + "xyz_nodes.txt", delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([nz/2, ny/2 , nx/2 + 10], dtype = int) 
dh = np.array([dx, dy, dz])

for it in iterations:
    for type in types:
        vp = readBinaryVolume(nz, nx, ny, vp_location + f"{type}_model_iteration_{it}.bin")
        vp = 1 / gaussian_filter(1 / vp, 2.0)
        check_model(vp, dh, slices, subplots)
        # check_geometry(vp, shots, nodes, dh, slices, subplots)
        plt.savefig(f"figures/{type}_model_iteration_{it}.png", dpi = 200)

labels = ["Podvin & Lecomte (1991)", "Jeong & Whitaker (2008)", "Noble, Gesret and Belayouni (2014)"] 

types = ["pod", "fim", "fsm"]

plt.figure(2, figsize = (15, 5))
for k, type in enumerate(types):
    convergency = np.loadtxt(f"../../../outputs/convergence/{type}_convergency.txt")
    plt.plot(convergency, "o-", label = labels[k])

plt.xticks(np.linspace(0,15,16, dtype=int),np.linspace(0,15,16, dtype=int))
plt.xlim([-0.1,15.1])
plt.title("Convergência", fontsize = 18)
plt.xlabel("Número de iterações", fontsize = 15)
plt.ylabel("Norma 2 da função objetivo", fontsize = 15)
plt.legend(loc = "upper right", fontsize = 15)
plt.tight_layout()
plt.savefig("figures/convergency.png", dpi = 200)
plt.show()