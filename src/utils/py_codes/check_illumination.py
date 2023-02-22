import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

filename = sys.argv[1]

nx = int(catch_parameter(filename, "x_samples"))
ny = int(catch_parameter(filename, "y_samples"))
nz = int(catch_parameter(filename, "z_samples"))

dx = float(catch_parameter(filename, "x_spacing"))
dy = float(catch_parameter(filename, "y_spacing"))
dz = float(catch_parameter(filename, "z_spacing"))

vp_location = catch_parameter(filename, "illumination_folder")[1:-1]

vp = readBinaryVolume(nz, nx, ny, vp_location + f"illumination_nz{nz}_nx{nx}_ny{ny}.bin")

geometry_folder = catch_parameter(filename, "geometry_folder")[1:-1]

shots = np.loadtxt(geometry_folder + "xyz_shots.txt", delimiter = ',')
nodes = np.loadtxt(geometry_folder + "xyz_nodes.txt", delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([nz/2, nx/2, ny/2], dtype = int) 
dh = np.array([dx, dy, dz])

unit_shape = (3,)

unit_shots = np.zeros((1,3))
unit_nodes = np.zeros((1,3))

if np.shape(shots) == unit_shape and np.shape(nodes) != unit_shape:

    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    check_illumination(vp, unit_shots, nodes, dh, slices, subplots)

elif np.shape(shots) != unit_shape and np.shape(nodes) == unit_shape:

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]
    
    check_illumination(vp, shots, unit_nodes, dh, slices, subplots)

elif np.shape(shots) == unit_shape and np.shape(nodes) == unit_shape:

    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]
    
    check_illumination(vp, unit_shots, unit_nodes, dh, slices, subplots)

else:

    check_illumination(vp, shots, nodes, dh, slices, subplots)

plt.show()



