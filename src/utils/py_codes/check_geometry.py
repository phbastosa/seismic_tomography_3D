import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

filename = sys.argv[1]

x_samples = int(catch_parameter(filename, "x_samples"))
y_samples = int(catch_parameter(filename, "y_samples"))
z_samples = int(catch_parameter(filename, "z_samples"))

x_spacing = float(catch_parameter(filename, "x_spacing"))
y_spacing = float(catch_parameter(filename, "y_spacing"))
z_spacing = float(catch_parameter(filename, "z_spacing"))

vp_location = catch_parameter(filename, "vp_file")[1:-1]

vp = readBinaryVolume(z_samples, x_samples, y_samples, vp_location)

geometry_folder = catch_parameter(filename, "geometry_folder")[1:-1]

shots = np.loadtxt(geometry_folder + "xyz_shots.txt", delimiter = ',')
nodes = np.loadtxt(geometry_folder + "xyz_nodes.txt", delimiter = ',')

subplots = np.array([1, 1], dtype = int)
slices = np.array([z_samples/2, x_samples/2, y_samples/2], dtype = int) # [xy, zy, zx]
dh = np.array([x_spacing, y_spacing, z_spacing])

unit_shape = (3,)

unit_shots = np.zeros((1,3))
unit_nodes = np.zeros((1,3))

if np.shape(shots) == unit_shape and np.shape(nodes) != unit_shape:

    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    check_geometry(vp, unit_shots, nodes, dh, slices, subplots)

elif np.shape(shots) != unit_shape and np.shape(nodes) == unit_shape:

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]
    
    check_geometry(vp, shots, unit_nodes, dh, slices, subplots)

elif np.shape(shots) == unit_shape and np.shape(nodes) == unit_shape:

    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]
    
    check_geometry(vp, unit_shots, unit_nodes, dh, slices, subplots)

else:

    check_geometry(vp, shots, nodes, dh, slices, subplots)

plt.show()



