import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

parameters = sys.argv[1]

nx = int(catch_parameter(parameters, "x_samples"))
ny = int(catch_parameter(parameters, "y_samples"))
nz = int(catch_parameter(parameters, "z_samples"))

dx = float(catch_parameter(parameters, "x_spacing"))
dy = float(catch_parameter(parameters, "y_spacing"))
dz = float(catch_parameter(parameters, "z_spacing"))

vp_location = catch_parameter(parameters, "vp_file")[1:-1]
geometry_folder = catch_parameter(parameters, "geometry_folder")[1:-1]
travel_time_file = "../outputs/travel_times/eikonal_nz105_nx401_ny401_shot_1.bin"

vp = readBinaryVolume(nz, nx, ny, vp_location)
tt = readBinaryVolume(nz,nx,ny,travel_time_file)

shots = np.loadtxt(geometry_folder + "xyz_shots.txt", delimiter = ',')
nodes = np.loadtxt(geometry_folder + "xyz_nodes.txt", delimiter = ',')

subplots = np.array([1, 1], dtype = int)
dh = np.array([dx, dy, dz])

unit_shape = (3,)

unit_shots = np.zeros((1,3))
unit_nodes = np.zeros((1,3))

if np.shape(shots) == unit_shape and np.shape(nodes) != unit_shape:
    
    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    slices = np.array([shots[2]/dz, shots[0]/dx, shots[1]/dy], dtype = int) 

    check_travel_time(vp, tt, unit_shots, nodes, dh, slices, subplots)

elif np.shape(shots) != unit_shape and np.shape(nodes) == unit_shape:

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]

    slices = np.array([shots[0,2]/dz, shots[0,0]/dx, shots[0,1]/dy], dtype = int) 

    check_travel_time(vp, tt, shots, unit_nodes, dh, slices, subplots)

elif np.shape(shots) == unit_shape and np.shape(nodes) == unit_shape:

    unit_shots[:,0] = shots[0]
    unit_shots[:,1] = shots[1]
    unit_shots[:,2] = shots[2]

    unit_nodes[:,0] = nodes[0]
    unit_nodes[:,1] = nodes[1]
    unit_nodes[:,2] = nodes[2]


    slices = np.array([shots[2]/dz, shots[0]/dx, shots[1]/dy], dtype = int) 

    check_travel_time(vp, tt, unit_shots, unit_nodes, dh, slices, subplots)

else:

    slices = np.array([shots[0,2]/dz, shots[0,0]/dx, shots[0,1]/dy], dtype = int) 

    check_travel_time(vp, tt, shots, nodes, dh, slices, subplots)

plt.show()
