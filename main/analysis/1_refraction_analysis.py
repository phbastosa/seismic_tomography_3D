import numpy as np
import matplotlib.pyplot as plt

from all_functions import boxPlot, readBinaryArray, readBinaryVolume, analiticalRefraction

nx = 481
ny = 41
nz = 51

dh = 12.5

v = np.array([2000,3500])
z = np.array([500]) 

refractiveModel = np.zeros((nz,nx,ny))

refractiveModel[:int(z[0]/dh),:,:] = v[0]
refractiveModel[int(z[0]/dh):,:,:] = v[1]

refractiveModel.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/refractiveModel_{nz}x{nx}x{ny}_{dh:.1f}m.bin")

# shots = np.array([500, 250, 0])

# nr = 26

# nodes = np.zeros((3, nr))

# nodes[0,:] = np.linspace(1000, 3500, nr)
# nodes[1,:] = np.ones(nr) * 250
# nodes[2,:] = np.zeros(nr)

# dh = np.array([dh,dh,dh])
# slices = np.array([int(nz/2), int(nx/2), int(ny/2)], dtype = int)

# n = 20319

# xRay = readBinaryArray(n, f"../../outputs/rayPosition/xRay_{n}_shot_1.bin")
# yRay = readBinaryArray(n, f"../../outputs/rayPosition/yRay_{n}_shot_1.bin")
# zRay = readBinaryArray(n, f"../../outputs/rayPosition/zRay_{n}_shot_1.bin")
# iRay = readBinaryArray(n, f"../../outputs/rayPosition/iRay_{n}_shot_1.bin")

# xzRay = np.zeros((2, n))

# xzRay[0,:] = xRay
# xzRay[1,:] = zRay

# eikonal = readBinaryVolume(nz,nx,ny,f"../../outputs/travelTimes/eikonal_nz{nz}_nx{nx}_ny{ny}_shot_1.bin")

# boxPlot(refractiveModel, dh, shots, nodes, slices, xzRay, eikonal, colorbarFix = 2.0)
# plt.show()

# #--------------------------------

# offsets = nodes[0,:] - shots[0]

# td, tr = analiticalRefraction(v,z,offsets)

# plt.plot(td)
# plt.plot(tr[0])

# plt.gca().invert_yaxis()
# plt.show()
