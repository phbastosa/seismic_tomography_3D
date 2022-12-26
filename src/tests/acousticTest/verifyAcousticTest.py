import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryArray(n, filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

def readBinaryMatrix(dim1, dim2, filename):
    data = np.fromfile(filename, dtype = np.float32, count = dim1*dim2)
    return np.reshape(data, [dim1, dim2], order = 'F')

def readBinaryVolume(dim1, dim2, dim3, filename):
    data = np.fromfile(filename, dtype = np.float32, count = dim1*dim2*dim3)
    return np.reshape(data, [dim1, dim2, dim3], order = 'F')

#------------------------------------------------------------------------------------
dt = 1e-3
tlag = 0.15

nt = int(2.0 * tlag / dt + 1)

t = np.arange(nt)*dt

source = readBinaryArray(nt, f"outputs/source_Nt{nt}.bin")

sourceFFT = np.fft.fft(np.insert(source, nt-1, np.zeros(10 * nt)))
f = np.fft.fftfreq(11 * nt, dt) 
mask = f >= 0

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15,5))

ax[0].plot(t - tlag, 1.0 / np.max(source) * source)
ax[0].set_xlim([-tlag, tlag])
ax[0].set_ylim([-0.6, 1.2])
ax[0].set_xlabel("Time [s]", fontsize = 15)
ax[0].set_ylabel("Normalized amplitude", fontsize = 15 )
ax[0].set_title("Ricker wavelet", fontsize = 20)

ax[1].plot(f[mask], 1.0 / np.max(np.abs(sourceFFT[mask])) * np.abs(sourceFFT[mask]))
ax[1].set_xlim([0, 30])
ax[1].set_ylim([0.0, 1.01])
ax[1].set_xlabel("Frequency [Hz]", fontsize = 15)
ax[1].set_ylabel("Normalized amplitude", fontsize = 15 )
ax[1].set_title("Wavelet spectrum", fontsize = 20)

plt.tight_layout()
plt.show()

#------------------------------------------------------------------------------------
nb = 50

damp1D = readBinaryArray(nb, f"outputs/damp1D_{nb}_samples.bin")
damp2D = readBinaryArray(nb**2, f"outputs/damp2D_{int(nb**2)}_samples.bin")
damp3D = readBinaryArray(nb**3, f"outputs/damp3D_{int(nb**3)}_samples.bin")

damp2D = np.reshape(damp2D, [nb,nb], order = 'F')
damp3D = np.reshape(damp3D, [nb,nb,nb], order = 'F')


plane1 = int(nb/2)
plane2 = int(2*nb/8)
plane3 = int(6*nb/8)

imgs = [damp3D[plane1,:,:], damp3D[:,plane1,:], damp3D[:,:,plane1],
        damp3D[plane2,:,:], damp3D[:,plane2,:], damp3D[:,:,plane2],
        damp3D[plane3,:,:], damp3D[:,plane3,:], damp3D[:,:,plane3]]

fig, ax = plt.subplots(3,3, figsize = (8, 6))

k = 0
for i in range(3):
    for j in range(3):
        ax[i,j].imshow(imgs[k], aspect = "auto")
        k += 1

fig.suptitle("Corner 3D for damping", fontsize = 20)

plt.tight_layout()
plt.show()



sx, sy, sz = np.loadtxt("outputs/shots.txt", delimiter = ",", unpack = True)
rx, ry, rz = np.loadtxt("outputs/nodes.txt", delimiter = ",", unpack = True)

nz, nx, ny = 51, 101, 101
dz = dx = dy = 12.5

model = readBinaryVolume(nz, nx, ny, "outputs/model_51x101x101.bin")

fig, ax = plt.subplots(1,3, figsize = (10, 5))

ax[0].imshow(model[0,:,:], vmin = np.min(model), vmax = np.max(model), cmap = "Greys")
ax[0].scatter(rx/dx, ry/dy, s = 2.0, c = "red")
ax[0].scatter(sx/dx, sy/dy, s = 5.0, c = "blue")
ax[0].set_xlabel("Distance x [m]", fontsize = 15)
ax[0].set_ylabel("Distance y [m]", fontsize = 15)
ax[0].set_title("XY plane", fontsize = 20)
ax[0].set_yticks(np.linspace(0, ny-1, 5))
ax[0].set_xticks(np.linspace(0, nx-1, 5))
ax[0].set_yticklabels(np.linspace(0, (ny-1)*dy, 5, dtype = int))
ax[0].set_xticklabels(np.linspace(0, (nx-1)*dx, 5, dtype = int))

ax[1].imshow(model[:,:,0], vmin = np.min(model), vmax = np.max(model), cmap = "Greys")
ax[1].scatter(rx/dx, rz/dz, s = 2.0, c = "red")
ax[1].scatter(sx/dx, sz/dz, s = 5.0, c = "blue")
ax[1].set_xlabel("Distance [m]", fontsize = 15)
ax[1].set_ylabel("Depth [m]", fontsize = 15)
ax[1].set_title("ZX plane", fontsize = 20)
ax[1].set_yticks(np.linspace(0, nz-1, 5))
ax[1].set_xticks(np.linspace(0, nx-1, 5))
ax[1].set_yticklabels(np.linspace(0, (nz-1)*dz, 5, dtype = int))
ax[1].set_xticklabels(np.linspace(0, (nx-1)*dx, 5, dtype = int))

ax[2].imshow(model[:,0,:], vmin = np.min(model), vmax = np.max(model), cmap = "Greys")
ax[2].scatter(ry/dy, rz/dz, s = 2.0, c = "red")
ax[2].scatter(sy/dy, sz/dz, s = 5.0, c = "blue")
ax[2].set_xlabel("Distance [m]", fontsize = 15)
ax[2].set_ylabel("Depth [m]", fontsize = 15)
ax[2].set_title("ZY plane", fontsize = 20)
ax[2].set_yticks(np.linspace(0, nz-1, 5))
ax[2].set_xticks(np.linspace(0, ny-1, 5))
ax[2].set_yticklabels(np.linspace(0, (nz-1)*dz, 5, dtype = int))
ax[2].set_xticklabels(np.linspace(0, (ny-1)*dy, 5, dtype = int))

plt.tight_layout()
plt.show()

nt = 1001
nx = 24**2

seismogram = readBinaryMatrix(nt,nx,"outputs/seismogram_shot_1.bin")

scale = 0.05 * np.std(seismogram)

plt.figure(figsize = (15, 5))

plt.imshow(seismogram, cmap = "Greys", aspect = "auto", vmin = -scale, vmax = scale)

plt.yticks(np.linspace(0,nt-1,11, dtype = int), np.linspace(0,nt-1,11, dtype = int))

plt.xlabel("Traces", fontsize = 15)
plt.ylabel("Time [ms]", fontsize = 15)
plt.title("Seismogram for shot 1", fontsize = 20)

plt.tight_layout()
plt.show()
