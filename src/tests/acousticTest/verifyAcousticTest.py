import numpy as np
import matplotlib.pyplot as plt

def readBinaryArray(n, filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

#------------------------------------------------------------------------------------
dt = 1e-3
tlag = 0.15

nt = int(2.0 * tlag / dt + 1)

t = np.arange(nt)*dt

source = readBinaryArray(nt, f"outputs/source_Nt{nt}.bin")

sourceFFT = np.fft.fft(source)
f = np.fft.fftfreq(nt, dt) 

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15,5))

ax[0].plot(t, source)

ax[1].plot(f, np.abs(sourceFFT))
ax[1].set_xlim([0, 30])

plt.tight_layout()
plt.show()

#------------------------------------------------------------------------------------
nb = 50

damp1D = readBinaryArray(nb, f"outputs/damp1D_{nb}_samples.bin")
damp2D = readBinaryArray(nb**2, f"outputs/damp2D_{int(nb**2)}_samples.bin")
damp3D = readBinaryArray(nb**3, f"outputs/damp3D_{int(nb**3)}_samples.bin")

damp2D = np.reshape(damp2D, [nb,nb], order = 'F')
damp3D = np.reshape(damp3D, [nb,nb,nb], order = 'F')

plt.figure(1)
plt.plot(damp1D)

plt.figure(2)
plt.imshow(damp2D)


plt.figure(3)
plt.subplot(331)
plt.imshow(damp3D[int(nb/2),:,:])

plt.subplot(332)
plt.imshow(damp3D[:,int(nb/2),:])

plt.subplot(333)
plt.imshow(damp3D[:,:,int(nb/2)])

plt.subplot(334)
plt.imshow(damp3D[int(2*nb/8),:,:])

plt.subplot(335)
plt.imshow(damp3D[:,int(2*nb/8),:])

plt.subplot(336)
plt.imshow(damp3D[:,:,int(2*nb/8)])

plt.subplot(337)
plt.imshow(damp3D[int(6*nb/8),:,:])

plt.subplot(338)
plt.imshow(damp3D[:,int(6*nb/8),:])

plt.subplot(339)
plt.imshow(damp3D[:,:,int(6*nb/8)])

plt.show()


sx, sy, sz = np.loadtxt("outputs/shots.txt", delimiter = ",", unpack = True)
rx, ry, rz = np.loadtxt("outputs/nodes.txt", delimiter = ",", unpack = True)

plt.figure(4)
plt.scatter(rx,ry)
plt.scatter(sx,sy)
plt.show()