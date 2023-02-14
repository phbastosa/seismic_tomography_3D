import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter

def readBinaryVolume(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order='F')
    return volume

nz = 81
nx = 301
ny = 301

model = readBinaryVolume(nz,nx,ny,"saltDome3D_downscale_z81_x301_y301.bin")

plt.figure(1,figsize=(18,6))

plt.subplot(131)
plt.imshow(model[:,:,int(ny/2)],aspect="auto",cmap="Greys",vmin=1500,vmax=5000)
cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.subplot(132)
plt.imshow(model[:,int(nx/2),:],aspect="auto",cmap="Greys",vmin=1500,vmax=5000)
cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.subplot(133)
plt.imshow(model[int(nz/2),:,:],aspect="auto",cmap="Greys",vmin=1500,vmax=5000)
cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.tight_layout()
plt.show()

