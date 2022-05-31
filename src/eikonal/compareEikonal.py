import sys
import numpy as np
import matplotlib.pyplot as plt

def readBinaryVolume(dim1,dim2,dim3,type,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order=type)
    return volume

def readBinaryVolume2(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

nx = 221
ny = 221
nz = 12

dh = 100

noble = readBinaryVolume(nz,nx,ny,"C","noble_travelTimes3D.bin")
jeong = readBinaryVolume2(nz,nx,ny,"fim_central_eikonal_nz12_nx221_ny221_shot_1.bin")
podvin = readBinaryVolume2(nz,nx,ny,"pod_central_eikonal_nz12_nx221_ny221_shot_1.bin")


plt.figure(1)
plt.subplot(311)
plt.contour(noble[:,:,int(ny/2)], levels=20)
plt.imshow(noble[:,:,int(ny/2)] - podvin[:,:,int(ny/2)], aspect="auto", cmap="Greys")

plt.subplot(312)
plt.contour(podvin[:,:,int(ny/2)], levels=20)
plt.imshow(podvin[:,:,int(ny/2)], aspect="auto", cmap="Greys")

plt.subplot(313)
plt.contour(jeong[:,:,int(ny/2)], levels=20)
plt.imshow(noble[:,:,int(ny/2)] - jeong[:,:,int(ny/2)], aspect="auto", cmap="Greys")

plt.show()