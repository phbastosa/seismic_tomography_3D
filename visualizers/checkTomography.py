import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.ndimage import gaussian_filter

def readBinaryVolume(dim1,dim2,dim3,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
    return np.reshape(data, [dim1,dim2,dim3], order='F')

def readBinaryArray(dim, filename):
    return np.fromfile(filename, dtype=np.int32, count=dim)

nx = 101
ny = 101
nz = 31

dh = 50.0

# model = readBinaryVolume(nz,nx,ny,f"inputs/models/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

# sx,sy,sz = np.loadtxt(f"outputs/geometry/shotsPosition.txt", delimiter = ",",unpack=True)
# rx,ry,rz = np.loadtxt(f"outputs/geometry/nodesPosition.txt", delimiter = ",",unpack=True)

# xzPlane = int(ny / 2)
# yzPlane = int(nx / 2)
# xyPlane = int(nz / 2)

# vmin = np.min(model)
# vmax = np.max(model) 

# xloc = np.linspace(0,nx-1,7,dtype=int)
# xlab = np.array(xloc * dh,dtype=int)

# yloc = np.linspace(0,ny-1,7,dtype=int)
# ylab = np.array(yloc * dh,dtype=int)

# zloc = np.linspace(0,nz-1,7,dtype=int)
# zlab = np.array(zloc * dh,dtype=int)

# plt.figure(1,figsize=(15, 5))

# plt.suptitle("Initial model", fontsize = 20)

# G = gridspec.GridSpec(2, 5)

# axes_1 = plt.subplot(G[:1,:3])
# plt.imshow(model[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
# plt.xticks(xloc,xlab)
# plt.yticks(zloc,zlab)
# plt.title(f"XZ plane", fontsize = 18)
# plt.xlabel("X axis [m]", fontsize = 15)
# plt.ylabel("Z axis [m]", fontsize = 15)

# plt.scatter(sx/dh,sz/dh, label = "Shots")
# plt.scatter(rx/dh,rz/dh, label = "Nodes") 

# plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
# plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

# plt.legend(loc = "lower left", fontsize = 10)

# axes_2 = plt.subplot(G[1:,:3])
# plt.imshow(model[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
# plt.xticks(yloc,ylab)
# plt.yticks(zloc,zlab)
# plt.title(f"YZ plane", fontsize = 18)
# plt.xlabel("Y axis [m]", fontsize = 15)
# plt.ylabel("Z axis [m]", fontsize = 15)

# plt.scatter(sy/dh,sz/dh, label = "Shots")
# plt.scatter(ry/dh,rz/dh, label = "Nodes")

# plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
# plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

# plt.legend(loc = "lower left", fontsize = 10)

# axes_3 = plt.subplot(G[:,3:])
# plt.imshow(model[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)

# plt.xticks(xloc,xlab)
# plt.yticks(yloc,ylab)
# plt.title(f"XY plane", fontsize = 18)
# plt.xlabel("X axis [m]", fontsize = 15)
# plt.ylabel("Y axis [m]", fontsize = 15)

# cbar = plt.colorbar()
# cbar.set_label("P wave velocity [m/s]",fontsize=15)

# plt.scatter(sx/dh,sy/dh, s = 5, label = "Shots")
# plt.scatter(rx/dh,ry/dh, s = 10, label = "Nodes")

# plt.plot(np.arange(nx),np.ones(nx)*xzPlane,'--r', markersize = 10)
# plt.plot(np.ones(ny)*yzPlane,np.arange(ny),'--r', markersize = 10)

# plt.legend(loc = "lower left", fontsize = 10)

# plt.gca().invert_yaxis()

# plt.tight_layout()

# plt.show()

trueModel = readBinaryVolume(nz,nx,ny,f"inputs/models/trueModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
initModel = readBinaryVolume(nz,nx,ny,f"inputs/models/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
predModel = readBinaryVolume(nz,nx,ny,f"outputs/models/estimatedModel_iteration_1.bin")

trueModel = trueModel[:,int(nx/2),int(ny/2)]
initModel = initModel[:,int(nx/2),int(ny/2)]
predModel = predModel[:,int(nx/2),int(ny/2)]

predModel = 1.0 / gaussian_filter(1.0 / predModel, 2)

depth = np.arange(nz) * dh

plt.figure(2, figsize=(5,10))
plt.plot(trueModel, depth)
plt.plot(initModel, depth)
plt.plot(predModel, depth)

plt.ylim([0,(nz-1)*dh])
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
