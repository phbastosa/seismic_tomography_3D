import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def readBinaryVolume(dim1,dim2,dim3,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
    return np.reshape(data, [dim1,dim2,dim3], order='F')

def readBinaryArray(dim, filename):
    return np.fromfile(filename, dtype=np.int32, count=dim)

nx = 101
ny = 101
nz = 31

dh = 50.0

eikonalShot = 0

model = readBinaryVolume(nz,nx,ny,f"../inputs/models/trueModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
# times = readBinaryVolume(nz,nx,ny,f"../outputs/travelTimes/eikonalTimes_shot_{eikonalShot+1}.bin")

shots = np.loadtxt(f"../outputs/geometry/shotsPosition.txt", delimiter = ",")
nodes = np.loadtxt(f"../outputs/geometry/nodesPosition.txt", delimiter = ",")

sx = shots[:,0]
sy = shots[:,1]
sz = shots[:,2]

rx = nodes[:,0]
ry = nodes[:,1]
rz = nodes[:,2]

xzPlane = int(ry[eikonalShot] / dh)
yzPlane = int(rx[eikonalShot] / dh)
xyPlane = int(rz[eikonalShot] / dh)

vmin = np.min(model)
vmax = np.max(model) 

xloc = np.linspace(0,nx-1,7,dtype=int)
xlab = np.array(xloc * dh,dtype=int)

yloc = np.linspace(0,ny-1,7,dtype=int)
ylab = np.array(yloc * dh,dtype=int)

zloc = np.linspace(0,nz-1,7,dtype=int)
zlab = np.array(zloc * dh,dtype=int)

plt.figure(1,figsize=(18, 6))

plt.suptitle("Benchmark model for the first arrival tomography")

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
# plt.contour(times[:,xzPlane,:], levels = 25)
plt.imshow(model[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"XZ plane")
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.xlabel("X axis [m]")
plt.ylabel("Z axis [m]")

plt.scatter(sx/dh,sz/dh)
plt.scatter(rx/dh,rz/dh)

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

axes_2 = plt.subplot(G[1:,:3])
# plt.contour(times[:,:,yzPlane], levels = 25)
plt.imshow(model[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"YZ plane")
plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.xlabel("Y axis [m]")
plt.ylabel("Z axis [m]")

plt.scatter(sy/dh,sz/dh)
plt.scatter(ry/dh,rz/dh)

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

axes_3 = plt.subplot(G[:,3:])
# plt.contour(times[xyPlane,:,:], levels = 25)
plt.imshow(model[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)

plt.title(f"XY plane")
plt.xticks(xloc,xlab)
plt.yticks(yloc,ylab)
plt.xlabel("X axis [m]")
plt.ylabel("Y axis [m]")

cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.scatter(sx/dh,sy/dh)
plt.scatter(rx/dh,ry/dh)

plt.plot(np.arange(nx),np.ones(nx)*xzPlane,'--r', markersize = 10)
plt.plot(np.ones(ny)*yzPlane,np.arange(ny),'--r', markersize = 10)

plt.gca().invert_yaxis()

plt.tight_layout()

plt.show()
