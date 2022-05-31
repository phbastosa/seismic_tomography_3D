import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def readBinaryVolume(dim1,dim2,dim3,type,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order=type)
    return volume

nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])

dx = float(sys.argv[4])
dy = float(sys.argv[5])
dz = float(sys.argv[6])

yzPlane = int(float(sys.argv[7]) / dx)
xzPlane = int(float(sys.argv[8]) / dy) 
xyPlane = int(float(sys.argv[9]) / dz)

model = readBinaryVolume(nz,ny,nx,"C",sys.argv[10])    
ttPaulo = readBinaryVolume(nz,ny,nx,"C",sys.argv[11])

vmin = np.min(model)
vmax = np.max(model) 

xloc = np.linspace(0,nx-1,7,dtype=int)
xlab = np.array(xloc * dx,dtype=int)

yloc = np.linspace(0,ny-1,7,dtype=int)
ylab = np.array(yloc * dy,dtype=int)

zloc = np.linspace(0,nz-1,7,dtype=int)
zlab = np.array(zloc * dz,dtype=int)

# Travel time with Podvin (1991) solver ------------------------------------------------------------------
plt.figure(2,figsize=(18, 6))

plt.suptitle("Travel times 3D - Fast Sweeping Method")

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
plt.contour(ttPaulo[:,xzPlane,:],levels=25)
plt.imshow(model[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"XZ plane")
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.xlabel("X axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

axes_2 = plt.subplot(G[1:,:3])
plt.contour(ttPaulo[:,:,yzPlane],levels=25)
plt.imshow(model[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"YZ plane")
plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.xlabel("Y axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

axes_3 = plt.subplot(G[:,3:])
plt.contour(ttPaulo[xyPlane,:,:],levels=20)
plt.imshow(model[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"XY plane")
plt.xticks(xloc,xlab)
plt.yticks(yloc,ylab)
plt.xlabel("X axis [m]")
plt.ylabel("Y axis [m]")

cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.plot(np.arange(nx),np.ones(nx)*xzPlane,'--r')
plt.plot(np.ones(ny)*yzPlane,np.arange(ny),'--r')

plt.gca().invert_yaxis()

plt.tight_layout()

plt.savefig("nobleFSM.png",dpi=200,bbox_inches="tight")

plt.show()