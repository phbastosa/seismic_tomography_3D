import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

nx = 101
ny = 101
nz = 101
dh = 10.0
nabc = 50
parameter = 0.0045

nxx = nx + 2*nabc
nyy = ny + 2*nabc
nzz = nz + 2*nabc

dampZX = np.ones((nzz,nxx))
dampZY = np.ones((nzz,nyy))
dampXY = np.ones((nxx,nyy))

damp3D = np.ones((nzz,nxx,nyy))

factor = np.zeros(nabc)

for k in range(nabc):   
    factor[k] = np.exp(-pow(parameter*(nabc-k),2))

# Left - Right damping
for i in range(nabc,nzz-nabc):
    dampZX[i,:nabc] = factor
    dampZX[i,nxx-nabc:nxx] = factor[::-1]
    
    dampZY[i,:nabc] = factor
    dampZY[i,nyy-nabc:nyy] = factor[::-1]

for i in range(nabc,nxx-nabc):
    dampXY[i,:nabc] = factor
    dampXY[i,nyy-nabc:nyy] = factor[::-1]

# Up - Bottom damping
for j in range(nabc,nxx-nabc):
    dampZX[:nabc,j] = factor
    dampZX[nzz-nabc:nzz,j] = factor[::-1]
    
for k in range(nabc,nyy-nabc):    
    dampZY[:nabc,k] = factor
    dampZY[nzz-nabc:nzz,k] = factor[::-1]

    dampXY[:nabc,k] = factor
    dampXY[nxx-nabc:nxx,k] = factor[::-1]

for i in range(nabc):
    # Up left corner
    dampZX[i:nabc,i] = factor[i]
    dampZY[i:nabc,i] = factor[i]
    dampXY[i:nabc,i] = factor[i]
    
    dampZX[i,i:nabc] = factor[i]    
    dampZY[i,i:nabc] = factor[i]
    dampXY[i,i:nabc] = factor[i]

    # Up right corner
    dampZX[i:nabc,nxx-i-1] = factor[i]
    dampZY[i:nabc,nyy-i-1] = factor[i]
    dampXY[nxx-i-1,i:nabc] = factor[i]
    dampZX[i,nxx-nabc:nxx-i-1] = factor[i]
    dampZY[i,nyy-nabc:nyy-i-1] = factor[i]
    dampXY[nxx-nabc:nxx-i-1,i] = factor[i]

    # Bottom left
    dampZX[nzz-i-1,i:nabc] = factor[i]    
    dampZY[nzz-i-1,i:nabc] = factor[i]
    dampXY[i:nabc,nyy-i-1] = factor[i]
    dampZX[nzz-nabc-1:nzz-i-1,i] = factor[i]
    dampZY[nzz-nabc-1:nzz-i-1,i] = factor[i]
    dampXY[i,nyy-nabc-1:nyy-i-1] = factor[i]

    # Bottom right
    dampZX[nzz-nabc-1:nzz-i,nxx-i-1] = factor[i]
    dampZX[nzz-i-1,nxx-nabc-1:nxx-i] = factor[i]   
    dampZY[nzz-nabc-1:nzz-i,nyy-i-1] = factor[i]
    dampZY[nzz-i-1,nyy-nabc-1:nyy-i] = factor[i]
    dampXY[nxx-i-1,nyy-nabc-1:nyy-i] = factor[i]
    dampXY[nxx-nabc-1:nxx-i,nyy-i-1] = factor[i]

    # 3D corners
    damp3D[i:nabc,i:nabc,i] = factor[i]
    damp3D[i:nabc,i,i:nabc] = factor[i]
    damp3D[i,i:nabc,i:nabc] = factor[i]
    damp3D[i:nabc,i:nabc,nyy-i-1] = factor[i]
    damp3D[i:nabc,nxx-i-1,i:nabc] = factor[i]
    damp3D[nzz-i-1,i:nabc,i:nabc] = factor[i]
    damp3D[i:nabc,nxx-nabc-1:nxx-i,i] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,i:nabc,i] = factor[i]
    damp3D[i:nabc,i,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,i,i:nabc] = factor[i]
    damp3D[i,i:nabc,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[i,nxx-nabc-1:nxx-i,i:nabc] = factor[i]
    damp3D[i:nabc,nxx-nabc-1:nxx-i,nyy-i-1] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,i:nabc,nyy-i-1] = factor[i]    
    damp3D[i:nabc,nxx-i-1,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,nxx-i-1,i:nabc] = factor[i]
    damp3D[nzz-i-1,i:nabc,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[nzz-i-1,nxx-nabc-1:nxx-i,i:nabc] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,nxx-nabc-1:nxx-i,i] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,i,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[i,nxx-nabc-1:nxx-i,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,nxx-nabc-1:nxx-i,nyy-i-1] = factor[i]
    damp3D[nzz-nabc-1:nzz-i,nxx-i-1,nyy-nabc-1:nyy-i] = factor[i]
    damp3D[nzz-i-1,nxx-nabc-1:nxx-i,nyy-nabc-1:nyy-i] = factor[i]

# Prisms XZ
for i in range(nabc,nyy-nabc):
    damp3D[:,:,i] = dampZX

# Prisms YZ
for i in range(nabc,nxx-nabc):
    damp3D[:,i,:] = dampZY

# Prisms XY
for i in range(nabc,nzz-nabc):
    damp3D[i,:,:] = dampXY

damp3D.flatten("F").astype("float32",order="F").tofile("inputs/models/cerjanVolumeABC.bin")

xzPlane = int(nyy/2)
yzPlane = int(nxx/2)
xyPlane = int(nzz/2)

vmin = 0
vmax = 1 

xloc = np.linspace(0,nxx-1,7,dtype=int)
xlab = np.array(xloc * dh,dtype=int)

yloc = np.linspace(0,nyy-1,7,dtype=int)
ylab = np.array(yloc * dh,dtype=int)

zloc = np.linspace(0,nzz-1,7,dtype=int)
zlab = np.array(zloc * dh,dtype=int)

# True model plot ------------------------------------------------------------------
plt.figure(1,figsize=(18, 6))

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
plt.imshow(damp3D[:,:,xzPlane],aspect="auto",cmap="Greys")
plt.title(f"XZ plane")
plt.xlabel("X axis [m]")
plt.ylabel("Z axis [m]")

plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)

axes_2 = plt.subplot(G[1:,:3])
plt.imshow(damp3D[:,yzPlane,:],aspect="auto",cmap="Greys")
plt.title(f"YZ plane")
plt.xlabel("Y axis [m]")
plt.ylabel("Z axis [m]")

plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)

axes_3 = plt.subplot(G[:,3:])
plt.imshow(damp3D[xyPlane,:,:].T,aspect="auto",cmap="Greys")
plt.title(f"XY plane")
plt.xlabel("X axis [m]")
plt.ylabel("Y axis [m]")

cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

plt.xticks(xloc,xlab)
plt.yticks(yloc,ylab)

plt.gca().invert_yaxis()

plt.tight_layout()

plt.show()