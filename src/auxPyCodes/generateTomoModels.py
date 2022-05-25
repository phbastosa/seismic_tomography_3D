import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.ndimage import gaussian_filter

def buildModel3D(nx,ny,nz,property,z):
    model = np.zeros((nz,nx,ny), dtype=float)

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1,len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model

def buildGradientModel3D(nx,ny,nz,dh,vi,zi,a,b,c):
    
    model = np.zeros((nz,nx,ny))

    wb = int(zi/dh)  # water bottom

    model[:wb,:,:] = vi

    for i in range(nz):
        if i >= wb:
            for j in range(nx):
                for k in range(ny):
                    model[i,j,k] = vi + a*(i*dh - zi) + b*dh*j + c*k*dh 

    return model

def createGaussianSurface(nx,ny,dx,dy,A,xc,yc,sigx,sigy):
    
    x,y = np.meshgrid(np.arange(nx)*dx,np.arange(ny)*dy)

    surface = A*np.exp(-((x - xc)**2/(2*sigx**2) + (y - yc)**2/(2*sigy**2)))

    return surface.T

nx = 101
ny = 101
nz = 31

dh = 50.0

v = np.array([1500, 1700, 1900, 2100, 2300, 2500, 2700, 2900, 3500],dtype=float)
z = np.array([ 200,  400,  600,  800, 1000, 1200, 1400, 1450,nz*dh],dtype=float) // dh

initModel = 1.0 / gaussian_filter(1.0 / buildModel3D(nx,ny,nz,v,z), 2.0)

vi = 1500
zi = 200

a = 0.8
b = 0.0
c = 0.0

trueModel = buildGradientModel3D(nx,ny,nz,dh,vi,zi,a,b,c)

A = -800

xc = 2500
yc = 2500

sigx = 1000
sigy = 1000

surface = createGaussianSurface(nx,ny,dh,dh,A,xc,yc,sigx,sigy) // dh + 29

for j in range(nx):
    for k in range(ny):
        for i in range(nz):
            if i > surface[j,k]:
                trueModel[i,j,k] = 3500.0

trueModel = 1.0 / gaussian_filter(1.0 / trueModel, 1.5)

trueModel.flatten("F").astype("float32",order="F").tofile(f"inputs/models/trueModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
initModel.flatten("F").astype("float32",order="F").tofile(f"inputs/models/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

xzPlane = int(nx/2)
yzPlane = int(ny/2)
xyPlane = int(nz/2)

vmin = np.min(trueModel)
vmax = np.max(trueModel) 

xloc = np.linspace(0,nx-1,7,dtype=int)
xlab = np.array(xloc * dh,dtype=int)

yloc = np.linspace(0,ny-1,7,dtype=int)
ylab = np.array(yloc * dh,dtype=int)

zloc = np.linspace(0,nz-1,7,dtype=int)
zlab = np.array(zloc * dh,dtype=int)

# True model plot ------------------------------------------------------------------
plt.figure(1,figsize=(18, 6))

plt.suptitle("Benchmark model for the first-arrival tomography")

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
plt.imshow(trueModel[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"XZ plane")
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.xlabel("X axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

axes_2 = plt.subplot(G[1:,:3])
plt.imshow(trueModel[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"YZ plane")
plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.xlabel("Y axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

axes_3 = plt.subplot(G[:,3:])
plt.imshow(trueModel[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
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

# Init model plot ------------------------------------------------------------------
plt.figure(2,figsize=(18, 6))

plt.suptitle("Benchmark model for the first-arrival tomography")

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
plt.imshow(initModel[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"XZ plane")
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.xlabel("X axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

axes_2 = plt.subplot(G[1:,:3])
plt.imshow(initModel[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.title(f"YZ plane")
plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.xlabel("Y axis [m]")
plt.ylabel("Z axis [m]")

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

axes_3 = plt.subplot(G[:,3:])
plt.imshow(initModel[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
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

plt.show()