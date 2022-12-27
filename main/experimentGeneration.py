import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------------------

def buildModel3D(nx,ny,nz,property,z):
    '''


    '''

    model = np.zeros((nz,nx,ny), dtype=float)

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1, len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model

def createGaussianSurface(nx,ny,dx,dy,A,xc,yc,sigx,sigy):
    '''

    '''

    x, y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)

    surface = A*np.exp(-((x - xc)**2/(2*sigx**2) + (y - yc)**2/(2*sigy**2)))

    return surface.T

def analiticalRefraction(v,z,x):
    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

#-------------------------------------------------------------------------------------------

ns = 34
nr = 14

sx = np.linspace(50,9950,ns)
sy = np.linspace(50,9950,ns)

rx = np.linspace(1100,8900,nr)
ry = np.linspace(1100,8900,nr)

sx,sy = np.meshgrid(sx,sy)
rx,ry = np.meshgrid(rx,ry)

rx = np.reshape(rx, nr**2)
ry = np.reshape(ry, nr**2)

sx = np.reshape(sx, ns**2)
sy = np.reshape(sy, ns**2)

cmpx = np.zeros((nr*ns)**2)
cmpy = np.zeros((nr*ns)**2)

offset = np.zeros((nr*ns)**2)

k = 0
for shot in range(ns**2):
    for rec in range(nr**2):
        cmpx[k] = 0.5 * (rx[rec] - sx[shot]) + sx[shot]   
        cmpy[k] = 0.5 * (ry[rec] - sy[shot]) + sy[shot]

        offset[k] = np.sqrt((sx[shot] - rx[rec])**2 + (sy[shot] - ry[rec])**2)

        k += 1

plt.figure(1, figsize = (8,6))
plt.scatter(sx, sy, s = 5)
plt.scatter(rx, ry, s = 10)
plt.scatter(cmpx,cmpy, s = 1)
plt.xlim([0,1e4])
plt.ylim([0,1e4])

plt.tight_layout()
plt.show(block = False)

#-------------------------------------------------------------------------------------------

v = np.array([1500,1700,2000])
z = np.array([200, 1050]) 

x = np.linspace(np.min(offset), np.max(offset), 101)

td, tr = analiticalRefraction(v,z,x)

plt.figure(2, figsize = (5,8))

plt.plot(td)
for i in range(len(tr)):
    plt.plot(tr[i])

plt.gca().invert_yaxis()
plt.tight_layout()
plt.show(block = False)

#-------------------------------------------------------------------------------------------

dx = dy = dz = 12.5
nx, ny, nz = 801, 801, 121

z = np.array([200,1250,nz*dz])

model = buildModel3D(nx,ny,nz,v,z//dz)

A = np.array([800, 700, -650, 800, 900])
xc = np.array([3500, 3500, 5000, 6500, 6500])
yc = np.array([3500, 6500, 5000, 3500, 6500])
sigx = np.array([2000, 1500, 1000, 1000, 2000])
sigy = np.array([2000, 1500, 1000, 1000, 2000])

surface = np.zeros_like(model[0,:,:])
for k in range(len(A)):
    surface += createGaussianSurface(nx,ny,dx,dy,A[k],xc[k],yc[k],sigx[k],sigy[k])

surface = 24 + np.array((np.max(surface)-surface)//dz, dtype=int)

plt.figure(3, figsize = (8,6))
plt.contour(surface, levels = 20)
plt.imshow(surface, aspect = "auto", cmap="Greys")

plt.tight_layout()
plt.show(block = False)

for j in range(nx):
    for k in range(ny):
        model[surface[j,k]:,j,k] = v[-1]

plt.figure(4, figsize = (10,5))
plt.subplot(131)
plt.imshow(model[40,:,:])

plt.subplot(132)
plt.imshow(model[:,:,500])

plt.subplot(133)
plt.imshow(model[:,300,:])

plt.tight_layout()
plt.show()

