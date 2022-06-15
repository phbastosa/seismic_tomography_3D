import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def analiticalRefraction2D(v,z,x):
    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

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
nz = 101

nt = 1001
dt = 1e-3

dh = 10.0
recSpacing = 50.0

vp = np.array([1500, 1600, 1700, 1800, 2000])

vs = vp / np.sqrt(3.0); vs[0] = 0.0 
rho = 310 * vp ** 0.25; rho[0] = 1000.0 

z = np.array([200, 200, 200, 200])
depth = np.array([ 200, 400, 600, 800, nz*dh],dtype=float) // dh

x = np.arange(0,nx*dh,recSpacing)

print(len(x))

td, t = analiticalRefraction2D(vp,z,x)

fb = np.zeros(len(x))

for i in range(len(x)):
    times = np.append(t[:,i], td[i])
    fb[i] = np.min(times)

maxDepth = int(depth[-1]*dh)
wellDepth = np.arange(maxDepth)
vpWell = np.ones(maxDepth) * vp[0]
vsWell = np.ones(maxDepth) * vs[0]
rhoWell = np.ones(maxDepth) * rho[0]

for i in range(1,len(depth)):
    var = slice(int(depth[i-1]*dh),int(depth[i]*dh))
    vpWell[var] = np.ones(int((depth[i]-depth[i-1])*dh)) * vp[i]
    vsWell[var] = np.ones(int((depth[i]-depth[i-1])*dh)) * vs[i]
    rhoWell[var] = np.ones(int((depth[i]-depth[i-1])*dh)) * rho[i]

vpModel = buildModel3D(nx,ny,nz,vp,depth)
vsModel = buildModel3D(nx,ny,nz,vs,depth)
rhoModel = buildModel3D(nx,ny,nz,rho,depth)

vpModel.flatten("F").astype("float32",order="F").tofile(f"inputs/models/vpModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
vsModel.flatten("F").astype("float32",order="F").tofile(f"inputs/models/vsModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
rhoModel.flatten("F").astype("float32",order="F").tofile(f"inputs/models/rhoModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")

xzPlane = int(ny/2)
yzPlane = int(nx/2)
xyPlane = int(nz/2)

xloc = np.linspace(0,nx-1,11,dtype=int)
xlab = np.array(xloc * dh / 1000,dtype=int)

yloc = np.linspace(0,ny-1,11,dtype=int)
ylab = np.array(yloc * dh / 1000,dtype=int)

zloc = np.linspace(0,nz-1,10,dtype=int)
zlab = np.around(zloc * dh / 1000, decimals = 1)

# True model plot ------------------------------------------------------------------
plt.figure(1,figsize=(15, 6))

G = gridspec.GridSpec(2, 7)

axes_1 = plt.subplot(G[:1,:3])
plt.title(f"2.5 D refractive geometry", fontsize = 20)
plt.xlabel("X axis [km]", fontsize = 15)
plt.ylabel("Y axis [m]", fontsize = 15)

xr = np.arange(0,nx*dh,recSpacing)
yr = np.ones(len(xr)) * (ny-1)*dh/2

plt.plot(xr,yr, "ok", markersize = 5, label = "Receivers position")
plt.plot(xr[0],yr[0], "*r", markersize = 10, label = "Shot position")

plt.xlim([0,(nx-1)*dh])
plt.ylim([0,(ny-1)*dh])
plt.xticks(xloc*dh,xlab)

plt.legend(loc = "upper right")

axes_2 = plt.subplot(G[1:,:3])
plt.imshow(vpModel[:,:,xzPlane],aspect="auto",cmap="Greys")
plt.title(f"Model section: xz plane", fontsize = 20)
plt.xlabel("X axis [km]", fontsize = 15)
plt.ylabel("Z axis [km]", fontsize = 15)
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)

axes_3 = plt.subplot(G[:,3:5])
plt.title("Property logs", fontsize=20)
plt.plot(vpWell, wellDepth, label = "P wave velocity [m/s]")
plt.plot(vsWell, wellDepth, label = "S wave velocity [m/s]")
plt.plot(rhoWell, wellDepth, label = "Density [kg/m³]")
plt.xlabel("Seismic property", fontsize = 15)
plt.ylabel("Depth [km]", fontsize = 15)
plt.legend(loc = "upper right")
plt.xlim([np.min(vs)-500,np.max(vp)+500])
plt.ylim([0,(nz-1)*dh])
plt.yticks(zloc*dh,zlab)
plt.gca().invert_yaxis()

axes_4 = plt.subplot(G[:,5:])

plt.plot(x,td)

# a = np.around(td,decimals=2)
# for i in range(len(z)):
#     b = np.around(t[i],decimals=1)

#     xc = np.where(a == b)[0][0]

#     plt.plot(x[xc:],t[i,xc:])

plt.plot(x,fb,"k")

plt.ylim([0,nt*dt])
plt.xlim([0,nx*dh])

plt.xticks(xloc*dh,xlab)
plt.title("Expected seismogram", fontsize = 20)
plt.xlabel("Offset [km]", fontsize = 15)
plt.ylabel("Two way time [s]", fontsize = 15)

plt.gca().invert_yaxis()

plt.tight_layout()
plt.show()