import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

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

    return surface

def analiticalRefraction(v,z,x):
    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

def boxPlot(model, dx, dy, dz, xySlice, zySlice, zxSlice):
    
    xyPlane = model[xySlice,:,:].T
    zxPlane = model[:,:,zxSlice]
    zyPlane = model[:,zySlice,:].T

    ticks = np.array([3,7,7], dtype = int)

    maxModelDistance = np.max(np.shape(model))
    minModelDistance = np.min(np.shape(model))

    #------------------------------------------------
    modelShape = np.array(np.shape(model))
    [nz,nx,ny] = modelShape
    [z, x, y] = 3.0 * (minModelDistance / maxModelDistance) * modelShape / maxModelDistance

    vmin = np.min(model)
    vmax = np.max(model)

    px = 1/plt.rcParams['figure.dpi']  
    fig = plt.figure(1, figsize=(1000*px, 800*px))

    xloc = np.linspace(0,nx,ticks[1], dtype = int)
    yloc = np.linspace(0,ny,ticks[2], dtype = int)
    zloc = np.linspace(0,nz,ticks[0], dtype = int)

    xlab = np.around(xloc * dx / 1000, decimals = 1)
    ylab = np.around(yloc * dy / 1000, decimals = 1)
    zlab = np.around(zloc * dz / 1000, decimals = 1)

    #------------------------------------------------
    ax1 = fig.add_axes([0.75 - x, 0.98 - y, x, y])
    ax1.imshow(xyPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax1.tick_params(direction = 'in', axis='x') 
    ax1.tick_params(direction = 'in', axis='y') 
    ax1.set_xticklabels([])

    ax1.plot(np.arange(modelShape[1]), np.ones(modelShape[1])*zxSlice, "--g")
    ax1.plot(np.ones(modelShape[2])*zySlice, np.arange(modelShape[2]), "--m")

    ax1.set_xticks(xloc)
    ax1.set_yticks(yloc[1:])
    ax1.set_yticklabels(ylab[1:])

    ax1.grid(color='w', linestyle='--', linewidth=0.3)

    ax1.set_ylabel("Y [km]")
    ax1.invert_yaxis()

    #------------------------------------------------
    ax2 = fig.add_axes([0.75, 0.98 - y, z, y])
    ax2.imshow(zyPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax2.tick_params(direction = 'out', axis='x') 
    ax2.tick_params(direction = 'in', axis='y') 
    ax2.set_yticklabels([])

    ax2.set_xlabel("Z [km]")

    ax2.plot(np.arange(modelShape[0]), np.ones(modelShape[0])*zxSlice, "--g")
    ax2.plot(np.ones(modelShape[2])*xySlice, np.arange(modelShape[2]), "--r")

    ax2.set_yticks(yloc)

    ax2.set_xticks(zloc[1:])
    ax2.set_xticklabels(zlab[1:])

    ax2.grid(color = 'w', linestyle = '--', linewidth = 0.3)

    #------------------------------------------------
    ax3 = fig.add_axes([0.75 - x, 0.98 - y - z  , x, z])
    im = ax3.imshow(zxPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)

    ax3.plot(np.arange(modelShape[1]), np.ones(modelShape[1])*xySlice, "--r")
    ax3.plot(np.ones(modelShape[0])*zySlice, np.arange(modelShape[0]), "--m")

    ax3.set_xticks(xloc)
    ax3.set_yticks(zloc)

    ax3.set_xticklabels(xlab)
    ax3.set_yticklabels(zlab)

    ax3.grid(color = 'w', linestyle = '--', linewidth = 0.3)

    ax3.set_xlabel("X [km]")
    ax3.set_ylabel("Z [km]")

    #------------------------------------------------
    ax4 = fig.add_axes([0.75 - x, 0.98 - y - 2*z, x, z])

    ax4.axis("off")

    cmap = mpl.colormaps["Greys"]
    norm = mpl.colors.Normalize(vmin*1e-3, vmax*1e-3)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("bottom", size="10%", pad=0)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(vmin*1e-3, vmax*1e-3, 5), orientation = "horizontal")
    cbar.ax.set_xticklabels(np.around(np.linspace(vmin*1e-3, vmax*1e-3, 5), decimals = 1))
    cbar.set_label("P wave velocity [km/s]")
    
    return None

def multiBoxPlot(models:np.ndarray, dh:np.ndarray, slices:np.ndarray, subplots:np.ndarray) -> None:

    modelShape = np.array(np.shape(models[0]))
    maxModelDistance = np.max(np.shape(models[0]))
    minModelDistance = np.min(np.shape(models[0]))
    
    [z, x, y] = 3.0 * (minModelDistance / maxModelDistance) * modelShape / maxModelDistance

    vmin = np.min(models[0])
    vmax = np.max(models[0])

    px = 1/plt.rcParams['figure.dpi']  
    ticks = np.array([2,7,7], dtype = int)

    fig = plt.figure(1, figsize=(600*px*subplots[1], 500*px*subplots[0]))

    xloc = np.linspace(0,nx-1,ticks[1], dtype = int)
    yloc = np.linspace(0,ny-1,ticks[2], dtype = int)
    zloc = np.linspace(0,nz-1,ticks[0], dtype = int)

    m2km = 1e-3

    xlab = np.around(xloc * dx * m2km, decimals = 1)
    ylab = np.around(yloc * dy * m2km, decimals = 1)
    zlab = np.around(zloc * dz * m2km, decimals = 1)

    axes = np.array([[0.75 - x, 0.98 - y      , x, y], 
                     [    0.75, 0.98 - y      , z, y],
                     [0.75 - x, 0.98 - y - z  , x, z],
                     [0.75 - x, 0.98 - y - 2.5*z, x, z]])

    xTickDirection = ['out', 'out', 'out']
    yTickDirection = ['in', 'in', 'in']

    xTickLock = [xloc, zloc[1:], xloc]
    yTickLock = [yloc, yloc, zloc[1:]]

    xTickLabel = [[], zlab[1:], xlab]
    yTickLabel = [ylab, [], zlab[1:]]

    xLabel = ["X [km]", "Z [km]", "X [km]"]
    yLabel = ["Y [km]", "      ", "Z [km]"]

    yInvert = [ True, False, False]

    xSlices = [[np.arange(modelShape[1]), np.ones(modelShape[1])*slices[1], "--g"],
               [np.arange(modelShape[0]), np.ones(modelShape[0])*slices[1], "--g"],
               [np.arange(modelShape[1]), np.ones(modelShape[1])*slices[0], "--r"]] 

    ySlices = [[np.ones(modelShape[2])*slices[2], np.arange(modelShape[2]), "--m"],
               [np.ones(modelShape[2])*slices[0], np.arange(modelShape[2]), "--r"],
               [np.ones(modelShape[0])*slices[2], np.arange(modelShape[0]), "--m"]]

    #--------------------------------------------------------------------------------    

    subfigs = fig.subfigures(subplots[0], subplots[1])
    
    for i in range(subplots[0]):
        for j in range(subplots[1]):

            ind = i*subplots[0] + j 

            ims = [models[ind, slices[0],:,:].T, models[ind,:,slices[2],:].T, models[ind,:,:,slices[1]]]
            
            for k, axs in enumerate(axes):

                # Adjusting acording subplot size      
                if subplots[0] == 1:
                    if subplots[1] == 1:
                        ax = subfigs.add_axes(axs)                         
                    else:
                        ax = subfigs[j].add_axes(axs)

                elif subplots[1] == 1:
                    if subplots[0] == 1:
                        ax = subfigs.add_axes(axs)        
                    else:    
                        ax = subfigs[i].add_axes(axs)
                
                else:
                    ax = subfigs[i,j].add_axes(axs)

                # Setting colorbar
                if k == 3:

                    ax.axis("off")

                    cmap = mpl.colormaps["Greys"]
                    norm = mpl.colors.Normalize(vmin*1e-3, vmax*1e-3)
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("bottom", size="10%", pad=0)
                    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(vmin*1e-3, vmax*1e-3, 5), orientation = "horizontal")
                    cbar.ax.set_xticklabels(np.around(np.linspace(vmin*1e-3, vmax*1e-3, 5), decimals = 1))
                    cbar.set_label("Velocidade [km/s]")
                
                # plotting model slices 
                else:
                    
                    ax.imshow(ims[k], aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
                    
                    ax.tick_params(direction = xTickDirection[k], axis='x') 
                    ax.tick_params(direction = yTickDirection[k], axis='y') 
                    
                    ax.set_xticks(xTickLock[k])
                    ax.set_yticks(yTickLock[k])

                    ax.set_xticklabels(xTickLabel[k])
                    ax.set_yticklabels(yTickLabel[k])
 
                    ax.set_xlabel(xLabel[k])
                    ax.set_ylabel(yLabel[k])

                    ax.plot(xSlices[k][0], xSlices[k][1], xSlices[k][2])
                    ax.plot(ySlices[k][0], ySlices[k][1], ySlices[k][2])
                    
                    if yInvert[k]:
                       ax.invert_yaxis()
    
    return None

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

# Geometry plot
plt.figure(10, figsize = (8,6))
plt.scatter(sx, sy, s = 5, label = "Receptores")
plt.scatter(rx, ry, s = 10, label = "Fontes")
plt.scatter(cmpx,cmpy, s = 1, label = "Pontos médios")
plt.xlim([0,1e4])
plt.ylim([0,1e4])

plt.title(f"Geometria com {ns**2} fontes e {nr**2} receptores", fontsize = 18)
plt.xlabel("X [m]", fontsize = 15)
plt.ylabel("Y [m]", fontsize = 15)

plt.legend(loc = "lower left", fontsize = 10)

plt.tight_layout()

plt.savefig("../figures/0_geometry.png", dpi = 200)
plt.show(block=False)

#-------------------------------------------------------------------------------------------

v = np.array([1500,1700,2000])
z = np.array([200, 1050]) 

x = np.linspace(np.min(offset), np.max(offset), 101)

td, tr = analiticalRefraction(v,z,x)

# Initial travel time plot 
plt.figure(11, figsize = (10,6))

plt.subplot(121)

vp = np.zeros(1500)
vp[:200] = 1500
vp[200:1250] = 1700
vp[1250:] = 2000

depth = np.arange(1500)

plt.xlim([1000,2500])
plt.ylim([0,1500])

plt.plot(vp, depth)
plt.gca().invert_yaxis()

plt.yticks(np.linspace(0,1500,11), np.linspace(0,1500,11, dtype = int))
plt.xticks(np.linspace(1000,2500,7), np.linspace(1000,2500,7)/1000)

plt.title("Modelo de velocidade inicial", fontsize = 18)
plt.xlabel("Velocidade [km/s]", fontsize = 15)
plt.ylabel("Profundidade [m]", fontsize = 15)

plt.subplot(122)
plt.plot(x, td)
for i in range(len(tr)):
    plt.plot(x, tr[i])

plt.xlim([x[0],x[-1]])
plt.ylim([0, 8])

plt.title("Tempo de trânsito analítico", fontsize = 18)
plt.xlabel("Offset [m]", fontsize = 15)
plt.ylabel("Tempo [s]", fontsize = 15)

plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("../figures/1_initialTravelTimes.png", dpi = 200)
plt.show(block=False)
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

i,j = np.where(-surface*dz < -1250)
surface[i,j] = 1250/dz

# Surface plot
fig, ax = plt.subplots(1,1, figsize = (8,8))

cmap = mpl.colormaps["coolwarm"]

ax.contour(-dz*surface, levels = 10)
ax.imshow(-dz*surface, aspect = "auto", cmap=cmap)

ax.invert_yaxis()

ax.set_title("Superfície gaussiana", fontsize = 18)
ax.set_xlabel("X [m]", fontsize = 15)
ax.set_ylabel("Y [m]", fontsize = 15)

smin = np.min(-dz*surface)
smax = np.max(-dz*surface)

norm = mpl.colors.Normalize(smin,smax)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.6)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(smin, smax, 5), orientation = "horizontal")
cbar.ax.set_xticklabels(np.around(np.linspace(smin, smax, 5), decimals = 1))
cbar.set_label("Profundidade [m]", fontsize = 15)

plt.tight_layout()
plt.savefig("../figures/2_gaussianSurface.png", dpi = 200)
plt.show(block=False)

#-------------------------------------------------------------------------------------------
initModel = model.copy()

for j in range(nx):
    for k in range(ny):
        model[surface[j,k]:,j,k] = v[-1]

models = np.zeros((2,nz,nx,ny))

models[0,:,:,:] = model
models[1,:,:,:] = initModel

dh = np.array([dz, dx, dy])
slices = np.array([int(nz / 2), int(ny / 2), int(nx / 2)], dtype = int) # [xy, zx, zy]
subplots = np.array([1, 2], dtype = int)

multiBoxPlot(models, dh, slices, subplots)
plt.savefig("../figures/3_models.png", dpi = 200)
plt.show(block=False)

