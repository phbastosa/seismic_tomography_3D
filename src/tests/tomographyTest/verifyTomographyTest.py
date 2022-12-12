import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from skimage import exposure 
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable

def perc(matrix,value):
    p = np.percentile(matrix,[.5, value])                     
    image = exposure.rescale_intensity(matrix, in_range=(p[0],p[1]), out_range=(0,255))                    

    return image       

def readBinaryVolume(dim1,dim2,dim3,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
    return np.reshape(data, [dim1,dim2,dim3], order='F')

def readBinaryArray(dim, filename):
    return np.fromfile(filename, dtype=np.int32, count=dim)

def boxPlot(model, dx, dy, dz, xySlice, zySlice, zxSlice, colorbarFix = 1.5):
    
    xyPlane = model[xySlice,:,:].T
    zxPlane = model[:,:,zxSlice]
    zyPlane = model[:,zySlice,:].T

    ticks = np.array([3,7,7], dtype = int)

    #------------------------------------------------
    axis = np.array(np.shape(model))
    [z, x, y] = axis * 0.6 / np.max(axis)

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
    ax1 = fig.add_axes([0.7 - x, 0.95 - y, 0 + x, 0 + y])
    ax1.imshow(xyPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax1.tick_params(direction = 'in', axis='x') 
    ax1.tick_params(direction = 'in', axis='y') 
    ax1.set_xticklabels([])

    ax1.plot(np.arange(axis[1]), np.ones(axis[1])*zxSlice, "--g")
    ax1.plot(np.ones(axis[2])*zySlice, np.arange(axis[2]), "--m")

    ax1.set_xticks(xloc)
    ax1.set_yticks(yloc[1:])
    ax1.set_yticklabels(ylab[1:])

    ax1.grid(color='w', linestyle='--', linewidth=0.3)

    ax1.set_ylabel("Y [km]")
    ax1.invert_yaxis()

    #------------------------------------------------
    ax2 = fig.add_axes([0.7, 0.95 - y, 0 + z, 0 + y])
    ax2.imshow(zyPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax2.tick_params(direction = 'out', axis='x') 
    ax2.tick_params(direction = 'in', axis='y') 
    ax2.set_yticklabels([])

    ax2.set_xlabel("Z [km]")

    ax2.plot(np.arange(axis[0]), np.ones(axis[0])*zxSlice, "--g")
    ax2.plot(np.ones(axis[2])*xySlice, np.arange(axis[2]), "--r")

    ax2.set_yticks(yloc)

    ax2.set_xticks(zloc[1:])
    ax2.set_xticklabels(zlab[1:])

    ax2.grid(color = 'w', linestyle = '--', linewidth = 0.3)

    #------------------------------------------------
    ax3 = fig.add_axes([0.7 - x, 0.95 - y - z, 0 + x, 0 + z])
    im = ax3.imshow(zxPlane, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)

    ax3.plot(np.arange(axis[1]), np.ones(axis[1])*xySlice, "--r")
    ax3.plot(np.ones(axis[0])*zySlice, np.arange(axis[0]), "--m")

    ax3.set_xticks(xloc)
    ax3.set_yticks(zloc)

    ax3.set_xticklabels(xlab)
    ax3.set_yticklabels(zlab)

    ax3.grid(color = 'w', linestyle = '--', linewidth = 0.3)

    ax3.set_xlabel("X [km]")
    ax3.set_ylabel("Z [km]")

    #------------------------------------------------
    ax4 = fig.add_axes([0.7 - x, 0.95 - y - colorbarFix*z, 0 + x, 0 + z])

    ax4.axis("off")

    cmap = mpl.colormaps["Greys"]
    norm = mpl.colors.Normalize(vmin*1e-3, vmax*1e-3)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("bottom", size="10%", pad=0)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(vmin*1e-3, vmax*1e-3, 5), orientation = "horizontal")
    cbar.ax.set_xticklabels(np.around(np.linspace(vmin*1e-3, vmax*1e-3, 5), decimals = 1))
    cbar.set_label("P wave velocity [km/s]")
    
    return None


nx = 101
ny = 101
nz = 31

dh = 50.0

# model = readBinaryVolume(nz,nx,ny,f"inputs/models/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
model = readBinaryVolume(nz,nx,ny,f"outputs/models/estimatedModel_iteration_10.bin")

model[-1,:,:] = model[-2,:,:]
model[:,-1,:] = model[:,-2,:]
model[:,:,-1] = model[:,:,-2]

modelPerc = gaussian_filter(model, 0)

# modelPerc[:4,:,:] = 0

sx,sy,sz = np.loadtxt(f"outputs/geometry/shotsPosition.txt", delimiter = ",",unpack=True)
rx,ry,rz = np.loadtxt(f"outputs/geometry/nodesPosition.txt", delimiter = ",",unpack=True)

xzPlane = int(ny / 2)
yzPlane = int(nx / 2)
xyPlane = int(nz / 2)

vmin = 1500
vmax = 3500 

xloc = np.linspace(0,nx-1,7,dtype=int)
xlab = np.array(xloc * dh,dtype=int)

yloc = np.linspace(0,ny-1,7,dtype=int)
ylab = np.array(yloc * dh,dtype=int)

zloc = np.linspace(0,nz-1,7,dtype=int)
zlab = np.array(zloc * dh,dtype=int)

plt.figure(1,figsize=(15, 5))

plt.suptitle("Illumination model 200 m", fontsize = 20)

G = gridspec.GridSpec(2, 5)

axes_1 = plt.subplot(G[:1,:3])
plt.imshow(modelPerc[:,xzPlane,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.title(f"XZ plane", fontsize = 18)
plt.xlabel("X axis [m]", fontsize = 15)
plt.ylabel("Z axis [m]", fontsize = 15)

# plt.scatter(sx/dh,sz/dh, label = "Shots")
# plt.scatter(rx/dh,rz/dh, label = "Nodes") 

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),'--r')
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,'--r')

# plt.legend(loc = "lower left", fontsize = 10)

axes_2 = plt.subplot(G[1:,:3])
plt.imshow(modelPerc[:,:,yzPlane],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)
plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.title(f"YZ plane", fontsize = 18)
plt.xlabel("Y axis [m]", fontsize = 15)
plt.ylabel("Z axis [m]", fontsize = 15)

# plt.scatter(sy/dh,sz/dh, label = "Shots")
# plt.scatter(ry/dh,rz/dh, label = "Nodes")

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),'--r')
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,'--r')

# plt.legend(loc = "lower left", fontsize = 10)

axes_3 = plt.subplot(G[:,3:])
plt.imshow(modelPerc[xyPlane,:,:],aspect="auto",cmap="Greys",vmin=vmin,vmax=vmax)

plt.xticks(xloc,xlab)
plt.yticks(yloc,ylab)
plt.title(f"XY plane", fontsize = 18)
plt.xlabel("X axis [m]", fontsize = 15)
plt.ylabel("Y axis [m]", fontsize = 15)

cbar = plt.colorbar()
cbar.set_label("P wave velocity [m/s]",fontsize=15)

# plt.scatter(sx/dh,sy/dh, s = 5, label = "Shots")
# plt.scatter(rx/dh,ry/dh, s = 10, label = "Nodes")

plt.plot(np.arange(nx),np.ones(nx)*xzPlane,'--r', markersize = 10)
plt.plot(np.ones(ny)*yzPlane,np.arange(ny),'--r', markersize = 10)

# plt.legend(loc = "lower left", fontsize = 10)

plt.gca().invert_yaxis()

plt.tight_layout()

plt.show()

trueModel = readBinaryVolume(nz,nx,ny,f"inputs/models/trueModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
initModel = readBinaryVolume(nz,nx,ny,f"inputs/models/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
brrnModel = readBinaryVolume(nz,nx,ny,f"outputs/models/estimatedModel_iteration_10.bin")
# tkv0Model = readBinaryVolume(nz,nx,ny,f"outputs/models/model_5it_tk0.bin")
# tkv1Model = readBinaryVolume(nz,nx,ny,f"outputs/models/model_5it_tk1.bin")
# tkv2Model = readBinaryVolume(nz,nx,ny,f"outputs/models/model_5it_tk2.bin")

trueModel = trueModel[:,int(nx/2),int(ny/2)]
initModel = initModel[:,int(nx/2),int(ny/2)]
brrnModel = brrnModel[:,int(nx/2),int(ny/2)]
# tkv0Model = tkv0Model[:,int(nx/2),int(ny/2)]
# tkv1Model = tkv1Model[:,int(nx/2),int(ny/2)]
# tkv2Model = tkv2Model[:,int(nx/2),int(ny/2)]

n = 2

brrnModel = 1.0 / gaussian_filter(1.0 / brrnModel, n)
# tkv0Model = 1.0 / gaussian_filter(1.0 / tkv0Model, n)
# tkv1Model = 1.0 / gaussian_filter(1.0 / tkv1Model, n)
# tkv2Model = 1.0 / gaussian_filter(1.0 / tkv2Model, n)

depth = np.arange(nz) * dh

ylab = np.linspace(0,depth[-1],11, dtype = int)

plt.figure(2, figsize=(15, 7))

G = gridspec.GridSpec(1, 4)

ax1 = plt.subplot(G[:,:1])

plt.plot(trueModel, depth, label = "True model")
plt.plot(initModel, depth, label = "Initial model")
plt.plot(brrnModel, depth, label = "Recovered model")

plt.ylim([0,(nz-1)*dh])

plt.title("Berriman", fontsize = 18)
plt.ylabel("Depth [m]", fontsize = 18)
plt.xlabel("P wave velocity [m/s]", fontsize = 15)

plt.legend(loc = "upper right")

plt.yticks(ylab,ylab)

plt.gca().invert_yaxis()

# ax2 = plt.subplot(G[:,1:2])

# plt.plot(trueModel, depth, label = "True model")
# plt.plot(initModel, depth, label = "Initial model")
# plt.plot(tkv0Model, depth, label = "Recovered model")

# plt.title("Zero order Tikhonov", fontsize = 18)
# plt.ylabel("Depth [m]", fontsize=18)
# plt.xlabel("P wave velocity [m/s]", fontsize = 15)

# plt.ylim([0,(nz-1)*dh])

# plt.legend(loc = "upper right")

# plt.gca().invert_yaxis()

# ax3 = plt.subplot(G[:,2:3])
# plt.plot(trueModel, depth, label = "True model")
# plt.plot(initModel, depth, label = "Initial model")
# plt.plot(tkv1Model, depth, label = "Recovered model")

# plt.ylim([0,(nz-1)*dh])

# plt.title("First order Tikhonov", fontsize = 18)
# plt.ylabel("Depth [m]", fontsize=18)
# plt.xlabel("P wave velocity [m/s]", fontsize = 15)

# plt.legend(loc = "upper right")

# plt.gca().invert_yaxis()

# ax4 = plt.subplot(G[:,3:])

# plt.plot(trueModel, depth, label = "True model")
# plt.plot(initModel, depth, label = "Initial model")
# plt.plot(tkv2Model, depth, label = "Recovered model")

# plt.ylim([0,(nz-1)*dh])

# plt.title("Second order Tikhonov", fontsize = 18)
# plt.ylabel("Depth [m]", fontsize=18)
# plt.xlabel("P wave velocity [m/s]", fontsize = 15)

# plt.legend(loc = "upper right")

# plt.gca().invert_yaxis()

plt.tight_layout()
plt.show()

# resBer = np.loadtxt("outputs/convergency/residuo_berriman.txt")
# resTk0 = np.loadtxt("outputs/convergency/residuo_tk0.txt")
# resTk1 = np.loadtxt("outputs/convergency/residuo_tk1.txt")
# resTk2 = np.loadtxt("outputs/convergency/residuo_tk2.txt")

# iteration = np.arange(6)

# plt.figure(1, figsize=(5,8))

# plt.plot(resBer, iteration, label = "Berriman")
# plt.plot(resTk0, iteration, label = "0th order Tikhonov")
# plt.plot(resTk1, iteration, label = "1st order Tikhonov")
# plt.plot(resTk2, iteration, label = "2nd order Tikhonov")

# plt.xlim([0, 5])
# plt.xlim([60,140])

# plt.title("Convergence curve", fontsize = 20)
# plt.ylabel("Iteration number", fontsize=18)
# plt.xlabel(r"Residuous norm $||d_{obs} - d_{cal}||_2^2$", fontsize = 18)

# plt.legend(loc = "lower right", fontsize = 15)

# plt.tight_layout()
# plt.gca().invert_yaxis()
# plt.show()


