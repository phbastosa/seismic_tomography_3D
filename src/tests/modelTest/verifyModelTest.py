import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryVolume(dim1,dim2,dim3,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
    return np.reshape(data, [dim1,dim2,dim3], order='F')

def boxPlot(model, xySlice, zySlice, zxSlice):
    
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
    ax4 = fig.add_axes([0.7 - x, 0.95 - y - 1.70*z, 0 + x, 0 + z])

    ax4.axis("off")

    cmap = mpl.colormaps["Greys"]
    norm = mpl.colors.Normalize(vmin*1e-3, vmax*1e-3)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("bottom", size="10%", pad=0)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(vmin*1e-3, vmax*1e-3, 5), orientation = "horizontal")
    cbar.ax.set_xticklabels(np.around(np.linspace(vmin*1e-3, vmax*1e-3, 5), decimals = 1))
    cbar.set_label("P wave velocity [km/s]")
    
    return None

nb = 50

nx = 801
ny = 801
nz = 187

dx = 25.0
dy = 25.0
dz = 25.0

reducedModel = readBinaryVolume(nz, nx, ny,"outputs/innerModel.bin")
expandedModel = readBinaryVolume(nz + 2*nb, nx + 2*nb, ny + 2*nb,"outputs/expandedModel.bin")

xySlice = int(nz/2)
zySlice = int(nx/2)
zxSlice = int(ny/2)


boxPlot(expandedModel, xySlice, zySlice, zxSlice)
plt.show()

boxPlot(reducedModel, xySlice, zySlice, zxSlice)
plt.show()
