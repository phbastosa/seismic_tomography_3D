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
    
def boxPlot(models:np.ndarray, dh:np.ndarray, slices:np.ndarray, subplots:np.ndarray) -> None:

    modelShape = np.array(np.shape(models[0]))
    maxModelDistance = np.max(np.shape(models[0]))
    minModelDistance = np.min(np.shape(models[0]))
    
    [z, x, y] = 3.0 * (minModelDistance / maxModelDistance) * modelShape / maxModelDistance

    vmin = np.min(models[0])
    vmax = np.max(models[0])

    px = 1/plt.rcParams['figure.dpi']  
    ticks = np.array([3,7,7], dtype = int)

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
                     [0.75 - x, 0.98 - y - 2*z, x, z]])

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
                    cbar.set_label("P wave velocity [km/s]")
                
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
    
    plt.savefig("testFigure.png", dpi = 200)
    plt.show(block = False)

    return None


################################################################################################

nx = 101
ny = 101
nz = 21

dx = 50.0
dy = 50.0
dz = 50.0

nIt = 7

# Importing geometry acquisition

shots = np.loadtxt(f"outputs/shotsPosition.txt", delimiter = ",",unpack=True)
nodes = np.loadtxt(f"outputs/nodesPosition.txt", delimiter = ",",unpack=True)

# Importing models
models = np.zeros((nIt+2, nz, nx, ny))

models[0,:,:,:] = readBinaryVolume(nz,nx,ny,f"outputs/initModel_{nz}x{nx}x{ny}_{dx:.0f}m.bin")

for i in range(1,nIt+1):
    models[i,:,:,:] = readBinaryVolume(nz,nx,ny,f"outputs/estimatedModel_iteration_{i}.bin")

models[-1,:,:,:] = readBinaryVolume(nz,nx,ny,f"outputs/trueModel_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
 
# Defining plane positions

dh = np.array([dz, dx, dy])
slices = np.array([int(nz / 2), int(ny / 2), int(nx / 2)], dtype = int) # [xy, zx, zy]
subplots = np.array([3, 3], dtype = int)

# Plotting models

boxPlot(models, dh, slices, subplots)

# Plotting convergency 

convergency = np.loadtxt("outputs/convergency.txt")

plt.figure(2, figsize = (10,5))
plt.plot(convergency)
plt.title("Convergency using Tikhonov first order with 1000 of regularization parameter")
plt.xlabel("Iterations number")
plt.ylabel("L2 residuous norm = sum((dobs - dcal)**2)")
plt.savefig("convergency.png", dpi = 200)
plt.show(block = False)