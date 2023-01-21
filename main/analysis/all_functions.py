import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from math import erf
from scipy.linalg import inv, solve, norm

from matplotlib import patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

def readBinaryArray(n,filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

def readBinaryMatrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)    
    return np.reshape(data, [n1,n2], order='F')

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

def analyticalRefraction(v,z,x):

    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

def fullBoxPlot(model, dh, shots, nodes, slices, xzRay, eikonal, colorbarFix = 1.5):
    
    xyModel = model[slices[2],:,:].T
    zxModel = model[:,:,slices[0]]
    zyModel = model[:,slices[1],:].T

    xyEikonal = eikonal[slices[2],:,:].T
    zxEikonal = eikonal[:,:,slices[0]]
    zyEikonal = eikonal[:,slices[1],:].T

    ticks = np.array([3,9,3], dtype = int)

    #------------------------------------------------
    axis = np.array(np.shape(model))
    [nz,nx,ny] = axis
    [z, x, y] = axis * 0.6 / np.max(axis)

    vmin = np.min(model)
    vmax = np.max(model)

    px = 1/plt.rcParams['figure.dpi']  
    fig = plt.figure(1, figsize=(500*px, 200*px))

    xloc = np.linspace(0,nx,ticks[1], dtype = int)
    yloc = np.linspace(0,ny,ticks[2], dtype = int)
    zloc = np.linspace(0,nz,ticks[0], dtype = int)

    xlab = np.around(xloc * dh[1] / 1000, decimals = 1)
    ylab = np.around(yloc * dh[2] / 1000, decimals = 1)
    zlab = np.around(zloc * dh[0] / 1000, decimals = 1)

    #------------------------------------------------
    ax1 = fig.add_axes([0.7 - x, 0.95 - y, 0 + x, 0 + y])
    ax1.contour(xyEikonal, levels = 10, cmap = "seismic")
    ax1.imshow(xyModel, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax1.tick_params(direction = 'in', axis='x') 
    ax1.tick_params(direction = 'out', axis='y') 
    ax1.set_xticklabels([])

    ax1.plot(np.arange(axis[1]), np.ones(axis[1])*slices[2], "--g", alpha = 0.3)
    ax1.plot(np.ones(axis[2])*slices[1], np.arange(axis[2]), "--m", alpha = 0.3)

    ax1.scatter(shots[0]/dh[0], shots[1]/dh[1])
    ax1.scatter(nodes[:,0]/dh[0], nodes[:,1]/dh[1])

    ax1.set_xticks(xloc)
    ax1.set_yticks(yloc[1:])
    ax1.set_yticklabels(ylab[1:])

    ax1.grid(color='w', linestyle='--', linewidth=0.3)

    ax1.set_ylabel("Y [km]")
    ax1.invert_yaxis()

    #------------------------------------------------
    ax2 = fig.add_axes([0.7, 0.95 - y, 0 + z, 0 + y])
    ax2.contour(zyEikonal, levels = 10, cmap = "seismic")
    ax2.imshow(zyModel, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)
    ax2.tick_params(direction = 'out', axis='x') 
    ax2.tick_params(direction = 'in', axis='y') 
    ax2.set_yticklabels([])

    ax2.set_xlabel("Z [km]")

    ax2.plot(np.arange(axis[0]), np.ones(axis[0])*slices[2], "--g", alpha = 0.3)
    ax2.plot(np.ones(axis[2])*slices[0], np.arange(axis[2]), "--r", alpha = 0.3)

    ax2.scatter(shots[2]/dh[2], shots[1]/dh[1])
    ax2.scatter(nodes[:,2]/dh[2], nodes[:,1]/dh[1])

    ax2.set_yticks(yloc)

    ax2.set_xticks(zloc[1:])
    ax2.set_xticklabels(zlab[1:])

    ax2.grid(color = 'w', linestyle = '--', linewidth = 0.3)

    #------------------------------------------------
    ax3 = fig.add_axes([0.7 - x, 0.95 - y - z, 0 + x, 0 + z])
    ax3.contour(zxEikonal, levels = 10, cmap = "seismic")
    ax3.imshow(zxModel, aspect = 'auto', cmap = "Greys", vmin = vmin, vmax = vmax)

    ax3.plot(np.arange(axis[1]), np.ones(axis[1])*slices[0], "--r", alpha = 0.3)
    ax3.plot(np.ones(axis[0])*slices[1], np.arange(axis[0]), "--m", alpha = 0.3)

    ax3.scatter(shots[0]/dh[0], shots[2]/dh[2])
    ax3.scatter(nodes[:,0]/dh[0], nodes[:,2]/dh[2])

    ax3.scatter(xzRay[0,:]/dh[0], xzRay[1,:]/dh[2], s = 0.1, c = "green")

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
    cbar.set_label("Velocidade P [km/s]")
    
    return None

def buildModel3D(nx,ny,nz,property,z):

    model = np.zeros((nz,nx,ny), dtype=float)

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1, len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model

def createGaussianSurface(nx,ny,dx,dy,A,xc,yc,sigx,sigy):

    x, y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)

    surface = A*np.exp(-((x - xc)**2/(2*sigx**2) + (y - yc)**2/(2*sigy**2)))

    return surface

def multiBoxPlot(models:np.ndarray, shots:np.ndarray, nodes:np.ndarray, dh:np.ndarray, slices:np.ndarray, subplots:np.ndarray) -> None:

    if np.sum(subplots) == 2:
        modelShape = np.array(np.shape(models))
        maxModelDistance = np.max(np.shape(models))
        minModelDistance = np.min(np.shape(models))
        
        vmin = np.min(models)
        vmax = np.max(models)

    else:
        modelShape = np.array(np.shape(models[0]))
        maxModelDistance = np.max(np.shape(models[0]))
        minModelDistance = np.min(np.shape(models[0]))
        
        vmin = np.min(models[0])
        vmax = np.max(models[0])

    nz, nx, ny = modelShape
    [z, x, y] = 2.0 * (minModelDistance / maxModelDistance) * modelShape / maxModelDistance

    px = 1/plt.rcParams['figure.dpi']  
    ticks = np.array([3,7,7], dtype = int)

    fig = plt.figure(1, figsize=(700*px*subplots[1], 600*px*subplots[0]))

    xloc = np.linspace(0,nx-1,ticks[1], dtype = int)
    yloc = np.linspace(0,ny-1,ticks[2], dtype = int)
    zloc = np.linspace(0,nz-1,ticks[0], dtype = int)

    m2km = 1e-3

    xlab = np.around(xloc * dh[0] * m2km, decimals = 1)
    ylab = np.around(yloc * dh[1] * m2km, decimals = 1)
    zlab = np.around(zloc * dh[2] * m2km, decimals = 1)

    axes = np.array([[0.75 - x, 0.98 - y      , x, y], 
                     [    0.75, 0.98 - y      , z, y],
                     [0.75 - x, 0.98 - y - z  , x, z],
                     [0.75 - x, 0.98 - y - 2.5*z, x, z]])

    xTickDirection = ['out', 'out', 'out']
    yTickDirection = ['out', 'in', 'out']

    xTickLock = [xloc, zloc[1:], xloc]
    yTickLock = [yloc, yloc, zloc[1:]]

    xTickLabel = [[], zlab[1:], xlab]
    yTickLabel = [ylab, [], zlab[1:]]

    xLabel = ["X [km]", "Z [km]", "X [km]"]
    yLabel = ["Y [km]", "      ", "Z [km]"]

    yInvert = [ True, True, False]

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

            if np.sum(subplots) == 2:
                ims = [models[slices[0],:,:].T, models[:,slices[2],:].T, models[:,:,slices[1]]]
            else:
                ims = [models[ind, slices[0],:,:].T, models[ind,:,slices[2],:].T, models[ind,:,:,slices[1]]]

            xshot = [shots[:,0]/dh[0], shots[:,2]/dh[2], shots[:,0]/dh[0]]
            yshot = [shots[:,1]/dh[1], shots[:,1]/dh[1], shots[:,2]/dh[2]]

            xnode = [nodes[:,0]/dh[0], nodes[:,2]/dh[2], nodes[:,0]/dh[0]]
            ynode = [nodes[:,1]/dh[1], nodes[:,1]/dh[1], nodes[:,2]/dh[2]]

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
                    
                    ax.scatter(xshot[k], yshot[k], s = 0.01)
                    ax.scatter(xnode[k], ynode[k], s = 5.0)

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

def irls(A, b, tolr, tolx, p, maxiter):
    ''' 
    Solve for the 1-norm solution

        Input Parameters:
            A       - Matrix of the system of equations.
            b       - Right hand side of the system of equations.
            tolr    - Tolerance below which residuals are ignored.
            tolx    - Stopping tolerance.  Stop when (norm(newx-x)/(1+norm(x)) < tolx)
            p       - Specifies which p-norm to use (most often, p=1.)
            maxiter - Limit on number of iterations of IRLS

        Output Parameters:
            x - Approximate L_p solution.
    '''

    # Find the size of the matrix A.
    [m,n] = np.shape(A)

    # Start the first iteration with R=I, and x=A\b (the least squares solution)
    R = np.eye(m)
    x = solve(A.T @ A, A.T @ b)

    #  Now loop up to maxiter iterations
    iter = 1
    while iter <= maxiter:
        
        iter += 1

        # compute the current residual 
        r = A @ x - b

        # for each row adjust the weighting factor r based on the residual
        for i in range(m):
            if np.abs(r[i] < tolr):
                r[i] = np.abs(tolr)**(p - 2)
            else:
                r[i] = np.abs(r[i])**(p - 2)

        # insert the weighting factors into R
        R = np.diag(r)

        # find the solution to the weighted problem
        newx = solve(A.T @ R @ A, A.T @ R @ b)

        # check for convergence
        if norm(newx - x) / (1 + norm(x)) < tolx:
            x = newx
            return x
        else:
            x = newx
