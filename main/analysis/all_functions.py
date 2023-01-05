import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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

def rickerGenerator(fmax,nsrc,dt):    
    ricker = np.zeros(nsrc)
    fc = fmax/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1 - 2 * np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
        aux2 = np.exp(-np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))
        ricker[i + s] = aux1 * aux2 

    return ricker

def zplane(z,p,k,subplot):

    ax = plt.subplot(subplot)
    uc = patches.Circle((0,0), radius=1, fill=False, color='black', ls='dashed')
    ax.add_patch(uc)
    
    # Plot the zeros and set marker properties    
    t1 = plt.plot(z.real, z.imag, 'go', ms=10)
    plt.setp( t1, markersize=10.0, markeredgewidth=1.0, markeredgecolor='k', markerfacecolor='g')

    # Plot the poles and set marker properties
    t2 = plt.plot(p.real, p.imag, 'rx', ms=10)
    plt.setp( t2, markersize=12.0, markeredgewidth=3.0, markeredgecolor='r', markerfacecolor='r')

    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_title("Diagrama de polos e zeros")
    
    # set the ticks
    r = 1.5; plt.axis('scaled'); plt.axis([-r, r, -r, r])
    ticks = [-1, -.5, .5, 1]; plt.xticks(ticks); plt.yticks(ticks)

def fourierInterpolation(signal, dt, new_dt):

    ns = len(signal)

    periods = int(dt / new_dt)

    freqs = np.fft.fftfreq(ns, dt)
    fft_signal = np.fft.fft(signal)

    new_fft_signal = np.insert(fft_signal, int(ns / 2), np.zeros(periods * ns)) 

    new_ns = periods * ns + ns

    scale_factor = new_ns / ns
    output_signal = scale_factor * np.real(np.fft.ifft(new_fft_signal))
    
    return output_signal

def analiticalRefraction(v,z,x):

    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

def boxPlot(model, dh, shots, nodes, slices, xzRay, eikonal, colorbarFix = 1.5):
    
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
    ax1.scatter(nodes[0]/dh[0], nodes[1]/dh[1])

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
    ax2.scatter(nodes[2]/dh[2], nodes[1]/dh[1])

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
    ax3.scatter(nodes[0]/dh[0], nodes[2]/dh[2])

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
    cbar.set_label("P wave velocity [km/s]")
    
    return None
