import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

from all_functions import buildModel3D, createGaussianSurface, multiBoxPlot
#------------------------------------------------------------------------------------

dx = 12.5
dy = 12.5
dz = 12.5

nx = 401
ny = 401
nz = 101

#----------------------------------------------------------------------------

waterBottom = np.zeros((ny,nx))

top = 200      
base = 300

A = np.array([10, 50, 100])
xc = np.array([1000, 2500, 4500])
yc = np.array([1000, 2500, 4500])
sigx = np.array([1000, 1000, 1000])
sigy = np.array([1000, 1000, 6000])

for k in range(len(A)):
    waterBottom += createGaussianSurface(nx,ny,dx,dy,A[k],xc[k],yc[k],sigx[k],sigy[k])

waterBottom = top//dz + np.array((np.max(waterBottom)-waterBottom)//dz, dtype=int)

i,j = np.where(-waterBottom*dz < -base)
waterBottom[i,j] = base/dz

# Defining geometry with topography

node_xline = 11
node_yline = 11
node_all = node_xline * node_yline

node_x = np.linspace(500, 4500, node_xline)
node_y = np.linspace(500, 4500, node_yline)
node_z = np.zeros(node_all)

for i in range(node_yline):
    for j in range(node_xline):
        node_z[i + node_yline*j] = dz*waterBottom[int(node_y[i]/dy),int(node_x[j]/dx)]

node_z.astype("float32", order = "F").tofile(f"../../inputs/geometry/nodesTopography.bin")

node_x, node_y = np.meshgrid(node_x, node_y)
node_x = np.reshape(node_x, [node_all])
node_y = np.reshape(node_y, [node_all])

shot_xline = 100
shot_yline = 100
shot_all = shot_xline * shot_yline

shot_x = np.linspace(25, 4975, shot_xline)
shot_y = np.linspace(25, 4975, shot_yline)
shot_z = np.ones(shot_all) * 25.0

shot_x, shot_y = np.meshgrid(shot_x, shot_y)
shot_x = np.reshape(shot_x, [shot_all])
shot_y = np.reshape(shot_y, [shot_all])

# water bottom surface plot
fig, ax = plt.subplots(1,1, figsize = (8,8))

cmap = mpl.colormaps["coolwarm"]
smin = np.min(-dz*waterBottom)
smax = np.max(-dz*waterBottom)

ax.contour(-dz*waterBottom, levels = 10)
ax.imshow(-dz*waterBottom, aspect = "auto", cmap=cmap)
ax.scatter(node_x/dx, node_y/dy, c = "black", s = 5.0, label = "Posição dos nodes")
ax.scatter(shot_x/dx, shot_y/dy, c = "green", s = 0.01, label = "Posição dos tiros")

ax.set_title("Topografia do fundo marinho", fontsize = 18)
ax.set_xlabel("X [m]", fontsize = 15)
ax.set_ylabel("Y [m]", fontsize = 15)
ax.invert_yaxis()

ax.set_xticks(np.linspace(0,nx-1,9, dtype = int))
ax.set_xticklabels(np.around(np.linspace(0,nx-1,9)*dx, decimals = 1))

ax.set_yticks(np.linspace(0,ny-1,9, dtype = int))
ax.set_yticklabels(np.around(np.linspace(0,ny-1,9)*dy, decimals = 1))

norm = mpl.colors.Normalize(smin,smax)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.6)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(smin, smax, 5), orientation = "horizontal")
cbar.ax.set_xticklabels(np.around(np.linspace(smin, smax, 5), decimals = 1))
cbar.set_label("Profundidade [m]", fontsize = 15)

ax.legend(loc = "upper right", fontsize = 12)

plt.tight_layout()
plt.savefig("../../figures/2_waterBottomSurface.png", dpi = 200)
plt.show()

#----------------------------------------------------------------------------

top = 500
base = 1000

model = np.zeros((nz,nx,ny))

surface = np.zeros((ny, nx))

A = np.array([500, 400, -300, 400, 500])
xc = np.array([2000, 3000, 2500, 2000, 3000])
yc = np.array([2000, 2000, 2500, 3000, 3000])
sigx = np.array([500, 400, 400, 500, 400])
sigy = np.array([400, 500, 400, 400, 500])

for k in range(len(A)):
    surface += createGaussianSurface(nx,ny,dx,dy,A[k],xc[k],yc[k],sigx[k],sigy[k])

surface = top//dz + np.array((np.max(surface)-surface)//dz, dtype=int)

i,j = np.where(-surface*dz < -base)
surface[i,j] = base/dz

# Surface plot
fig, ax = plt.subplots(1,1, figsize = (8,8))

cmap = mpl.colormaps["coolwarm"]
smin = np.min(-dz*surface)
smax = np.max(-dz*surface)

ax.contour(-dz*surface, levels = 10)
ax.imshow(-dz*surface, aspect = "auto", cmap=cmap)
ax.scatter(node_x/dx, node_y/dy, c = "black", s = 5.0, label = "Posição dos nodes")
ax.scatter(shot_x/dx, shot_y/dy, c = "green", s = 0.01, label = "Posição dos tiros")

ax.set_xticks(np.linspace(0,nx-1,9, dtype = int))
ax.set_xticklabels(np.around(np.linspace(0,nx-1,9)*dx, decimals = 1))

ax.set_yticks(np.linspace(0,ny-1,9, dtype = int))
ax.set_yticklabels(np.around(np.linspace(0,ny-1,9)*dy, decimals = 1))

ax.invert_yaxis()

ax.set_title("Topografia da anomalia", fontsize = 18)
ax.set_xlabel("X [m]", fontsize = 15)
ax.set_ylabel("Y [m]", fontsize = 15)

norm = mpl.colors.Normalize(smin,smax)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.6)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(smin, smax, 5), orientation = "horizontal")
cbar.ax.set_xticklabels(np.around(np.linspace(smin, smax, 5), decimals = 1))
cbar.set_label("Profundidade [m]", fontsize = 15)
ax.legend(loc = "upper right", fontsize = 12)

plt.tight_layout()
plt.savefig("../../figures/3_anomalySurface.png", dpi = 200)
plt.show()

#----------------------------------------------------------------------------
v0 = 1550.0 # m/s
dv = 25.0   # m/s

wb = 250

v = np.array([1500,2000,3500])
z = np.array([wb,base,nz*dz])

initModel = buildModel3D(nx,ny,nz,v,z//dz)

for j in range(nx):
    for k in range(ny):
        nGradient = int(surface[k,j]) - int(waterBottom[k,j])
        gradient = slice(int(waterBottom[k,j]), int(surface[k,j]))

        model[int(surface[k,j]):,j,k] = v[-1]
        model[:int(waterBottom[k,j]), j, k] = v[0]
        model[gradient, j, k] = v0 + dv*np.arange(nGradient)

        initModel[int(base/dz):,j,k] = v[-1]
        initModel[:int(wb/dz), j, k] = v[0]
        initModel[int(wb/dz):int(base/dz), j, k] = v0 + dv*np.arange((int(base/dz)-int(wb/dz)))
        
#----------------------------------------------------------------------------

models = np.zeros((2,nz,nx,ny))

models[0,:,:,:] = model
models[1,:,:,:] = initModel

dh = np.array([dx, dy, dz])

shots = np.zeros((shot_all, 3))
nodes = np.zeros((node_all, 3))

shots[:, 0] = shot_x
shots[:, 1] = shot_y
shots[:, 2] = shot_z

nodes[:, 0] = node_x
nodes[:, 1] = node_y
nodes[:, 2] = node_z

slices = np.array([int(nz / 2), int(ny / 2), int(nx / 2)], dtype = int) # [xy, zx, zy]
subplots = np.array([1, 2], dtype = int)

multiBoxPlot(models, shots, nodes, dh, slices, subplots)
plt.savefig("../../figures/4_benchmark_models.png")
plt.show()

# Low frequency initial model
initModel.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/initModel_{nz}x{nx}x{ny}_{dx:.1f}m.bin")
model.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/trueModel_{nz}x{nx}x{ny}_{dx:.1f}m.bin")

