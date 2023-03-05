import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

from functions import *
#------------------------------------------------------------------------------------

dx = 10
dy = 10
dz = 10

nx = 701
ny = 501
nz = 101

#----------------------------------------------------------------------------

waterBottom = np.zeros((ny,nx))

top = 200      
base = 300

A = np.array([100, -50, -50])
xc = np.array([100, 5000, 5000])
yc = np.array([2500, 100, 4900])
sigx = np.array([2000, 2000, 2000])
sigy = np.array([100000, 1000, 1000])

for k in range(len(A)):
    waterBottom += createGaussianSurface(nx,ny,dx,dy,A[k],xc[k],yc[k],sigx[k],sigy[k])

waterBottom = top//dz + np.array((np.max(waterBottom)-waterBottom)//dz, dtype=int)

i,j = np.where(-waterBottom*dz < -base)
waterBottom[i,j] = base/dz

# Defining geometry with topography

node_xline = 16
node_yline = 11
node_all = node_xline * node_yline

node_x = np.linspace(500, 6500, node_xline)
node_y = np.linspace(500, 4500, node_yline)

node_z = np.zeros(node_all)

for i in range(node_yline):
    for j in range(node_xline):
        node_z[i + node_yline*j] = dz*waterBottom[int(node_y[i]/dy),int(node_x[j]/dx)]

topo_file = open("../../../inputs/geometry/nodes_topography.txt", 'w')

for k in range(len(node_z)):
    topo_file.write(f"{node_z[k]}\n")
    
topo_file.close()

node_x, node_y = np.meshgrid(node_x, node_y)

node_x = np.reshape(node_x, [node_all], order = "F")
node_y = np.reshape(node_y, [node_all], order = "F")

shot_xline = 140
shot_yline = 100

shot_all = shot_xline * shot_yline

shot_x = np.linspace(25, 6975, shot_xline)
shot_y = np.linspace(25, 4975, shot_yline)

shot_z = np.ones(shot_all) * 25.0

shot_x, shot_y = np.meshgrid(shot_x, shot_y)

shot_x = np.reshape(shot_x, [shot_all], order = "F")
shot_y = np.reshape(shot_y, [shot_all], order = "F")

plt.figure(1, figsize = (10, 6))
plt.scatter(shot_x, shot_y, s = 1.0, label = "Posição dos tiros")
plt.scatter(node_x, node_y, s = 20.0, label = "Posição dos receptores")

plt.title("Geometria de aquisição", fontsize = 18)
plt.xlabel("X [m]", fontsize = 15)
plt.ylabel("Y [m]", fontsize = 15)

plt.legend(loc = "lower left")
plt.xlim(0,7000)
plt.ylim(0,5000)
plt.tight_layout()
plt.savefig("complete_geometry.png", dpi = 200)
plt.show()

# water bottom surface plot
fig, ax = plt.subplots(1,1, figsize = (10,8))

cmap = mpl.colormaps["coolwarm"]
smin = np.min(-dz*waterBottom)
smax = np.max(-dz*waterBottom)

ax.contour(-dz*waterBottom, levels = 10)
ax.imshow(-dz*waterBottom, aspect = "auto", cmap=cmap)
ax.scatter(node_x/dx, node_y/dy, c = "black", s = 8.0, label = "Posição dos receptores")
# ax.scatter(shot_x/dx, shot_y/dy, c = "green", s = 0.1, label = "Shots position")
# ax.scatter(xc/dx, yc/dy, c = "black", s = 30, label = "Centros das funções gaussianas")

ax.set_title("Topografia do fundo marinho", fontsize = 18)
ax.set_xlabel("X [m]", fontsize = 15)
ax.set_ylabel("Y [m]", fontsize = 15)
ax.invert_yaxis()

ax.set_xticks(np.linspace(0,nx-1,9, dtype = int))
ax.set_xticklabels(np.array(np.linspace(0,nx-1,9)*dx, dtype = int))

ax.set_yticks(np.linspace(0,ny-1,9, dtype = int))
ax.set_yticklabels(np.array(np.linspace(0,ny-1,9)*dy, dtype = int))

norm = mpl.colors.Normalize(smin,smax)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.6)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(smin, smax, 5), orientation = "horizontal")
cbar.ax.set_xticklabels(np.around(np.linspace(smin, smax, 5), decimals = 1))
cbar.set_label("Profundidade [m]", fontsize = 15)

ax.legend(loc = "lower left", fontsize = 12)

plt.tight_layout()
plt.savefig("water_bottom_surface_with_nodes.png", dpi = 200)
plt.show()
#----------------------------------------------------------------------------

top = 350
base = 950

model = np.zeros((nz,nx,ny))

surface = np.zeros((ny, nx))

A = np.array([550, 600, -300, 450, 400])

xc = np.array([2500, 2500, 3500, 4500, 4500])
yc = np.array([2000, 3000, 2500, 2000, 3000])

sigx = np.array([800, 800, 500, 800, 800])
sigy = np.array([500, 500, 300, 600, 600])

for k in range(len(A)):
    surface += createGaussianSurface(nx,ny,dx,dy,A[k],xc[k],yc[k],sigx[k],sigy[k])

surface = top//dz + np.array((np.max(surface)-surface)//dz, dtype=int)

i,j = np.where(-surface*dz < -base)
surface[i,j] = base/dz

# Surface plot
fig, ax = plt.subplots(1,1, figsize = (10,8))

cmap = mpl.colormaps["coolwarm"]
smin = np.min(-dz*surface)
smax = np.max(-dz*surface)

ax.contour(-dz*surface, levels = 10)
ax.imshow(-dz*surface, aspect = "auto", cmap=cmap)
# ax.scatter(node_x/dx, node_y/dy, c = "black", s = 5.0, label = "Posição dos receptores")
# ax.scatter(shot_x/dx, shot_y/dy, c = "green", s = 0.1, label = "Posição dos tiros")
ax.scatter(xc/dx, yc/dy, c = "black", s = 30, label = "Centros das funções gaussianas")

ax.set_xticks(np.linspace(0,nx-1,9, dtype = int))
ax.set_xticklabels(np.array(np.linspace(0,nx-1,9)*dx, dtype = int))

ax.set_yticks(np.linspace(0,ny-1,9, dtype = int))
ax.set_yticklabels(np.array(np.linspace(0,ny-1,9)*dy, dtype = int))

ax.invert_yaxis()

ax.set_title("Topografia da superfície alvo", fontsize = 18)
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
plt.savefig("target_surface_gaussians.png", dpi = 200)
plt.show()

# sys.exit()
#----------------------------------------------------------------------------
v0 = 1650.0 # m/s
dv = 10.0   # m/s

wb = 250

v = np.array([1500, 2000, 3500])
z = np.array([wb, base, nz*dz])

initModel = buildModel3D(nx, ny, nz, v, z//dz)

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

dh = np.array([dx, dy, dz])

shots = np.zeros((shot_all, 3))
nodes = np.zeros((node_all, 3))

shots[:, 0] = shot_x
shots[:, 1] = shot_y
shots[:, 2] = shot_z

nodes[:, 0] = node_x
nodes[:, 1] = node_y
nodes[:, 2] = node_z

i,j,k = np.where(model <= 1500)

vp = model.copy()
vs = model / 1.7
rho = 310 * model ** 0.25

vs[i,j,k] = 0.0
rho[i,j,k] = 1000

plt.figure(5, figsize = (18, 5))
plt.subplot(131)
plt.imshow(vp[:,:,250], aspect = "auto")

plt.subplot(132)
plt.imshow(vs[:,:,250], aspect = "auto")

plt.subplot(133)
plt.imshow(rho[:,:,250], aspect = "auto")

plt.tight_layout()
plt.show()

vp.flatten("F").astype("float32", order = "F").tofile(f"../../../inputs/models/vp_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
vs.flatten("F").astype("float32", order = "F").tofile(f"../../../inputs/models/vs_{nz}x{nx}x{ny}_{dx:.0f}m.bin")
rho.flatten("F").astype("float32", order = "F").tofile(f"../../../inputs/models/rho_{nz}x{nx}x{ny}_{dx:.0f}m.bin")

initModel.flatten("F").astype("float32", order = "F").tofile(f"../../../inputs/models/initModel_{nz}x{nx}x{ny}_{dx:.0f}m.bin")


#slices = np.array([int(600/dz), int(2500/dy), int(3300/dx)], dtype = int) # xy, zx, zy
#subplots = np.array([1, 1], dtype = int)

#check_geometry(model, shots, nodes, dh, slices, subplots)
#plt.savefig("true_model.png")
#plt.show()

#check_geometry(initModel, shots, nodes, dh, slices, subplots)
#plt.savefig("init_model.png")
#plt.show()

# Low frequency initial model
#model.flatten("F").astype("float32", order = "F").tofile(f"../../../inputs/models/trueModel_{nz}x{nx}x{ny}_{dx:.1f}m.bin")

