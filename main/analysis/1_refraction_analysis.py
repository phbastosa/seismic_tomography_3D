import numpy as np
import matplotlib.pyplot as plt

from all_functions import boxPlot, readBinaryArray, readBinaryMatrix, readBinaryVolume, analyticalRefraction

nx = 481
ny = 41
nz = 51

dx = 12.5
dy = 12.5
dz = 12.5

v = np.array([2000,3500])
z = np.array([500]) 

refractiveModel = np.zeros((nz,nx,ny))

refractiveModel[:int(z[0]/dz),:,:] = v[0]
refractiveModel[int(z[0]/dz):,:,:] = v[1]

refractiveModel.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/refractiveModel_{nz}x{nx}x{ny}_{dz:.1f}m.bin")

shots = np.loadtxt("../../outputs/geometry/shots.txt", delimiter = ",")
nodes = np.loadtxt("../../outputs/geometry/nodes.txt", delimiter = ",")

dh = np.array([dx,dy,dz])
slices = np.array([int(nz/2), int(nx/2), int(ny/2)], dtype = int)

n = 112272

xRay = readBinaryArray(n, f"../../outputs/rayPosition/xRay_{n}_shot_1.bin")
yRay = readBinaryArray(n, f"../../outputs/rayPosition/yRay_{n}_shot_1.bin")
zRay = readBinaryArray(n, f"../../outputs/rayPosition/zRay_{n}_shot_1.bin")
iRay = readBinaryArray(n, f"../../outputs/rayPosition/iRay_{n}_shot_1.bin")

xzRay = np.zeros((2, n))

xzRay[0,:] = xRay
xzRay[1,:] = zRay

eikonal = readBinaryVolume(nz,nx,ny,f"../../outputs/travelTimes/eikonal_nz{nz}_nx{nx}_ny{ny}_shot_1.bin")

# boxPlot(refractiveModel, dh, shots, nodes, slices, xzRay, eikonal, colorbarFix = 2.0)
# # print it
# plt.show()

#--------------------------------

offsets = nodes[:,0] - shots[0]

td, tr = analyticalRefraction(v,z,offsets)

analyticalArrivals = np.zeros(len(offsets))

for i in range(len(offsets)):
    if td[i] < tr[0,i]:
        analyticalArrivals[i] = td[i]
    else:
        analyticalArrivals[i] = tr[0,i]
        
dt = 1e-3        
nt = 2501

cutTime = 2.0
delayTime = 0.15 

seismic = readBinaryMatrix(nt,len(offsets), f"../../inputs/seismograms/seismogram_{nt}x{len(offsets)}_shot_1.bin")

seismic = seismic[int(delayTime/dt-1):,:] # removing anti causal 
seismic = seismic[:int(cutTime/dt+1),:]   # cutting in 2 seconds

perc = 0.5 * np.std(seismic)

time = np.arange(len(seismic))

nodeFar = 80
nodeNear = 20

ywells = time
xwells = [np.ones(len(time)) * 20, np.ones(len(time)) * 80]
ampWells = [int(analyticalArrivals[nodeNear]/dt), int(analyticalArrivals[nodeFar]/dt)]

#-----------------------------------------------
fig, ax = plt.subplots(1, 3, figsize = (15, 6))

ax[0].imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -perc, vmax = perc)
ax[0].plot(analyticalArrivals/dt, "-r", label = "Tempo analítico")
ax[0].plot(xwells[0],ywells, label = "Offset 2000 m")
ax[0].plot(xwells[1],ywells, label = "Offset 5000 m")
ax[0].set_yticks(np.linspace(0, len(seismic)-1, 11, dtype = int))
ax[0].set_yticklabels(np.around(np.linspace(0, len(seismic)-1, 11)*dt, decimals = 1))
ax[0].set_xticks(np.linspace(0, len(seismic[0])-1, 7, dtype = int))
ax[0].set_xticklabels(np.linspace(1000, 5500, 7, dtype = int))
ax[0].set_title("Sismograma", fontsize = 18)
ax[0].set_xlabel("Offset [m]", fontsize = 15)
ax[0].set_ylabel("Tempo [s]", fontsize = 15)
ax[0].legend(fontsize = 12, loc = "upper right")

ax[1].plot(seismic[:,nodeNear], time*dt)
ax[1].plot(seismic[ampWells[0],nodeNear], ampWells[0]*dt, "ok", label = f"Tempo analítico = {ampWells[0]*dt:.3f} s")
ax[1].set_yticks(np.linspace(0, len(seismic)-1, 11, dtype = int))
ax[1].set_yticklabels(np.around(np.linspace(0, len(seismic)-1, 11)*dt, decimals = 1))
ax[1].set_xticks(np.linspace(-0.03, 0.03, 5))
ax[1].set_xticklabels(np.around(np.linspace(-0.03, 0.03, 5), decimals = 2))
ax[1].set_title("Traço do offset 2000 m", fontsize = 18)
ax[1].set_xlabel("Amplitude [Pa]", fontsize = 15)
ax[1].set_ylabel("Tempo [s]", fontsize = 15)
ax[1].set_xlim([-0.03,0.03])
ax[1].set_ylim([0, 2])
ax[1].invert_yaxis()
ax[1].legend(fontsize = 12, loc = "upper right")

tpickMax = np.where(seismic[:,nodeFar] == np.max(seismic[:,nodeFar]))[0][0]
tpickMin = np.where(seismic[:,nodeFar] == np.min(seismic[:,nodeFar]))[0][0]

ax[2].plot(seismic[:,nodeFar], time*dt)
ax[2].plot(seismic[ampWells[1],nodeFar], ampWells[1]*dt, "ok", label = f"Tempo analítico = {ampWells[1]*dt:.3f} s")
ax[2].plot(seismic[tpickMax,nodeFar], tpickMax*dt, "or", label = f"Margem máxima = {(tpickMax - ampWells[1])*dt:.3f} s")
ax[2].plot(seismic[tpickMin,nodeFar], tpickMin*dt, "og", label = f"Margem mínima = {(tpickMin - ampWells[1])*dt:.3f} s")
ax[2].set_yticks(np.linspace(0, len(seismic)-1, 11, dtype = int))
ax[2].set_yticklabels(np.around(np.linspace(0, len(seismic)-1, 11)*dt, decimals = 1))
ax[2].set_xticks(np.linspace(-0.00015, 0.00015, 5))
ax[2].set_xticklabels(np.around(np.linspace(-0.00015, 0.00015, 5), decimals = 6))
ax[2].set_title("Traço do offset 5000 m", fontsize = 18)
ax[2].set_xlabel("Amplitude [Pa]", fontsize = 15)
ax[2].set_ylabel("Tempo [s]", fontsize = 15)
ax[2].set_xlim([-0.00015,0.00015])
ax[2].set_ylim([0, 2])
ax[2].invert_yaxis()
ax[2].legend(loc = "upper right", fontsize = 12)

plt.tight_layout()

plt.savefig("../../figures/1_refractionAnalysis.png", dpi = 200)
plt.show()
