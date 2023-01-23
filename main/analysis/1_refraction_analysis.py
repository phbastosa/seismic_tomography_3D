import numpy as np
import matplotlib.pyplot as plt

from all_functions import readBinaryArray, readBinaryMatrix 
from all_functions import readBinaryVolume, analyticalRefraction
from scipy.ndimage import gaussian_filter

nx = 401
ny = 401
nz = 21

dx = 12.5
dy = 12.5
dz = 12.5

v = np.array([2000,3500])
z = np.array([200]) 

refractiveModel = np.zeros((nz,nx,ny))

refractiveModel[:int(z[0]/dz),:,:] = v[0]
refractiveModel[int(z[0]/dz):,:,:] = v[1]

refractiveModel.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/refractiveModel_{nz}x{nx}x{ny}_{dz}m.bin")

sx,sy,sz = np.loadtxt("../../inputs/geometry/shots.txt", delimiter = ",",unpack=True)
rx,ry,rz = np.loadtxt("../../inputs/geometry/nodes.txt", delimiter = ",",unpack=True)

#--------------------------------

dt = 1e-3
nt = 3001
traces = 161

seismic = readBinaryMatrix(nt, traces, "../../inputs/seismograms/seismogram_3001x161_shot_1.bin") 

pod = readBinaryArray(traces, "../../outputs/dcal/refracStudy_pod_times_nr161_shot_1.bin")
fim = readBinaryArray(traces, "../../outputs/dcal/refracStudy_fim_times_nr161_shot_1.bin")
fsm = readBinaryArray(traces, "../../outputs/dcal/refracStudy_fsm_times_nr161_shot_1.bin")

offset = np.sqrt((sx - rx)**2 + (sy - ry)**2)

td, tr = analyticalRefraction(v, z, offset)

delay = 0.1 
tcut = 2.5

seismic = seismic[int(delay/dt-1):,:]
seismic = seismic[:int(tcut/dt+1),:]

nt_updated = len(seismic)

for k in range(traces):
    seismic[int(tr[0,k]/dt + 50):, k] = np.zeros(nt_updated - int(tr[0,k]/dt + 50))
    seismic[:,k] = gaussian_filter(seismic[:,k], 5.0) 

perc = 0.05 * np.std(seismic)

time = np.arange(len(seismic))

plt.figure(2, figsize = (8, 9))

plt.imshow(seismic, aspect="auto", cmap="Greys", vmin=-perc, vmax=perc)
plt.plot(tr[0]/dt)
plt.plot(pod/dt)
plt.plot(fim/dt)
plt.plot(fsm/dt)

plt.tight_layout()
plt.show()
