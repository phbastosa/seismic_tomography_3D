import numpy as np
import matplotlib.pyplot as plt

from all_functions import *

from scipy.ndimage import gaussian_filter

# General recovered models plot 

shots = np.loadtxt("../../inputs/geometry/shots.txt", delimiter=",")
nodes = np.loadtxt("../../inputs/geometry/nodes.txt", delimiter=",")

trueModel = readBinaryVolume(105,401,401,f"../../inputs/models/trueModel_105x401x401_12.5m.bin")
trueModel = trueModel[::2,::2,::2]

nz = 53
nx = 201
ny = 201 

dh = 25

mask_z_beg = 200
mask_z_end = 1000

initModel = readBinaryVolume(nz,nx,ny,f"../../inputs/models/initModel_{nz}x{nx}x{ny}_{dh}m.bin")

pod_model = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/pod_estimatedModel_iteration_20.bin")
pod_model[:int(mask_z_beg/dh),:,:] = initModel[:int(mask_z_beg/dh),:,:]
pod_model[int(mask_z_end/dh):,:,:] = initModel[int(mask_z_end/dh):,:,:]
pod_model = 1.0 / gaussian_filter(1.0 / pod_model, 4.0)
pod_model.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/pod_initModel_{nz}x{nx}x{ny}_{dh}m.bin")

fim_model = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/fim_estimatedModel_iteration_20.bin")
fim_model[:int(mask_z_beg/dh),:,:] = initModel[:int(mask_z_beg/dh),:,:]
fim_model[int(mask_z_end/dh):,:,:] = initModel[int(mask_z_end/dh):,:,:]
fim_model[:,:,:] = 1.0 / gaussian_filter(1.0 / fim_model, 4.0)
fim_model.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/fim_initModel_{nz}x{nx}x{ny}_{dh}m.bin")

fsm_model = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/fsm_estimatedModel_iteration_20.bin")
fsm_model[:int(mask_z_beg/dh),:,:] = initModel[:int(mask_z_beg/dh),:,:]
fsm_model[int(mask_z_end/dh):,:,:] = initModel[int(mask_z_end/dh):,:,:]
fsm_model[:,:,:] = 1.0 / gaussian_filter(1.0 / fsm_model, 4.0)
fsm_model.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/fsm_initModel_{nz}x{nx}x{ny}_{dh}m.bin")

slices = np.array([31,100,100]) # XY, ZX, ZY
subplots = np.array([1,1])

multiBoxPlot(pod_model,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/podvin_result.png", dpi = 200)
plt.show()

multiBoxPlot(fim_model,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/jeong_result.png", dpi = 200)
plt.show()

multiBoxPlot(fsm_model,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/noble_result.png", dpi = 200)
plt.show()

#---------------------------------------------------------------------------------------------
# Trace analysis

npoints = 5
ntraces = 5

xc = np.array([2000, 3000, 2500, 2000, 3000]) / dh 
yc = np.array([2000, 2000, 2500, 3000, 3000]) / dh

depth = np.arange(nz) * dh
v_logs = np.zeros((npoints, ntraces, nz))

for point in range(npoints):
    v_logs[point, 0, :] = trueModel[:,int(xc[point]),int(yc[point])]
    v_logs[point, 1, :] = initModel[:,int(xc[point]),int(yc[point])]
    v_logs[point, 2, :] = pod_model[:,int(xc[point]),int(yc[point])]
    v_logs[point, 3, :] = fim_model[:,int(xc[point]),int(yc[point])]
    v_logs[point, 4, :] = fsm_model[:,int(xc[point]),int(yc[point])]

fig, ax = plt.subplots(1, npoints, figsize = (15,6))

labels = ["Reference", "Initial", "Classical", "Fast Iterative", "Fast Sweeping"]

for i in range(ntraces):
    for j in range(npoints):
        ax[i].plot(v_logs[i, j, :], depth, label = labels[j])
    
    ax[i].set_title(f"Gaussian {i+1}", fontsize = 18)
    ax[i].set_xlabel("Velocity [m/s]", fontsize = 15)
    ax[i].set_ylabel("Depth [m]", fontsize = 15)
    ax[i].set_xlim([1000,4000])
    ax[i].set_ylim([0,1300])

    ax[i].set_yticks(np.linspace(0,1300,11))
    ax[i].set_yticklabels(np.linspace(0,1300,11, dtype = int))

    ax[i].legend()
    ax[i].invert_yaxis()    

plt.tight_layout()
plt.savefig("../../figures/traceAnalysis.png", dpi = 200)
plt.show()

#---------------------------------------------------------------------------------------------

pod_convergency = np.loadtxt("../../outputs/convergency/pod_convergency.txt")
fim_convergency = np.loadtxt("../../outputs/convergency/fim_convergency.txt")
fsm_convergency = np.loadtxt("../../outputs/convergency/fsm_convergency.txt")

plt.figure(4, figsize = (10, 5))
plt.plot(pod_convergency, "-ob", label = "Classical kernel")
plt.plot(fim_convergency, "-oy", label = "Fast iterative kernel")
plt.plot(fsm_convergency, "-og", label = "Fast sweeping kernel")

plt.title("Convergence", fontsize = 20)
plt.xlabel("Iterations", fontsize = 15)
plt.ylabel(r"$||d_{obs} - d_{cal}||^2_2$", fontsize = 15)

plt.ylim([0,200])

plt.legend(loc = "upper right", fontsize = 15)
plt.tight_layout()
plt.savefig("../../figures/convergency.png", dpi = 200)
plt.show()




