import numpy as np
import matplotlib.pyplot as plt

from all_functions import *

from scipy.ndimage import gaussian_filter

# General recovered models plot 
# Figure with 12 subfigures 

shots = np.loadtxt("../../inputs/geometry/shots.txt", delimiter=",")
nodes = np.loadtxt("../../inputs/geometry/nodes.txt", delimiter=",")

trueModel = readBinaryVolume(105,401,401,f"../../inputs/models/trueModel_105x401x401_12.5m.bin")
trueModel = trueModel[::2,::2,::2]

nz = 53
nx = 201
ny = 201 

dh = 25

initModel = readBinaryVolume(nz,nx,ny,f"../../inputs/models/initModel_{nz}x{nx}x{ny}_{dh}m.bin")

pod_models = np.zeros((3,nz,nx,ny))
pod_models[0,:,:] = initModel
pod_models[-1,:,:] = trueModel

fim_models =  np.zeros((3,nz,nx,ny))
fim_models[0,:,:] = initModel
fim_models[-1,:,:] = trueModel

fsm_models =  np.zeros((3,nz,nx,ny))
fsm_models[0,:,:] = initModel
fsm_models[-1,:,:] = trueModel

pod_models[1,:,:,:] = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/pod_estimatedModel_iteration_10.bin")
pod_models[1,:,:,:] = 1.0 / gaussian_filter(1.0 / pod_models[1,:,:,:], 3.0)
 
fim_models[1,:,:,:] = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/fim_estimatedModel_iteration_10.bin")
fim_models[1,:,:,:] = 1.0 / gaussian_filter(1.0 / fim_models[1,:,:,:], 3.0)

fsm_models[1,:,:,:] = readBinaryVolume(nz,nx,ny,f"../../outputs/recoveredModels/fsm_estimatedModel_iteration_10.bin")
fsm_models[1,:,:,:] = 1.0 / gaussian_filter(1.0 / fsm_models[1,:,:,:], 3.0)

slices = np.array([31,100,100]) # XY, ZX, ZY
subplots = np.array([1,3])

multiBoxPlot(pod_models,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/podvin_result.png", dpi = 200)
plt.show()

multiBoxPlot(fim_models,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/jeong_result.png", dpi = 200)
plt.show()

multiBoxPlot(fsm_models,shots,nodes,dh,slices,subplots)
plt.savefig("../../figures/noble_result.png", dpi = 200)
plt.show()

#---------------------------------------------------------------------------------------------







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
plt.savefig("../../figures/convergencia.png", dpi = 200)
plt.show()




