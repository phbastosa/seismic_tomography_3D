import numpy as np
import matplotlib.pyplot as plt

files = [["pointShot", "pointNode"], ["xLineShots", "xLineNodes"], 
         ["yLineShots", "yLineNodes"], ["carpetShots", "carpetNodes"], 
         ["carpetShotsReciprocity", "carpetNodesReciprocity"],  
         ["circularShots", "circularNodes"]]

fig, axs = plt.subplots(2,3, figsize = (13, 8))

ind = 0
for i in range(len(axs)):
    for j in range(len(axs[0])):
        [xs, ys, zs] = np.loadtxt(f"outputs/{files[ind][0]}.txt", delimiter = ",", unpack = True)    
        [xr, yr, zr] = np.loadtxt(f"outputs/{files[ind][1]}.txt", delimiter = ",", unpack = True)    

        axs[i,j].scatter(xs, ys, c = 'green')
        axs[i,j].scatter(xr, yr, c = 'black')
    
        axs[i,j].set_xlim([0,10000])
        axs[i,j].set_ylim([0,10000])

        ind += 1

plt.tight_layout()
plt.show()