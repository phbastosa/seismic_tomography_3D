import timeit
import numpy as np
import matplotlib.pyplot as plt

from functions import *

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

#--------------------------------------------------------------------

nt = 3001
traces = 140

nodes_all = 11 * 16 
shots_all = 100 * 140

pick_folder = "../../../inputs/picks/"
data_folder = "../../../inputs/seismic_data/"

dt = 0.001

fw = 11

gather_outlier = [27,110]

# start = timeit.default_timer()

# for node in range(nodes_all): 

#     rawPicks_all = readBinaryArray(shots_all,f"{pick_folder}rawPicks_shot_{node+1}_{shots_all}_samples.bin")
#     output_picks = np.zeros(shots_all)

#     if node in gather_outlier:
#         diff = rawPicks_all[1:] - rawPicks_all[:-1]  
#         outliers = np.where(diff > 1.0)[0]
#         for k in outliers:
#             rawPicks_all[k-2:k+2] = median_filter(rawPicks_all[k-2:k+2], 5)

#     for line in range(traces):

#         window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

#         rawPicks = rawPicks_all[window].copy()    

#         output_picks[window] = savgol_filter(rawPicks, fw, 2)     

#     print(f"File {pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin was written successfully.")
#     output_picks.astype("float32", order = "F").tofile(f"{pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin")

# stop = timeit.default_timer()
# print('Run time: ', stop - start)  

# Quality control 

tloc = np.linspace(0, nt-1, 11, dtype = int)
tlab = np.around(tloc * dt, decimals = 1)

xloc = np.linspace(0, traces-1, 11, dtype = int)
xlab = np.linspace(0, traces, 11, dtype = int)

gather = 1 # 1 - 176
group  = 10 # 1 - 140

seismic = readBinaryMatrix(nt, shots_all, f"{data_folder}seismogram_{nt}x{shots_all}_shot_{gather}.bin") 
input_picks = readBinaryArray(shots_all, f"{pick_folder}rawPicks_shot_{gather}_{shots_all}_samples.bin")
output_picks = readBinaryArray(shots_all, f"{pick_folder}obsData_{shots_all}_samples_shot_{gather}.bin")

scale = 5.0 * np.std(seismic)

fig, ax =  plt.subplots(1,1, figsize = (20, 7))

ax.imshow(seismic[:,(group-1)*100:group*100], aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax.plot(input_picks[(group-1)*100:group*100]/dt, "o")
ax.plot(output_picks[(group-1)*100:group*100]/dt, "o")

# ax.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
# ax.plot(input_picks/dt, "o")
# ax.plot(output_picks/dt, "o")

ax.set_title(f"Sismograma {gather} - linha {group} de {traces}", fontsize = 18)
ax.set_ylabel("Tempo [s]", fontsize = 15)
ax.set_xlabel("Índice do traço", fontsize = 15)
ax.set_yticks(tloc)
ax.set_yticklabels(tlab)

ax.set_xticks(np.linspace(0,99, 7, dtype = int))
ax.set_xticklabels(np.linspace((group-1)*100,group*100, 7, dtype = int))

plt.tight_layout()
plt.show()
