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

start = timeit.default_timer()

for node in range(nodes_all): 

    rawPicks_all = readBinaryArray(shots_all,f"{pick_folder}rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = np.zeros(shots_all)

    diff = rawPicks_all[1:] - rawPicks_all[:-1]  

    outliers = np.where(diff > 0.2)[0]

    for k in outliers:
        rawPicks_all[k-2:k+2] = median_filter(rawPicks_all[k-2:k+2], 5)

    for line in range(traces):

        window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

        rawPicks = rawPicks_all[window].copy()    

        output_picks[window] = savgol_filter(rawPicks, fw, 2)     

    print(f"File {pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin was written successfully.")
    output_picks.astype("float32", order = "F").tofile(f"{pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin")

stop = timeit.default_timer()
print('Run time: ', stop - start)  

# Quality control 

tloc = np.linspace(0, nt-1, 11, dtype = int)
tlab = np.around(tloc * dt, decimals = 1)

xloc = np.linspace(0, traces-1, 11, dtype = int)
xlab = np.linspace(0, traces, 11, dtype = int)

gather = 100

seismic = readBinaryMatrix(nt, shots_all, f"{data_folder}seismogram_{nt}x{shots_all}_shot_{gather+1}.bin") 
input_picks = readBinaryArray(shots_all, f"{pick_folder}rawPicks_shot_{gather+1}_{shots_all}_samples.bin")
output_picks = readBinaryArray(shots_all, f"{pick_folder}obsData_{shots_all}_samples_shot_{gather+1}.bin")

scale = 0.001 * np.std(seismic)

fig, ax =  plt.subplots(1,1, figsize = (15, 5))

ax.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
ax.plot(input_picks/dt, "o")
ax.plot(output_picks/dt, "o")
ax.set_title(f"Picked seismogram for node {gather+1}", fontsize= 18)
ax.set_ylabel("Time [s]", fontsize = 15)
ax.set_xlabel("Trace index", fontsize = 15)
ax.set_yticks(tloc)
ax.set_yticklabels(tlab)
ax.set_xlim([0, shots_all])

plt.tight_layout()
plt.show()
