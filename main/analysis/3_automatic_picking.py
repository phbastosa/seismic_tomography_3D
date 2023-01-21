import timeit
import numpy as np
import matplotlib.pyplot as plt

from all_functions import readBinaryMatrix, readBinaryArray, irls

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

#--------------------------------------------------------------------
start = timeit.default_timer()

nt = 3001
traces = 100
nodes_all = 121 
shots_all = 10000

dt = 1e-3

fw = 11

for node in range(nodes_all): 

    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = np.zeros(shots_all)

    for line in range(traces):

        window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

        rawPicks = rawPicks_all[window].copy()    

        rawPicks[int(fw/2):-int(fw/2)] = median_filter(rawPicks[int(fw/2):-int(fw/2)], fw)

        output_picks[window] = savgol_filter(rawPicks, fw, 2)     

    print(f"File ../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin was written successfully.")
    output_picks.astype("float32", order = "F").tofile(f"../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin")

stop = timeit.default_timer()
print('Run time: ', stop - start)  

# Quality control 

tlag = 0.15
tcut = 4.50

updated_nt = 4501

times = slice(int(tlag/dt), int((tlag+tcut)/dt))

tloc = np.linspace(0, updated_nt-1, 11, dtype = int)
tlab = np.around(np.linspace(0, updated_nt-1, 11) * dt, decimals = 1)

xloc = np.linspace(0, traces-1, 11, dtype = int)
xlab = np.linspace(0, traces, 11, dtype = int)

for node in range(nodes_all):
    seismic_all = readBinaryMatrix(nt, shots_all, f"../../inputs/seismograms/seismogram_{nt}x{shots_all}_shot_{node+1}.bin") 
    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = readBinaryArray(shots_all, f"../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin")

    seismic = seismic_all[times, :]    

    scale = 0.5 * np.std(seismic)

    plt.figure(1, figsize = (20, 5))
    plt.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
    plt.plot(rawPicks_all/dt, "o")
    plt.plot(output_picks/dt, "o")
    plt.title(f"{node+1}")
    plt.xlim([0, shots_all])
    plt.tight_layout()
    plt.show()

