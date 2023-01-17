import timeit
import numpy as np
import matplotlib.pyplot as plt

from all_functions import readBinaryMatrix, readBinaryArray, irls

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

#--------------------------------------------------------------------
start = timeit.default_timer()

nt = 5001
traces = 177
nodes_all = 225 
shots_all = 31329

dt = 1e-3

fw = 11

for node in range(38, 39): 

    seismic_all = readBinaryMatrix(nt, shots_all, f"../../inputs/seismograms/seismogram_{nt}x{shots_all}_shot_{node+1}.bin")
    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = np.zeros(shots_all)

    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")

    for line in range(118, 119):

        window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

        rawPicks = rawPicks_all[window].copy()    

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

    # seismic_all = readBinaryMatrix(nt, shots_all, f"../../inputs/seismograms/seismogram_{nt}x{shots_all}_shot_{node+1}.bin")
    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = readBinaryArray(shots_all, f"../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin")

    # scale = 0.9 * np.std(seismic_all)

    plt.figure(1, figsize = (18, 5))
    plt.plot(rawPicks_all)
    plt.plot(output_picks)
    plt.title(f"{node+1}")
    plt.ylim([0, 4.5])
    plt.gca().invert_yaxis()
    plt.show()

    # for line in range(traces):
        
    #     trace = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

    #     seismic = seismic_all[times, trace]

    #     plt.figure(1, figsize = (10, 7))

    #     plt.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
    #     plt.plot(rawPicks_all[trace]/dt, label = "Original picks")
    #     plt.plot(output_picks[trace]/dt, label = "Output picks")

    #     plt.xlabel("Traces", fontsize = 15)
    #     plt.ylabel("Time [s]", fontsize = 15)
    #     plt.title(f"Receiver gather {node+1} offset {line}", fontsize = 18)

    #     plt.legend(loc = "upper right", fontsize = 12)

    #     plt.yticks(tloc, tlab)
    #     plt.xticks(xloc, xlab)

    #     plt.tight_layout()
    #     plt.show()
