import timeit
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

def readBinaryArray(n,filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

#--------------------------------------------------------------------
start = timeit.default_timer()

nt = 3101
traces = 100
nodes_all = 441 
shots_all = 10000

dt = 1e-3

fw = 11

for node in range(nodes_all): 

    rawPicks_all = readBinaryArray(shots_all, f"../../inputs/picks/rawPicks_shot_{node+1}_{shots_all}_samples.bin")
    output_picks = np.zeros(shots_all)

    for line in range(traces):

        window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

        rawPicks = rawPicks_all[window].copy()    

        rawPicks[int(fw):-int(fw)] = median_filter(rawPicks[int(fw):-int(fw)], fw)

        output_picks[window] = savgol_filter(rawPicks, fw, 2)     

    print(f"File ../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin was written successfully.")
    output_picks.astype("float32", order = "F").tofile(f"../../inputs/picks/obsData_{shots_all}_samples_shot_{node+1}.bin")

stop = timeit.default_timer()
print('Run time: ', stop - start)  

# Quality control 

tlag = 0.10
tcut = 2.80

updated_nt = 2.80

times = slice(int(tlag/dt), int((tlag+tcut)/dt))

tloc = np.linspace(0, int(updated_nt/dt)-1, 11, dtype = int)
tlab = np.around(tloc * dt, decimals = 1)

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
    plt.title(f"Picked seismogram for node {node+1}", fontsize= 18)
    plt.ylabel("Time [s]", fontsize = 15)
    plt.xlabel("Trace index", fontsize = 15)
    plt.yticks(tloc,tlab)
    plt.xlim([0, shots_all])
    plt.tight_layout()
    plt.show()

