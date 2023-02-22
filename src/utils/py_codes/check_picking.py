import sys
import timeit
import numpy as np
import matplotlib.pyplot as plt

from functions import *

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

#--------------------------------------------------------------------

parameters = sys.argv[1]    

nt = int(catch_parameter(parameters,"nt"))
traces = int(catch_parameter(parameters,"gather_traces"))

nodes_all = int(catch_parameter(parameters,"nodes_all")) 
shots_all = int(catch_parameter(parameters,"shots_all"))

pick_folder = catch_parameter(parameters,"pick_folder")[1:-1]
data_folder = catch_parameter(parameters,"data_folder")[1:-1]

dt = float(catch_parameter(parameters,"dt"))

fw = int(catch_parameter(parameters,"qc_filter_window"))

enable_qc = eval(catch_parameter(parameters,"enable_qc"))

if enable_qc:
    start = timeit.default_timer()

    for node in range(nodes_all): 

        rawPicks_all = readBinaryArray(shots_all,f"{pick_folder}rawPicks_shot_{node+1}_{shots_all}_samples.bin")
        output_picks = np.zeros(shots_all)

        for line in range(traces):

            window = slice(int(line*shots_all/traces),int((line+1)*shots_all/traces))

            rawPicks = rawPicks_all[window].copy()    

            rawPicks[int(fw):-int(fw)] = median_filter(rawPicks[int(fw):-int(fw)], fw)

            output_picks[window] = savgol_filter(rawPicks, fw, 2)     

        print(f"File {pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin was written successfully.")
        output_picks.astype("float32", order = "F").tofile(f"{pick_folder}obsData_{shots_all}_samples_shot_{node+1}.bin")

    stop = timeit.default_timer()
    print('Run time: ', stop - start)  

# Quality control 

check_data = eval(catch_parameter(parameters,"check_data"))

if check_data: 

    tlag = float(catch_parameter(parameters,"tlag"))
    tcut = float(catch_parameter(parameters,"tcut"))

    updated_nt = int(tcut/dt)+1

    times = slice(int(tlag/dt), int((tlag+tcut)/dt))

    tloc = np.linspace(0, updated_nt-1, 11, dtype = int)
    tlab = np.around(tloc * dt, decimals = 1)

    xloc = np.linspace(0, traces-1, 11, dtype = int)
    xlab = np.linspace(0, traces, 11, dtype = int)

    igather = int(catch_parameter(parameters,"igather"))
    fgather = int(catch_parameter(parameters,"fgather"))
    dgather = int(catch_parameter(parameters,"dgather"))

    gathers = np.arange(igather, fgather+dgather, dgather, dtype = int)

    for gather in gathers:
        seismic_all = readBinaryMatrix(nt, shots_all, f"{data_folder}seismogram_{nt}x{shots_all}_shot_{gather}.bin") 
        rawPicks_all = readBinaryArray(shots_all, f"{pick_folder}rawPicks_shot_{gather}_{shots_all}_samples.bin")
        output_picks = readBinaryArray(shots_all, f"{pick_folder}obsData_{shots_all}_samples_shot_{gather}.bin")

        seismic = seismic_all[times, :]    

        scale = 0.5 * np.std(seismic)

        plt.figure(1, figsize = (15, 5))
        plt.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
        plt.plot(rawPicks_all/dt, "o")
        plt.plot(output_picks/dt, "o")
        plt.title(f"Picked seismogram for node {gather}", fontsize= 18)
        plt.ylabel("Time [s]", fontsize = 15)
        plt.xlabel("Trace index", fontsize = 15)
        plt.yticks(tloc,tlab)
        plt.xlim([0, shots_all])
        plt.tight_layout()
        plt.show()

