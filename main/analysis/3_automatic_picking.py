import timeit
import numpy as np
import matplotlib.pyplot as plt

from all_functions import readBinaryMatrix

from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter

#--------------------------------------------------------------------
start = timeit.default_timer()

nt = 5001
dt = 1e-3
tlag = 0.15
tcut = 4.50

sline = 177
shots_all = 31329
nodes_all = 64 

node = 1
iline = 177

seismic = readBinaryMatrix(nt, shots_all, f"../../inputs/seismograms/seismogram_{nt}x{shots_all}_shot_{node}.bin")

traces = slice(int((iline-1)*(shots_all/sline)), int(iline*(shots_all/sline)))

times = slice(int(tlag/dt), int((tcut+tlag)/dt + 1))

scale = 0.9 * np.std(seismic)

seismic = seismic[times,traces]

picks = np.zeros_like(seismic[0])

window = 0.100
iw = int(window/dt)

for trace in range(len(seismic[0])):

    T = seismic[:, trace]
    A = np.zeros(len(seismic)-2*iw)
    B = np.zeros(len(seismic)-2*iw)

    for i in range(iw,len(seismic)-iw):
        A[i-iw] = np.sum(T[i-iw:i]) + 0.000000001
        B[i-iw] = np.sum(T[i:i+iw]) + 0.000000001
    
    S = np.abs((B/A) * (B-A))
    S *= (B + A) * (B / A) * (B - A) / (iw + i)**2

    S = np.abs(S)
    S *= 1.0 / np.max(S)

    index = np.where(S > 0.01)[0][0]

    picks[trace] = index + 2*(iw + 10)



picks_gauss = gaussian_filter(picks, 1.5)
picks_savgol = savgol_filter(picks, 15, 2)
picks_mean = (picks_gauss + picks_savgol) / 2

picks_mean[0] = picks_savgol[0]
picks_mean[-1] = picks_savgol[-1]

stop = timeit.default_timer()
print('Run time: ', stop - start)  

tloc = np.linspace(0, len(seismic)-1, 11, dtype = int)
tlab = np.around(np.linspace(0, len(seismic)-1, 11) * dt, decimals = 1)

xloc = np.linspace(0, len(seismic[0])-1, 11, dtype = int)
xlab = np.linspace(0, len(seismic[0]), 11, dtype = int)

plt.figure(1, figsize = (10, 7))
plt.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)
plt.plot(picks, label = "Original pick")
plt.plot(picks_gauss, label = "Gaussian smoothing")
plt.plot(picks_savgol, label = "Savgol smoothing")
plt.plot(picks_mean, label = "Mean using smoothed picks")

plt.xlabel("Traces", fontsize = 15)
plt.ylabel("Time [s]", fontsize = 15)

plt.legend(loc = "upper right", fontsize = 12)

plt.yticks(tloc, tlab)
plt.xticks(xloc, xlab)

plt.tight_layout()
plt.show()
