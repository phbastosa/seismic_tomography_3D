import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

def readBinaryMatrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)    
    return np.reshape(data, [n1,n2], order='F')

def pick_first_break(window, dt, trace):

    iw = int(window / dt)

    A = np.zeros(len(trace))
    B = np.zeros(len(trace))
    S = np.zeros(len(trace))

    for t in range(iw, len(trace)-iw):
        A[t] = np.sum(trace[t - int(iw/2):t])
        B[t] = np.sum(trace[t:t + int(iw/2)])

        if A[t] != 0:
            S[t] = np.abs((B[t] / A[t]) * (B[t] - A[t])) 
            S[t] *= (B[t] + A[t]) * (B[t] - A[t]) /  (iw + t)**2

    S *= 1 / np.max(S)    

    return S

#---------------------------
# plot dos sismogramas 

lines = 10000
nt = 3001
dt = 1e-3

seismic = readBinaryMatrix(nt, lines,"seismogram_3001x10000_shot_1.bin")

perc = 5.0 * np.std(seismic)

npick = 100

panels = [seismic[:,:npick],
          seismic[:,30*npick:31*npick],
          seismic[:,60*npick:61*npick]]

titles = ["Offset curto", "Offset médio", "Offset longo"]

xloc = np.linspace(0, npick-1, 9)
xlab = np.linspace(0, npick, 9, dtype = int)

fig, ax = plt.subplots(1,3, figsize = (15,6))
for k in range(len(panels)):
    ax[k].imshow(panels[k], aspect = "auto", cmap = "Greys", vmin = -perc, vmax = perc)
    ax[k].plot(np.ones(nt)*npick*0.5, np.arange(nt), "-b", label = "Traço analisado")    
    ax[k].set_title(titles[k], fontsize = 18)
    ax[k].set_xlabel("Índice do traço", fontsize = 15)
    ax[k].set_ylabel("Tempo [ms]", fontsize = 15)

    ax[k].set_xticks(xloc)
    ax[k].set_xticklabels(xlab)

    ax[k].legend()

plt.tight_layout()
plt.show()

# Teste de pick nos traços

traces = [panels[0][:,51].copy(), 
          panels[1][:,51].copy(), 
          panels[2][:,51].copy()]

times = np.arange(nt) * dt

windows = np.arange(0.005, 0.030, 0.005)

fig, ax = plt.subplots(1, 3, figsize = (12,7))

picks = np.zeros(len(traces))

for k in range(len(traces)):
    traces[k] *= 1 / np.max(traces[k])
    
    ax[k].plot(traces[k], times)

    for w in windows:
        S = pick_first_break(w, dt, traces[k])
        
        for s in range(len(S)):
            if S[s] > 1e-10:
                picks[k] = s
                break

        ax[k].plot(traces[k][int(picks[k])], dt*int(picks[k]), "o", label = f"Janela de {w*1e3:.1f} ms")

    ax[k].plot(np.zeros(nt), np.arange(nt)*dt, "k")    
    ax[k].set_ylim(0,3)
    ax[k].invert_yaxis()
    ax[k].legend(loc = "lower left")

    ax[k].set_title(titles[k], fontsize = 18)
    ax[k].set_xlabel("Amplitude", fontsize = 15)
    ax[k].set_ylabel("Tempo [s]", fontsize = 15)

plt.tight_layout()
plt.show()

# pick completo

median = 3
savgol = 11

final_pick = []

perc = 5 * np.std(seismic)

fig, ax = plt.subplots(1,3, figsize = (15,6))

for k in range(len(panels)):
    
    picks = np.zeros(npick)
    for t in range(npick):
        trace = panels[k][:,t] * 1 / np.max(panels[k][:,t])

        S = pick_first_break(0.025, dt, trace)   

        for s in range(len(S)):
            if S[s] > 1e-6:
                picks[t] = s
                break

    ax[k].imshow(panels[k], aspect = "auto", cmap = "Greys", vmin = -perc, vmax = perc)
    ax[k].plot(picks, label = "Pick bruto")

    # quality control

    picks[median:-median] = median_filter(picks[median:-median], median)

    ax[k].plot(picks, label = "Pick com filtro mediana")
    
    picks = savgol_filter(picks.copy(), savgol, 2)
    picks = savgol_filter(picks.copy(), savgol, 2)

    ax[k].plot(picks, label = "Pick suavizado")

    ax[k].set_title(titles[k], fontsize = 18)
    ax[k].set_xlabel("Índice do traço", fontsize = 15)
    ax[k].set_ylabel("Tempo [ms]", fontsize = 15)

    ax[k].set_xticks(xloc)
    ax[k].set_xticklabels(xlab)

    final_pick.append(picks)

    ax[k].legend()

plt.tight_layout()
plt.show()

