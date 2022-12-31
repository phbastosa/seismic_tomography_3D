import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy import signal
from matplotlib import patches
from matplotlib.figure import Figure
from matplotlib import rcParams

#-------------------------------------------------------------------------------------------

def rickerGenerator(fmax,nsrc,dt):    
    ricker = np.zeros(nsrc)
    fc = fmax/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1 - 2 * np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
        aux2 = np.exp(-np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))
        ricker[i + s] = aux1 * aux2 

    return ricker

def zplane(z,p,k,subplot):

    ax = plt.subplot(subplot)
    uc = patches.Circle((0,0), radius=1, fill=False, color='black', ls='dashed')
    ax.add_patch(uc)
    
    # Plot the zeros and set marker properties    
    t1 = plt.plot(z.real, z.imag, 'go', ms=10)
    plt.setp( t1, markersize=10.0, markeredgewidth=1.0, markeredgecolor='k', markerfacecolor='g')

    # Plot the poles and set marker properties
    t2 = plt.plot(p.real, p.imag, 'rx', ms=10)
    plt.setp( t2, markersize=12.0, markeredgewidth=3.0, markeredgecolor='r', markerfacecolor='r')

    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_title("Diagrama de polos e zeros")
    
    # set the ticks
    r = 1.5; plt.axis('scaled'); plt.axis([-r, r, -r, r])
    ticks = [-1, -.5, .5, 1]; plt.xticks(ticks); plt.yticks(ticks)

def fourierInterpolation(signal, dt, new_dt):

    ns = len(signal)

    periods = int(dt / new_dt)

    freqs = np.fft.fftfreq(ns, dt)
    fft_signal = np.fft.fft(signal)

    new_fft_signal = np.insert(fft_signal, int(ns / 2), np.zeros(periods * ns)) 

    new_ns = periods * ns + ns

    scale_factor = new_ns / ns
    output_signal = scale_factor * np.real(np.fft.ifft(new_fft_signal))
    
    return output_signal

#-------------------------------------------------------------------------------------------

nsrc = 30
dt = 0.01
fmax = 30

ricker = rickerGenerator(fmax, nsrc, dt)

#-------------------------------------------------------------------------------------------

a = [1]
b = ricker

z, p, k = signal.tf2zpk(b, a)

#-------------------------------------------------------------------------------------------

j = complex(0,1)

pn = np.arange(nsrc,dtype=int) 

fftRicker = np.fft.fft(ricker)
fftRicker *= 1 / np.max(np.abs(fftRicker))

freqs = np.fft.fftfreq(nsrc,dt)
mask = freqs >= 0

#-------------------------------------------------------------------------------------------

distance = np.abs(z)
angle = np.angle(z)

position = np.where(distance > 1)

distance[position] = 1/distance[position]

z_in = np.zeros(len(z),dtype=complex)

z_in[position] = distance[position]*np.cos(angle[position]) + j*distance[position]*np.sin(angle[position]) 

#-------------------------------------------------------------------------------------------

b_in, a_in = signal.zpk2tf(z_in, p, k)

system = signal.TransferFunction(b_in, a_in, dt=dt)

_ , reconst = signal.dfreqresp(system,n=int(nsrc/2))

fftRickerReconst = np.zeros(nsrc,dtype=complex)

fftRickerReconst[int(nsrc/2):] = reconst[::-1]
fftRickerReconst[:int(nsrc/2)] = reconst[:]

rickerReconst = np.real(np.fft.ifft(reconst,n=nsrc))

rickerReconst *= 1 / np.max(np.abs(rickerReconst))
fftRickerReconst *= 1 / np.max(np.abs(fftRickerReconst))

#-------------------------------------------------------------------------------------------

plt.figure(1,figsize=(12,6))
xmask = np.linspace(pn[0],pn[-1]+1,11)

plt.subplot(221)
plt.stem(pn,ricker)
plt.xticks(xmask,(xmask-(pn[-1]+1)/2)*dt)
plt.title("Sinal de fase zero", fontsize = 18)
plt.xlabel("Tempo [s]", fontsize = 15)
plt.ylabel("Amplitude", fontsize = 15)
plt.xlim([1, nsrc-1])

zplane(z, p, k, 222)
plt.title("Diagrama de polos e zeros", fontsize = 18)

plt.subplot(223)
plt.stem(freqs[mask], np.abs(fftRicker[mask]))
plt.xlim([0,fmax])
plt.title("Espectro de amplitudes", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 15)
plt.ylabel("Amplitude", fontsize = 15)

plt.subplot(224)
plt.plot(freqs[mask],np.angle(fftRicker[mask],deg=True),"o")
plt.xlim([0,fmax])
plt.ylim([-180,180])
plt.title("Espectro de fase", fontsize = 15)
plt.xlabel("Frequência [Hz]", fontsize = 15)
plt.ylabel("Ângulo [°]", fontsize = 15)

plt.tight_layout()
plt.savefig("../figures/4_zeroPhaseRicker.png", dpi = 200)
plt.show(block=False)

#--------------------------------------------------------------

plt.figure(2,figsize=(12,6))
plt.subplot(221)
plt.stem(pn,rickerReconst)
plt.xticks(xmask,xmask*dt)
plt.title("Sinal de fase mínima", fontsize = 18)
plt.xlabel("Tempo [s]", fontsize = 15)
plt.ylabel("Amplitude", fontsize = 15)
plt.xlim([0, 18])

zplane(z_in, p, k, 222)
plt.title("Diagrama de polos e zeros", fontsize = 18)

plt.subplot(223)
plt.stem(freqs[mask],np.abs(fftRickerReconst[mask]))
plt.xlim([0,fmax])
plt.title("Espectro de amplitudes", fontsize = 15)
plt.xlabel("Frequência [Hz]", fontsize = 15)
plt.ylabel("Amplitude", fontsize = 15)

plt.subplot(224)
plt.plot(freqs[mask],np.angle(fftRickerReconst[mask],deg=True),"o")
plt.xlim([0,fmax])
plt.ylim([-180,180])
plt.title("Espectro de fase", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 15)
plt.ylabel("Ângulo [°]", fontsize = 15)

plt.tight_layout()
plt.savefig("../figures/5_minPhaseRicker.png", dpi = 200)
 
plt.show(block=False)

#-------------------------------------------------------------------------------------------

new_dt = 0.001

zRicker = fourierInterpolation(ricker, dt, new_dt)
mRicker = fourierInterpolation(rickerReconst, dt, new_dt)

new_ns = len(zRicker)

zRicker.astype("float32", order = "F").tofile(f"../inputs/wavelets/sourceZeroPhase_{new_ns}_{fmax:.0f}Hz_{new_dt*1000:.0f}ms.bin")
mRicker.astype("float32", order = "F").tofile(f"../inputs/wavelets/sourceMinPhase_{new_ns}_{fmax:.0f}Hz_{new_dt*1000:.0f}ms.bin")

new_ns += 3000

full_zRicker = np.zeros(new_ns)
full_mRicker = np.zeros(new_ns)

full_zRicker[:len(zRicker)] = zRicker
full_mRicker[:len(mRicker)] = mRicker

t = np.arange(new_ns) * new_dt

zfRicker = np.fft.fft(full_zRicker)
mfRicker = np.fft.fft(full_mRicker)
freq = np.fft.fftfreq(new_ns, new_dt)
mask = freq >= 0

fig, ax = plt.subplots(2,2, figsize = (12, 6))

ax[0][0].plot(t[:len(zRicker)], zRicker)
ax[0][0].set_xticks(np.linspace(t[0],t[len(zRicker)-1],11))
ax[0][0].set_xticklabels(np.around(np.linspace(t[0]-0.16,t[len(zRicker)-1]-0.17,11), decimals = 2))
ax[0][0].set_xlim([0,t[len(zRicker)-1]])
ax[0][0].set_ylim([-1.5,1.5])
ax[0][0].set_title("Fonte fase zero", fontsize = 18)
ax[0][0].set_xlabel("Tempo [s]", fontsize = 15)
ax[0][0].set_ylabel("Amplitude", fontsize = 15)

ax[0][1].plot(freq[mask], np.abs(zfRicker[mask]))
ax[0][1].set_xlim([0,fmax])
ax[0][1].set_title("Espectro da fonte fase zero", fontsize = 18)
ax[0][1].set_xlabel("Frequência [Hz]", fontsize = 15)
ax[0][1].set_ylabel("Amplitude", fontsize = 15)

ax[1][0].plot(t[:len(zRicker)], mRicker)
ax[1][0].set_xticks(np.linspace(t[0],t[len(zRicker)-1],11))
ax[1][0].set_xticklabels(np.around(np.linspace(t[0],t[len(zRicker)-1],11), decimals = 2))
ax[1][0].set_xlim([0,t[len(zRicker)-1]])
ax[1][0].set_ylim([-1.5,1.5])
ax[1][0].set_title("Fonte fase mínima", fontsize = 18)
ax[1][0].set_xlabel("Tempo [s]", fontsize = 15)
ax[1][0].set_ylabel("Amplitude", fontsize = 15)

ax[1][1].plot(freq[mask], np.abs(mfRicker[mask]))
ax[1][1].set_xlim([0,fmax])
ax[1][1].set_title("Espectro da fonte fase mínima", fontsize = 18)
ax[1][1].set_xlabel("Frequência [Hz]", fontsize = 15)
ax[1][1].set_ylabel("Amplitude", fontsize = 15)

plt.tight_layout()
plt.savefig("../figures/7_targetWavelets.png", dpi = 200)
plt.show(block = False)