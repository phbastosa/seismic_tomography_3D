import numpy as np
import matplotlib.pyplot as plt

from scipy import signal
from matplotlib import patches

def rickerGenerator(fcut,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fcut/(3*np.sqrt(np.pi))
    
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

def dtft(n, x):
    
	j = complex(0,1)

	w = np.linspace(-np.pi, np.pi, int(1e4))

	X = np.zeros(len(w), dtype = complex) 

	for i in range(len(n)):
		X += x[i] * np.exp(-j*w*n[i])

	return w, X

def fourierInterpolation(signal, dt, new_dt):
    
    ns = len(signal)

    periods = int(dt / new_dt)

    fft_signal = np.fft.fft(signal)

    new_fft_signal = np.insert(fft_signal, int(ns / 2), np.zeros(periods * ns)) 

    new_ns = periods * ns + ns

    scale_factor = new_ns / ns
    output_signal = scale_factor * np.real(np.fft.ifft(new_fft_signal))
    
    return output_signal

def min_phase_filter(wavelet, dt):

    nsrc = int(0.1 * len(wavelet)) + 1
    down_dt = 0.1 * dt 

    a = [1]
    b = wavelet[::10]
    j = complex(0,1)

    z, p, k = signal.tf2zpk(b, a)

    angle = np.angle(z)
    distance = np.abs(z)

    position = np.where(distance > 1.0)

    distance[position] = 1.0 / distance[position]

    z_in = np.zeros(len(z), dtype = complex)

    z_in[position] = distance[position]*np.cos(angle[position]) + j*distance[position]*np.sin(angle[position]) 

    b_in, a_in = signal.zpk2tf(z_in, p, k)

    system = signal.TransferFunction(b_in, a_in, dt = down_dt)

    _ , reconstruction = signal.dfreqresp(system, n = int(nsrc / 2))

    output_signal = np.real(np.fft.ifft(reconstruction, n = nsrc))
    output_signal[0] = 0

    output_signal = fourierInterpolation(output_signal, dt, down_dt)

    output_signal *= 1.0 / np.max(np.abs(output_signal))

    output_signal = output_signal[:len(wavelet)]

    output_signal[-int(len(wavelet)//5):] = 0

    w = 11

    output_signal[w:-w] = signal.savgol_filter(output_signal[w:-w], w, 1)

    return output_signal

nt = 151
dt = 1e-3
fmax = 40

z_ricker = rickerGenerator(fmax, nt, dt)
m_ricker = min_phase_filter(z_ricker, dt)

output = np.insert(m_ricker, nt-1, np.zeros(3001 - nt))

output *= 1e6

output.astype("float32", order = "F").tofile(f"ricker_min_phase_{len(output)}_{dt*1e3:.0f}ms.bin")

t = np.arange(nt) * dt
lag = 0.5 * (nt-1) * dt

plt.figure(1,figsize=(18,8))

plt.subplot(221)
plt.plot(t - lag, z_ricker)
plt.xlim([-lag, lag])
plt.title("Ricker de fase zero", fontsize = 18)
plt.xlabel("Tempo [s]", fontsize = 14)
plt.ylabel("Amplitude", fontsize = 14)

z, p, k = signal.tf2zpk(z_ricker[::10], [1])
zplane(z, p, k, 222)
plt.title("Diagrama de polos e zeros", fontsize = 18)

plt.subplot(223)
w, dft_z_ricker = dtft((t - lag)/dt, z_ricker)
plt.plot(w, np.abs(dft_z_ricker))
plt.xlim([0,fmax/180])
plt.xticks(np.linspace(0,fmax/180, 7), np.linspace(0,fmax, 7, dtype = int))
plt.title("Espectro de amplitudes", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 14)
plt.ylabel("Amplitude", fontsize = 14)

plt.subplot(224)
plt.plot(w, np.angle(dft_z_ricker, deg = True), "o")
plt.xlim([0, fmax/180])
plt.yticks(np.linspace(-180, 180, 11), np.linspace(-180, 180, 11, dtype = int))
plt.xticks(np.linspace(0,fmax/180, 7), np.linspace(0,fmax, 7, dtype = int))
plt.ylim([-180, 180])
plt.title("Espectro de fase", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 14)
plt.ylabel("Ângulo [°]", fontsize = 14)

plt.tight_layout()

#--------------------------------------------------------------

plt.figure(2,figsize=(18,8))

plt.subplot(221)
plt.plot(t, m_ricker)
plt.xlim([0, 2*lag])
plt.title("Ricker de fase mínima", fontsize = 18)
plt.xlabel("Tempo [s]", fontsize = 14)
plt.ylabel("Amplitude", fontsize = 14)

z, p, k = signal.tf2zpk(m_ricker[::10], [1])
zplane(z, p, k, 222)
plt.title("Diagrama de polos e zeros", fontsize = 18)

plt.subplot(223)
w, dft_m_ricker = dtft((t - lag)/dt, m_ricker)
plt.plot(w, np.abs(dft_m_ricker))
plt.xlim([0,fmax/180])
plt.xticks(np.linspace(0,fmax/180, 7), np.linspace(0,fmax, 7, dtype = int))
plt.title("Espectro de amplitudes", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 14)
plt.ylabel("Amplitude", fontsize = 14)

plt.subplot(224)
plt.plot(w, np.angle(dft_m_ricker, deg = True), "o")
plt.xlim([0, fmax/180])
plt.yticks(np.linspace(-180, 180, 11), np.linspace(-180, 180, 11, dtype = int))
plt.xticks(np.linspace(0,fmax/180, 7), np.linspace(0,fmax, 7, dtype = int))
plt.ylim([-180, 180])
plt.title("Espectro de fase", fontsize = 18)
plt.xlabel("Frequência [Hz]", fontsize = 14)
plt.ylabel("Ângulo [°]", fontsize = 14)

plt.tight_layout()
plt.show()
