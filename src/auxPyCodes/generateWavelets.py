import numpy as np
import matplotlib.pyplot as plt

from scipy import signal, interpolate
from matplotlib import patches

def rickerGenerator(fmax,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fmax/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1.0-2.0*np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0)
        aux2 = np.exp(-np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0))
        ricker[i + s] = aux1 * aux2 

    return ricker

def intRickerGenerator(fmax, nsrc, dt):    
    
    ricker = rickerGenerator(fmax, nsrc, dt)    
    
    intRicker = np.zeros(len(ricker))
     
    for i in range(len(ricker)):
        intRicker[i] += np.sum(ricker[:i+1]) 

    return intRicker

def zplane(z, p, subplot, title):

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
    
    ax.set_title(title)

    # set the ticks
    r = 1.5; plt.axis('scaled'); plt.axis([-r, r, -r, r])
    ticks = [-1, -.5, .5, 1]; plt.xticks(ticks); plt.yticks(ticks)

def minPhasing(nsrc, dt, fmax, sg):

    j = complex(0,1)

    dt_ref = dt * 10 
    nsrc_ref = int(nsrc / 10) + 1 

    ricker = sg(fmax,nsrc_ref,dt_ref)

    a = [1]
    b = ricker

    z, p, k = signal.tf2zpk(b, a)

    dis = np.abs(z)
    ang = np.angle(z)

    pos = np.where(dis > 1)

    dis[pos] = 1/dis[pos]

    z_in = np.zeros(len(z),dtype=complex)

    z_in[pos] = dis[pos]*np.cos(ang[pos]) + j*dis[pos]*np.sin(ang[pos]) 

    b_in, a_in = signal.zpk2tf(z_in, p, k)

    system = signal.TransferFunction(b_in, a_in, dt=dt_ref)

    _ , reconst = signal.dfreqresp(system,n=int(nsrc_ref/2))

    rickerReconst = np.real(np.fft.ifft(reconst,n=nsrc_ref))
    t_ref = np.arange(nsrc_ref) * dt_ref
    t = np.arange(nsrc) * dt

    interp = interpolate.interp1d(t_ref, rickerReconst, kind = 'cubic')

    return interp(t), z, z_in, p

nsrc = 201
dt = 0.001
fmax = 30

ilag = int(nsrc / 2)

sourceZeroPhase = intRickerGenerator(fmax,nsrc,dt)
sourceMinPhase, z, z_in, p = minPhasing(nsrc, dt, fmax, intRickerGenerator)

# normalization

sourceMinPhase *= 1.0 / np.max(np.abs(sourceMinPhase)) 
sourceZeroPhase *= 1.0 / np.max(np.abs(sourceZeroPhase)) 

sourceZeroPhase[0] = 0.0

sourceMinPhase.astype("float32", order = "F").tofile(f"inputs/wavelets/sourceMinPhase_{nsrc}_{dt*1e4:.0f}ms.bin")
sourceZeroPhase.astype("float32", order = "F").tofile(f"inputs/wavelets/sourceZeroPhase_{nsrc}_{dt*1e4:.0f}ms.bin")

tmp = np.arange(nsrc) * dt           # t min phase
tzp = (np.arange(nsrc) - ilag) * dt  # t zero phase

plt.figure(1, figsize = (12,8))

plt.subplot(321)
plt.plot(tzp, sourceZeroPhase)
plt.ylim(np.min(sourceZeroPhase)-0.1,np.max(sourceZeroPhase)+0.1)
plt.title("Zero phase wavelet - first gaussian derivative")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude normalized")

plt.subplot(322)

samp = int(6 * nsrc)
huge = np.zeros(samp)

huge[int(samp/2):int(samp/2 + nsrc)] = sourceZeroPhase

szp_fft = np.fft.fft(huge)
szp_fft *= 1 / np.max(np.abs(szp_fft))

freqs = np.fft.fftfreq(samp, dt)
mask = freqs >= 0

plt.plot(freqs[mask], np.abs(szp_fft[mask]))
plt.xlim([0,30])

plt.title("Zero phase wavelet spectra")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude normalized")

plt.subplot(323)
plt.plot(tmp, sourceMinPhase)
plt.ylim(np.min(sourceMinPhase)-0.1,np.max(sourceMinPhase)+0.2)
plt.title("Min phase wavelet - first gaussian derivative")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude normalized")

plt.subplot(324)

samp = int(6 * nsrc)
huge = np.zeros(samp)

huge[int(samp/2):int(samp/2 + nsrc)] = sourceMinPhase

smp_fft = np.fft.fft(huge)
smp_fft *= 1 / np.max(np.abs(smp_fft))

freqs = np.fft.fftfreq(samp, dt)
mask = freqs >= 0

plt.plot(freqs[mask], np.abs(smp_fft[mask]))
plt.xlim([0,30])

plt.title("Min phase wavelet spectra")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude normalized")

zplane(z, p, 325, "Zero phase wavelet diagram")

zplane(z_in, p, 326, "Min phase wavelet diagram")

plt.tight_layout()
plt.show()







# plt.figure(1,figsize=(18,8))
# xmask = np.arange(pn[0],pn[-1]+1,dtype=int)

# plt.subplot(221)
# plt.stem(pn,ricker)
# plt.xticks(xmask,xmask)
# plt.title("Sinal discreto")
# plt.xlabel("Amostras")
# plt.ylabel("Amplitude")

# zplane(z, p, 222)
# plt.title("Diagrama de polos e zeros")

# plt.subplot(223)
# plt.stem(freqs[mask],np.abs(fftRicker[mask]))
# plt.xlim([0,fmax])
# plt.title("Espectro de amplitudes")
# plt.xlabel("Frequência [Hz]")
# plt.ylabel("Amplitude")

# plt.subplot(224)
# plt.plot(freqs[mask],np.angle(fftRicker[mask],deg=True))
# plt.xlim([0,fmax])
# plt.ylim([-180,180])
# plt.title("Espectro de fase")
# plt.xlabel("Frequência [Hz]")
# plt.ylabel("Ângulo")

# plt.tight_layout()

# #--------------------------------------------------------------

# plt.figure(2,figsize=(18,8))
# plt.subplot(221)
# plt.stem(pn,rickerReconst)
# plt.xticks(xmask,xmask)
# plt.title("Sinal discreto")
# plt.xlabel("Amostras")
# plt.ylabel("Amplitude")

# zplane(z_in, p, 222)
# plt.title("Diagrama de polos e zeros")

# plt.subplot(223)
# plt.stem(freqs[mask],np.abs(fftRickerReconst[mask]))
# plt.xlim([0,fmax])
# plt.title("Espectro de amplitudes")
# plt.xlabel("Frequência [Hz]")
# plt.ylabel("Amplitude")

# plt.subplot(224)
# plt.plot(freqs[mask],np.angle(fftRickerReconst[mask],deg=True))
# plt.xlim([0,fmax])
# plt.ylim([-180,180])
# plt.title("Espectro de fase")
# plt.xlabel("Frequência [Hz]")
# plt.ylabel("Ângulo")

# plt.tight_layout()
 
# plt.show()
