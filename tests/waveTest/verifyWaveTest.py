import numpy as np
import matplotlib.pyplot as plt

def readBinaryMatrix(dim1,dim2,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2)
    return np.reshape(data, [dim1,dim2], order='F')

def readBinaryArray(dim,filename):
    return np.fromfile(filename, dtype=np.float32, count=dim)

nr = 961
nt = 1001

waveletPath = "outputs/ricker_30Hz_1001.bin"
seismogramPath = "outputs/seismogram_1001x961_shot_1.bin"

wavelet = readBinaryArray(nt, waveletPath)
seismogram = readBinaryMatrix(nt, nr, seismogramPath)

scale = 0.5 * np.std(seismogram)

plt.figure(1, figsize = (18, 6))
plt.imshow(seismogram, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)

plt.tight_layout()
plt.show()



plt.figure(2, figsize = (10,5))

plt.plot(wavelet)

plt.show()