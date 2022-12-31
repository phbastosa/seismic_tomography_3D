import numpy as np
import matplotlib.pyplot as plt

def readBinaryVolume(dim1,dim2,dim3,filename):
    data = np.fromfile(filename, dtype=np.float32, count=dim1*dim2*dim3)
    return np.reshape(data, [dim1,dim2,dim3], order='F')

nx = 31
ny = 31
nz = 31

mvSmooth = readBinaryVolume(nz, nx, ny, "outputs/mvSmooth.bin")
gnSmooth = readBinaryVolume(nz, nx, ny, "outputs/gnSmooth.bin")

sharpVol = np.zeros((nz,nx,ny))
sharpVol[15,15,15] = 100

images = [sharpVol[:,:,15],sharpVol[:,15,:],sharpVol[15,:,:],
          mvSmooth[:,:,15],mvSmooth[:,15,:],mvSmooth[15,:,:],
          gnSmooth[:,:,15],gnSmooth[:,15,:],gnSmooth[15,:,:]]

fig, axs = plt.subplots(3, 3, figsize = (12, 8))

ind = 0
for i in range(len(axs)):
    for j in range(len(axs[0])):
        axs[i,j].imshow(images[ind])
        ind += 1

        if i == 0:
            axs[i,j].set_title("sharp volume")
        elif i == 1:
            axs[i,j].set_title("moving average")
        elif i == 2:
            axs[i,j].set_title("gaussian filter")

        if j == 0:
            axs[i,j].set_xlabel("x")
            axs[i,j].set_ylabel("z")
        elif j == 1:
            axs[i,j].set_xlabel("y")
            axs[i,j].set_ylabel("z")
        elif j == 2:
            axs[i,j].set_xlabel("x")
            axs[i,j].set_ylabel("y")

plt.tight_layout()
plt.savefig("smoothing.png", dpi = 200)
plt.show(block=False)

# Plotting function interpolated

h = np.linspace(0, 1, 101)
z,x,y = np.meshgrid(h,h,h)

f = 1500 + 90*z + 30*x - 30*y

fig, axs = plt.subplots(1, 3, figsize = (10, 5))

axs[0].contour(f[:,:,50])
axs[0].imshow(f[:,:,50], cmap = "Greys")
axs[0].set_title("XZ plane")
axs[0].set_xlabel("X")
axs[0].set_ylabel("Z")

axs[1].contour(f[:,50,:])
axs[1].imshow(f[:,50,:], cmap = "Greys")
axs[1].set_title("YZ plane")
axs[1].set_xlabel("Y")
axs[1].set_ylabel("Z")

axs[2].contour(f[50,:,:])
axs[2].imshow(f[50,:,:], cmap = "Greys")
axs[2].set_title("XY plane")
axs[2].set_xlabel("X")
axs[2].set_ylabel("Y")

plt.tight_layout()
plt.savefig("triLinearInterpolation.png", dpi = 200)
plt.show(block = False)
