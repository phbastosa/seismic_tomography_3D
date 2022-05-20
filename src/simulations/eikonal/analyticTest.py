import numpy as np
import matplotlib.pyplot as plt

def readBinaryArray(n,filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

def analiticalRefraction2D(v,z,x):
    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

nx = 201
ny = 201
nz = 45
dh = 25.0
nr = 481

interface = 1000 / dh + 1

model = np.zeros((nz,nx,ny))

model[:int(interface),:,:] = 2000
model[int(interface):,:,:] = 3000

model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz}x{nx}x{ny}_{int(dh)}m.bin")

times = readBinaryVolume(nz,nx,ny,"eikonal_nz45_nx201_ny201_shot_1.bin")
arrivals = readBinaryArray(nr,"times_nr481_shot_1.bin")

shots = np.loadtxt("shots.txt", delimiter=",")
nodes = np.loadtxt("nodes.txt", delimiter=",")

x = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

v = np.array([2000, 3000])
z = np.array([1000])

td, t = analiticalRefraction2D(v,z,x)

offset = np.arange(nr) * 10 - 2400

plt.figure(1, figsize=(13,6))

plt.subplot(211)
plt.plot(offset, arrivals)
plt.plot(offset, t[0])
plt.plot(offset, td)

plt.subplot(212)
diffT = t[0] - arrivals

plt.plot(offset,diffT)

print(diffT[0], diffT[-1])
print(x[0] == x[-1])

plt.show()

