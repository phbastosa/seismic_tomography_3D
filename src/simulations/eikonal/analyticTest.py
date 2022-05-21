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

nx = 441
ny = 441
nz = 23
dh = 50.0

ns = 1
nr = 1256

interface = 1000 / dh 

model = np.zeros((nz,nx,ny))

model[:int(interface),:,:] = 2000
model[int(interface):,:,:] = 3000

model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz}x{nx}x{ny}_{int(dh)}m.bin")

tPod = readBinaryArray(nr,f"podvin_times_nr{nr}_shot_{ns}.bin")
tFIM = readBinaryArray(nr,f"fim_times_nr{nr}_shot_{ns}.bin")

shots = np.loadtxt(f"shots_n{ns}.txt", delimiter=",")
nodes = np.loadtxt(f"nodes_n{nr}.txt", delimiter=",")

plt.figure(2)
plt.scatter(shots[0], shots[1], label = "Shots")
plt.scatter(nodes[:,0], nodes[:,1], label = "Nodes")
plt.title("Geometry acquisition")
plt.legend()

x = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

v = np.array([2000, 3000])
z = np.array([1000])

td, t = analiticalRefraction2D(v,z,x)

shotId = np.arange(nr)

plt.figure(1, figsize=(13,8))

plt.subplot(211)
plt.plot(shotId, tPod, label = "Podvin")
plt.plot(shotId, tFIM, label = "FIM")
plt.plot(shotId, t[0], label = "Analytic")
plt.gca().invert_yaxis()
plt.title("Travel times")
plt.ylabel("Times [s]")
plt.xlabel("Shot index")
plt.legend()

plt.subplot(212)
plt.plot(shotId,np.abs(t[0] - tPod), label = "Podvin erros")
plt.plot(shotId,np.abs(t[0] - tFIM), label = "FIM erros")
plt.title("Absolute erros")
plt.ylabel("abs(Ta - Tc) [s]")
plt.xlabel("Shot index")

plt.ylim([0,0.06])

plt.tight_layout()
plt.show()

