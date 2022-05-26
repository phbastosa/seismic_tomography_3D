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

nx = np.array([441, 881, 1761], dtype = int)
ny = np.array([441, 881, 1761], dtype = int)
nz = np.array([23, 45, 89], dtype = int) 

dh = np.array([50, 25, 12.5], dtype = float)

# Creating models 
# for i in range(len(nx)):
#     interface = int(1000 / dh[i]) 
#     model = np.zeros((nz[i],nx[i],ny[i]))
#     model[:interface,:,:] = 1500
#     model[interface:,:,:] = 2000
    
#     if dh[i] == 12.5:
#         model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz[i]}x{nx[i]}x{ny[i]}_{dh[i]:.1f}m.bin")
#     else:
#         model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz[i]}x{nx[i]}x{ny[i]}_{dh[i]:.0f}m.bin")

n = 0
nr = 1256

tPod = readBinaryArray(nr,f"pod_{n+1}_{dh[n]:.0f}m.bin")
tFIM = readBinaryArray(nr,f"fim_{n+1}_{dh[n]:.0f}m.bin")

# plt.figure(2)
# plt.scatter(shots[0], shots[1], label = "Shots")
# plt.scatter(nodes[:,0], nodes[:,1], label = "Nodes")
# plt.title("Geometry acquisition")
# plt.legend()

# x = np.sqrt((shots[0] - nodes[:,0])**2 + (shots[1] - nodes[:,1])**2 + (shots[2] - nodes[:,2])**2)

# v = np.array([2000, 3000])
# z = np.array([1000])

# td, t = analiticalRefraction2D(v,z,x)

# shotId = np.arange(nr)

# plt.figure(1, figsize=(13,8))

# plt.subplot(211)
# plt.plot(shotId, tPod, label = "Podvin")
# plt.plot(shotId, tFIM, label = "FIM")
# plt.plot(shotId, t[0], label = "Analytic")
# plt.gca().invert_yaxis()
# plt.title("Travel times")
# plt.ylabel("Times [s]")
# plt.xlabel("Shot index")
# plt.legend()

# plt.subplot(212)
# plt.plot(shotId,np.abs(t[0] - tPod), label = "Podvin erros")
# plt.plot(shotId,np.abs(t[0] - tFIM), label = "FIM erros")
# plt.title("Absolute erros")
# plt.ylabel("abs(Ta - Tc) [s]")
# plt.xlabel("Shot index")

# plt.ylim([0,0.06])

# plt.tight_layout()
# plt.show()

