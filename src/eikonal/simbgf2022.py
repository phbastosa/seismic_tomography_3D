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

nx = np.array([221, 441, 881], dtype = int)
ny = np.array([221, 441, 881], dtype = int)
nz = np.array([12, 23, 45], dtype = int) 

dh = np.array([100, 50, 25], dtype = float)

# # Creating models 
# for i in range(len(nx)):
#     interface = int(1000 / dh[i]) 
#     model = np.zeros((nz[i],nx[i],ny[i]))
#     model[:interface,:,:] = 1500
#     model[interface:,:,:] = 2000
#     model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz[i]}x{nx[i]}x{ny[i]}_{dh[i]:.0f}m.bin")




