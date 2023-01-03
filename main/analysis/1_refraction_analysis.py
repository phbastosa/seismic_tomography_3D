import numpy as np
import matplotlib.pyplot as plt

def analiticalRefraction(v,z,x):

    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

nx = 801
ny = 161
nz = 201

dh = 12.5

v = np.array([2000,3500])
z = np.array([1200]) 

refractiveModel = np.zeros((nz,nx,ny))

refractiveModel[:int(z[0]/dh),:,:] = v[0]
refractiveModel[int(z[0]/dh):,:,:] = v[1]

refractiveModel.flatten("F").astype("float32", order = "F").tofile(f"../../inputs/models/refractiveModel_{nz}x{nx}x{ny}_{dh:.1f}m.bin")

x = np.linspace(1000, 9500, 321)

td, tr = analiticalRefraction(v,z,x)

plt.plot(td)
plt.plot(tr[0])

plt.gca().invert_yaxis()
plt.show()
