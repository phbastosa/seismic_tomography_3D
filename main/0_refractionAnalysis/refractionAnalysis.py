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

v = np.array([2000,3000,4000])
z = np.array([800,300]) 

x = np.linspace(500, 9500, 721)

td, tr = analiticalRefraction(v,z,x)

plt.plot(td)
plt.plot(tr[0])
plt.plot(tr[1])

plt.gca().invert_yaxis()
plt.show()
