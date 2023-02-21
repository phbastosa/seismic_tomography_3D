import numpy as np

nx = np.array([221, 441, 881], dtype = int)
ny = np.array([221, 441, 881], dtype = int)
nz = np.array([12, 23, 45], dtype = int) 

dh = np.array([100, 50, 25], dtype = float)

for i in range(len(dh)):
    
    interface = int(1000 / dh[i]) 
    
    model = np.zeros((nz[i],nx[i],ny[i]))
    
    model[:interface,:,:] = 1500
    model[interface:,:,:] = 2000

    model.flatten("F").astype("float32",order="F").tofile(f"outputs/refractiveModel_{nz[i]}x{nx[i]}x{ny[i]}_{dh[i]:.0f}m.bin")

    print(f"Model with dh = {dh[i]:.0f} m was written succesfully!")

print(" ")    
