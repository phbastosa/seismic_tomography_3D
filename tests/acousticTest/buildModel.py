import numpy as np
import matplotlib.pyplot as plt

def buildModel3D(nx, ny, nz, property, z):
    model = np.zeros((nz,nx,ny), dtype=float)

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1,len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model

nx = 101
ny = 101
nz = 51

dh = 12.5

v = np.linspace(1500, 2000, 7)
z = np.linspace(200, nz*dh, 7)//dh

model = buildModel3D(nx, ny, nz, v, z)

model.flatten('F').astype('float32', order = 'F').tofile(f'outputs/model_{nz}x{nx}x{ny}.bin')
