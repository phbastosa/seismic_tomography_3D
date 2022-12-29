import numpy as np

def buildModel3D(nx,ny,nz,property,z):
    model = np.zeros((nz,nx,ny), dtype=float)

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1,len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model

nx = 101
ny = 101
nz = 21

dh = 50.0

v = np.array([1500, 1700, 2000],dtype=float)
z = np.array([ 200, 800, nz*dh],dtype=float) // dh

initModel = buildModel3D(nx, ny, nz, v, z)

trueModel = initModel.copy()

trueModel[10:16,40:60,40:60] += 300

trueModel.flatten("F").astype("float32",order="F").tofile(f"outputs/trueModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")
initModel.flatten("F").astype("float32",order="F").tofile(f"outputs/initModel_{nz}x{nx}x{ny}_{dh:.0f}m.bin")