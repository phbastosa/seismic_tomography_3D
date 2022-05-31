import sys
import numpy as np

def buildModel3D(nx,ny,nz,property,z):
    model = np.zeros((nz,nx,ny))

    model[:int(z[0]),:,:] = np.ones((int(z[0]),nx,ny)) * property[0]

    for depth in range(1,len(z)):
        
        layer = slice(int(z[depth-1]),int(z[depth]))
    
        model[layer,:,:] = np.ones((int(z[depth] - z[depth-1]),nx,ny)) * property[depth]

    return model    

nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])

dh = float(sys.argv[4])

vp = np.array(sys.argv[5][1:-1].split(sep=","),dtype=float)

z = np.array(sys.argv[6][1:-1].split(sep=","),dtype=float) / dh

vpModel = buildModel3D(nx,ny,nz,vp,z)

vpModel.astype("float32",order="C").tofile(sys.argv[7])