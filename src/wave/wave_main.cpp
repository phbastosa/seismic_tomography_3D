# include <chrono>
# include <iostream>

# include "acoustic.hpp"

int main (int argc, char**argv)
{
    auto ti = std::chrono::system_clock::now();

    // Velocity model

    int nb = 50;

    int nx = 701;
    int ny = 501;
    int nz = 101;

    float dx = 10.0f;
    float dy = 10.0f;
    float dz = 10.0f;    
    
    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    int nPoints = nxx*nyy*nzz;

    float * V = readBinaryFloat("model/vp_101x701x501_10m.bin", nx*ny*nz);
    
    V = expand(V, nx, ny, nz, nb);

    // Wavefield setting

    float * U_pas = new float[nPoints]();
    float * U_pre = new float[nPoints]();
    float * U_fut = new float[nPoints]();

    // Absorbing boundary condition - Perfectly Matched Layers (PML)

    float factor = 0.0045;

    float * damp1D = new float[nb]();
    float * damp2D = new float[nb*nb]();
    float * damp3D = new float[nb*nb*nb]();

    dampers(damp1D, damp2D, damp3D, factor, nb);

    // Geometry

    int n_shots = 441;
    int n_nodes = 10000; 

    float * sx = readBinaryFloat("geometry/x_nodes_441_positions.bin", n_shots);
    float * sy = readBinaryFloat("geometry/y_nodes_441_positions.bin", n_shots);
    float * sz = readBinaryFloat("geometry/z_nodes_441_positions.bin", n_shots);
    
    float * rx = readBinaryFloat("geometry/x_shots_10000_positions.bin", n_nodes);
    float * ry = readBinaryFloat("geometry/y_shots_10000_positions.bin", n_nodes);
    float * rz = readBinaryFloat("geometry/z_shots_10000_positions.bin", n_nodes);

    // Boundaries compensation on integer geometry grid points

    int * rIdx = new int[n_nodes];
    int * rIdy = new int[n_nodes];
    int * rIdz = new int[n_nodes];

    for (int node = 0; node < n_nodes; node++)
    {
        rIdx[node] = (int)(rx[node] / dx) + nb;
        rIdy[node] = (int)(ry[node] / dy) + nb;
        rIdz[node] = (int)(rz[node] / dz) + nb;
    }

    // Wavelet and time variables

    int nt = 3001;
    float dt = 0.001f;
    // float fmax = 30;
    // float delay = 0.1f;

    // float * wavelet = rickerGeneration(delay, fmax, dt, nt);

    float * wavelet = readBinaryFloat("src/ricker_min_phase_3001_1ms.bin", nt);

    float * seismogram = new float[nt * n_nodes]();
    
    for (int shotId = 250; shotId < 251; shotId++)
    {
        float source_x = sx[shotId];
        float source_y = sy[shotId];
        float source_z = sz[shotId];

        int sIdx = (int)(source_x / dx) + nb;
        int sIdy = (int)(source_y / dy) + nb;    
        int sIdz = (int)(source_z / dz) + nb;    

        int sId = sIdz + sIdx*nzz + sIdy*nzz*nxx;

        setWaveField(U_pas,U_pre,U_fut,nPoints);

        # pragma acc enter data copyin(wavelet[0:nt])
        # pragma acc enter data copyin(damp1D[0:nb])
        # pragma acc enter data copyin(V[0:nPoints])
        # pragma acc enter data copyin(U_pas[0:nPoints])
        # pragma acc enter data copyin(U_pre[0:nPoints])
        # pragma acc enter data copyin(U_fut[0:nPoints])
        # pragma acc enter data copyin(damp2D[0:nb*nb])
        # pragma acc enter data copyin(damp3D[0:nb*nb*nb])
        # pragma acc enter data copyin(seismogram[0:nt*n_nodes])
        {    
            for (int timeStep = 0; timeStep < nt; timeStep++)
            {
                progressMessage(timeStep,nt,dt,shotId,source_x,source_y,source_z);
                
                applyWavelet(U_pre,wavelet,timeStep,sId);

                wavePropagation(V,U_pas,U_pre,U_fut,damp1D,damp2D,damp3D,nb,nxx,nyy,nzz,dx,dy,dz,dt);

                wavefieldUpdate(U_pas,U_pre,U_fut,nPoints);

                buildSeismogram(U_fut,nxx,nyy,nzz,seismogram,timeStep,nt,rIdx,rIdy,rIdz,n_nodes);
            }
        }
        # pragma acc exit data delete(wavelet[0:nt])
        # pragma acc exit data delete(damp1D[0:nb])
        # pragma acc exit data delete(V[0:nPoints])
        # pragma acc exit data delete(U_pas[0:nPoints])
        # pragma acc exit data delete(U_pre[0:nPoints])
        # pragma acc exit data copyout(U_fut[0:nPoints])
        # pragma acc exit data delete(damp2D[0:nb*nb])
        # pragma acc exit data delete(damp3D[0:nb*nb*nb])
        # pragma acc exit data copyout(seismogram[0:nt*n_nodes])

        writeBinaryFloat("seismogram_" + std::to_string(nt) + "x" + std::to_string(n_nodes) + "_shot_" + std::to_string(shotId+1) + ".bin", seismogram, nt*n_nodes);
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "Run time: " << elapsed_seconds.count() << " s\n";

    return 0;
}
