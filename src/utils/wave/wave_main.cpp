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

    float * V = new float[nx*ny*nz]();

    read_binary_float("../../../inputs/models/initModel_101x701x501_10m.bin", V, nx*ny*nz);
    
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

    int n_shots = 11*16;
    int n_nodes = 140*100;

    std::vector<std::string> splitted;
    std::vector<std::string> elements;

    read_text_file("../../../inputs/geometry/xyz_nodes.txt", elements);

    float * sx = new float[n_shots];
    float * sy = new float[n_shots];
    float * sz = new float[n_shots];

    for (int i = 0; i < n_shots; i++)
    {
        splitted = split(elements[i],',');

        sx[i] = std::stof(splitted[0]);
        sy[i] = std::stof(splitted[1]);
        sz[i] = std::stof(splitted[2]);

        std::vector<std::string>().swap(splitted);
    }    

    std::vector<std::string>().swap(elements);

    read_text_file("../../../inputs/geometry/xyz_shots.txt", elements);

    float * rx = new float[n_nodes];
    float * ry = new float[n_nodes];
    float * rz = new float[n_nodes];

    for (int i = 0; i < n_nodes; i++)
    {
        splitted = split(elements[i],',');

        rx[i] = std::stof(splitted[0]);
        ry[i] = std::stof(splitted[1]);
        rz[i] = std::stof(splitted[2]);

        std::vector<std::string>().swap(splitted);
    }    

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

    float * wavelet = new float[nt];

    read_binary_float("ricker_min_phase_3001_1ms.bin", wavelet, nt);

    float * seismogram = new float[nt * n_nodes]();
    
    for (int shotId = 0; shotId < 1; shotId++)
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

        write_binary_float("../../../inputs/seismic_data/seismogram_" + std::to_string(nt) + "x" + std::to_string(n_nodes) + "_shot_" + std::to_string(shotId+1) + ".bin", seismogram, nt*n_nodes);
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "Run time: " << elapsed_seconds.count() << " s\n";

    return 0;
}
