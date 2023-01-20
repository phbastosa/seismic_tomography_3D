# include <chrono>
# include <iostream>

# include "../../src/wave/wave.hpp"
# include "../../src/essentials/model.hpp"
# include "../../src/essentials/utils.hpp"
# include "../../src/essentials/geometry.hpp"

int main (int argc, char**argv)
{
    auto ti = std::chrono::system_clock::now();

    auto utils = Utils();
    auto model = Model();
    auto geometry = Geometry();

    // Velocity model

    model.nb = 100;

    model.nx = 201;
    model.ny = 201;
    model.nz = 201;
    
    model.dx = 12.5f;
    model.dy = 12.5f;
    model.dz = 12.5f;

    model.initialize();

    // Wavefield setting

    float * U_pas = new float[model.nPointsB]();
    float * U_pre = new float[model.nPointsB]();
    float * U_fut = new float[model.nPointsB]();

    // Absorbing boundary condition - Perfectly Matched Layers (PML)

    float factor = 0.0045f;

    float * damp1D = new float[model.nb]();
    float * damp2D = new float[model.nb*model.nb]();
    float * damp3D = new float[model.nb*model.nb*model.nb]();

    pml_dampers(damp1D, damp2D, damp3D, factor, model.nb);

    float * V = new float[model.nPointsB];

    for (int i = 0; i < model.nPointsB; ++i) 
    {
        V[i] = 2000.0f;        
    }

    // Generate central shot

    geometry.geometryFolder = "outputs/";

    float x = static_cast<float>(model.nx-1)*model.dx / 2.0f;
    float y = static_cast<float>(model.ny-1)*model.dy / 2.0f;

    geometry.set_SW(x,y);     
    geometry.set_NW(x,y);    
    geometry.set_SE(x,y);    

    geometry.shots.n_xline = 1; 
    geometry.shots.n_yline = 1; 
    geometry.shots.elevation = static_cast<float>(model.nz-1)*model.dz / 2.0f;
    
    geometry.setGridGeometry(geometry.shots);

    // Generate carpet nodes

    geometry.set_SW(500.0f, 500.0f);     
    geometry.set_NW(500.0f, 2000.0f);    
    geometry.set_SE(2000.0f, 500.0f);    

    geometry.nodes.n_xline = 31; 
    geometry.nodes.n_yline = 31; 
    geometry.nodes.elevation = static_cast<float>(model.nz-1)*model.dz / 2.0f;
    
    geometry.setGridGeometry(geometry.nodes);

    // Boundaries compensation on integer geometry grid points

    int * rx = new int[geometry.nodes.all];
    int * ry = new int[geometry.nodes.all];
    int * rz = new int[geometry.nodes.all];

    for (int node = 0; node < geometry.nodes.all; node++)
    {
        rx[node] = static_cast<int>(geometry.nodes.x[node] / model.dx) + model.nb;
        ry[node] = static_cast<int>(geometry.nodes.y[node] / model.dy) + model.nb;
        rz[node] = static_cast<int>(geometry.nodes.z[node] / model.dz) + model.nb;
    }

    // Wavelet and time variables

    int nt = 1001;
    float dt = 0.001f;
    float fmax = 30.0f;
    float delay = 0.3f;

    std::string seismFolder = "outputs/";
    std::string waveletFolder = "outputs/";

    float * wavelet = rickerGeneration(delay, fmax, dt, nt);

    utils.writeBinaryFloat(waveletFolder + "ricker_" + std::to_string(static_cast<int>(fmax)) + "Hz_" + std::to_string(nt) + ".bin", wavelet, nt);

    float * seismogram = new float[nt * geometry.nodes.all]();
    
    for (int shotId = 0; shotId < geometry.shots.all; shotId++)
    {
        float sx = geometry.shots.x[shotId];
        float sy = geometry.shots.y[shotId];
        float sz = geometry.shots.z[shotId];
        
        int sIdx = static_cast<int>(sx / model.dx) + model.nb;
        int sIdy = static_cast<int>(sy / model.dy) + model.nb;    
        int sIdz = static_cast<int>(sz / model.dz) + model.nb;    

        int sId = sIdz + sIdx*model.nzz + sIdy*model.nzz*model.nxx;

        setWaveField(U_pas,U_pre,U_fut,model.nPointsB);

        # pragma acc enter data copyin(wavelet[0:nt])
        # pragma acc enter data copyin(damp1D[0:model.nb])
        # pragma acc enter data copyin(V[0:model.nPointsB])
        # pragma acc enter data copyin(U_pas[0:model.nPointsB])
        # pragma acc enter data copyin(U_pre[0:model.nPointsB])
        # pragma acc enter data copyin(U_fut[0:model.nPointsB])
        # pragma acc enter data copyin(damp2D[0:model.nb*model.nb])
        # pragma acc enter data copyin(damp3D[0:model.nb*model.nb*model.nb])
        # pragma acc enter data copyin(seismogram[0:nt*geometry.nodes.all])
        {    
            for (int timeStep = 0; timeStep < nt; timeStep++)
            {
                progressMessage(timeStep,nt,dt,shotId,sx,sy,sz);
                
                applyWavelet(U_pre,wavelet,timeStep,sId);

                pml_wavePropagation(V,U_pas,U_pre,U_fut,damp1D,damp2D,damp3D,model.nb,model.nxx,model.nyy,model.nzz,model.dx,model.dy,model.dz,dt);

                wavefieldUpdate(U_pas,U_pre,U_fut,model.nPointsB);

                buildSeismogram(U_fut,model.nxx,model.nyy,model.nzz,seismogram,timeStep,nt,rx,ry,rz,geometry.nodes.all);
            }
        }
        # pragma acc exit data delete(wavelet[0:nt])
        # pragma acc exit data delete(damp1D[0:model.nb])
        # pragma acc exit data delete(V[0:model.nPointsB])
        # pragma acc exit data delete(U_pas[0:model.nPointsB])
        # pragma acc exit data delete(U_pre[0:model.nPointsB])
        # pragma acc exit data copyout(U_fut[0:model.nPointsB])
        # pragma acc exit data delete(damp2D[0:model.nb*model.nb])
        # pragma acc exit data delete(damp3D[0:model.nb*model.nb*model.nb])
        # pragma acc exit data copyout(seismogram[0:nt*geometry.nodes.all])

        utils.writeBinaryFloat(seismFolder + "seismogram_" + std::to_string(nt) + "x" + std::to_string(geometry.nodes.all) + "_shot_" + std::to_string(shotId+1) + ".bin", seismogram, nt*geometry.nodes.all);
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "Run time: " << elapsed_seconds.count() << " s\n";

    return 0;
}
