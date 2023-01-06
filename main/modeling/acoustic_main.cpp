# include <chrono>
# include <iostream>

# include "../../src/essentials/wave.hpp"
# include "../../src/essentials/model.hpp"
# include "../../src/essentials/utils.hpp"
# include "../../src/essentials/geometry.hpp"

int main (int argc, char**argv)
{
    auto ti = std::chrono::system_clock::now();

    auto utils = Utils();
    auto model = Model();
    auto geometry = Geometry();

    std::string parameters = argv[1];

    // Velocity model

    model.nb = std::stoi(utils.catchParameter("nb", parameters));

    model.nx = std::stoi(utils.catchParameter("nx", parameters));
    model.ny = std::stoi(utils.catchParameter("ny", parameters));
    model.nz = std::stoi(utils.catchParameter("nz", parameters));
    
    model.dx = std::stof(utils.catchParameter("dx", parameters));
    model.dy = std::stof(utils.catchParameter("dy", parameters));
    model.dz = std::stof(utils.catchParameter("dz", parameters));

    model.initialize();

    std::string modelPath = utils.catchParameter("modelPath", parameters);

    float * V = model.expand(utils.readBinaryFloat(modelPath, model.nPoints));

    // Wavefield setting

    float * U_pas = new float[model.nPointsB]();
    float * U_pre = new float[model.nPointsB]();
    float * U_fut = new float[model.nPointsB]();

    // Absorbing boundary condition - Perfectly Matched Layers (PML)

    float factor = std::stof(utils.catchParameter("factor", parameters));

    float * damp1D = new float[model.nb]();
    float * damp2D = new float[model.nb*model.nb]();
    float * damp3D = new float[model.nb*model.nb*model.nb]();

    pml_dampers(damp1D, damp2D, damp3D, factor, model.nb);

    // Geometry

    geometry.reciprocity = utils.str2bool(utils.catchParameter("reciprocity", parameters));
    geometry.saveGeometry = utils.str2bool(utils.catchParameter("saveGeometry", parameters));

    geometry.shots.elevation = std::stof(utils.catchParameter("shotsElevation", parameters));
    geometry.nodes.elevation = std::stof(utils.catchParameter("nodesElevation", parameters));

    geometry.shots.n_xline = std::stoi(utils.catchParameter("xShotNumber", parameters));
    geometry.shots.n_yline = std::stoi(utils.catchParameter("yShotNumber", parameters));
    
    geometry.shotsPath = utils.catchParameter("shotsPath", parameters);
    geometry.nodesPath = utils.catchParameter("nodesPath", parameters);

    bool shotsTopography = utils.str2bool(utils.catchParameter("shotsTopography", parameters));
    bool nodesTopography = utils.str2bool(utils.catchParameter("nodesTopography", parameters));

    std::vector<std::string> splitted;

    // Setting grid shots
    
    splitted = utils.split(utils.catchParameter("shotSW", parameters),',');
    geometry.set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("shotNW", parameters),',');
    geometry.set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("shotSE", parameters),',');
    geometry.set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    geometry.setGridGeometry(geometry.shots);

    if (shotsTopography) 
    {
        std::string shotsTopographyPath = utils.catchParameter("shotsTopographyPath", parameters);
        geometry.shots.z = utils.readBinaryFloat(shotsTopographyPath, geometry.shots.all);
    }

    // Setting grid nodes

    geometry.nodes.n_xline = std::stoi(utils.catchParameter("xNodeNumber", parameters));
    geometry.nodes.n_yline = std::stoi(utils.catchParameter("yNodeNumber", parameters));

    splitted = utils.split(utils.catchParameter("nodeSW", parameters),',');
    geometry.set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("nodeNW", parameters),',');
    geometry.set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("nodeSE", parameters),',');
    geometry.set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    geometry.setGridGeometry(geometry.nodes);

    if (nodesTopography) 
    {
        std::string nodesTopographyPath = utils.catchParameter("nodesTopographyPath", parameters);
        geometry.nodes.z = utils.readBinaryFloat(nodesTopographyPath, geometry.nodes.all);
    }

    if (geometry.saveGeometry) 
        geometry.exportPositions();

    if (geometry.reciprocity)
        geometry.setReciprocity();

    // Boundaries compensation on integer geometry grid points

    int * rx = new int[geometry.nodes.all];
    int * ry = new int[geometry.nodes.all];
    int * rz = new int[geometry.nodes.all];

    for (int node = 0; node < geometry.nodes.all; node++)
    {
        rx[node] = (int)(geometry.nodes.x[node] / model.dx) + model.nb;
        ry[node] = (int)(geometry.nodes.y[node] / model.dy) + model.nb;
        rz[node] = (int)(geometry.nodes.z[node] / model.dz) + model.nb;
    }

    // Wavelet and time variables

    int nt = std::stoi(utils.catchParameter("nt", parameters));
    float dt = std::stof(utils.catchParameter("dt", parameters));
    float fmax = std::stof(utils.catchParameter("fmax", parameters));
    float delay = std::stof(utils.catchParameter("delay", parameters));

    std::string seismFolder = utils.catchParameter("seismFolder", parameters);
    std::string waveletFolder = utils.catchParameter("waveletFolder", parameters);

    float * wavelet = rickerGeneration(delay, fmax, dt, nt);

    utils.writeBinaryFloat(waveletFolder + "ricker_" + std::to_string((int)fmax) + "Hz_" + std::to_string(nt) + ".bin", wavelet, nt);

    float * seismogram = new float[nt * geometry.nodes.all]();

    // std::cout<<nt<<std::endl;    
    // std::cout<<model.nb<<std::endl;    
    // std::cout<<model.nb*model.nb<<std::endl;    
    // std::cout<<model.nb*model.nb*model.nb<<std::endl;    
    // std::cout<<model.nPointsB<<std::endl;    
    // std::cout<<model.nPointsB<<std::endl;    
    // std::cout<<model.nPointsB<<std::endl;    
    // std::cout<<model.nPointsB<<std::endl;    
    // std::cout<<nt*geometry.nodes.all<<std::endl;    
    
    for (int shotId = 0; shotId < 1; shotId++)
    {
        float sx = geometry.shots.x[shotId];
        float sy = geometry.shots.y[shotId];
        float sz = geometry.shots.z[shotId];
        
        int sIdx = (int)(sx / model.dx) + model.nb;
        int sIdy = (int)(sy / model.dy) + model.nb;    
        int sIdz = (int)(sz / model.dz) + model.nb;    

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

    std::cout << "Run time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
