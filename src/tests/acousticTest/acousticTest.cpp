# include <chrono>
# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../acoustic/acoustic.hpp"

int main(int argc, char**argv)
{
    auto acoustic = Acoustic();

    auto start = std::chrono::steady_clock::now();
 
    // Testing source 

    acoustic.nt = 1001;
    acoustic.dt = 0.001f;
    acoustic.fcut = 30.0f;
    acoustic.tlag = 0.15f;

    acoustic.sourceGenerator();

    acoustic.writeBinaryFloat("outputs/source_Nt"+ std::to_string(acoustic.nsrc) + ".bin", acoustic.source, acoustic.nsrc);

    // Testing Cerjan factor attenuation

    acoustic.nb = 50;
    acoustic.factor = 0.0045;

    acoustic.dampingGenerator();

    acoustic.writeBinaryFloat("outputs/damp1D_" + std::to_string(acoustic.nb) + "_samples.bin", acoustic.damp1D, acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp2D_" + std::to_string(acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp2D, acoustic.nb*acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp3D_" + std::to_string(acoustic.nb*acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp3D, acoustic.nb*acoustic.nb*acoustic.nb);

    // Defining homogeneous model

    acoustic.nx = 101;
    acoustic.ny = 101;    
    acoustic.nz = 101;

    acoustic.dx = 12.5f;
    acoustic.dy = 12.5f;
    acoustic.dz = 12.5f;

    acoustic.initialize();

    acoustic.V = new float[acoustic.nPointsB];

    for (int index = 0; index < acoustic.nPointsB; index++) 
        acoustic.V[index] = 2500.0f;
        
    // Applying geometry acquisition 

    // Shots
    acoustic.set_SW(625.0f, 625.0f);     
    acoustic.set_NW(625.0f, 625.0f);    
    acoustic.set_SE(625.0f, 625.0f);    

    acoustic.shots.n_xline = 1; 
    acoustic.shots.n_yline = 1; 
    acoustic.shots.elevation = 625.0f;
    
    acoustic.setGridGeometry(acoustic.shots);

    // Receivers
    acoustic.set_SW(500.0f, 500.0f);     
    acoustic.set_NW(500.0f, 4500.0f);    
    acoustic.set_SE(4500.0f, 500.0f);    

    acoustic.nodes.n_xline = 21; 
    acoustic.nodes.n_yline = 21; 
    acoustic.nodes.elevation = 1250.0f;
    
    acoustic.setGridGeometry(acoustic.nodes);

    acoustic.shotsPath = "outputs/shotsPosition.txt";
    acoustic.nodesPath = "outputs/nodesPosition.txt";

    acoustic.exportPositions();

    // Preparing shots loop

    acoustic.forwardModeling();

    auto end = std::chrono::steady_clock::now();

    std::cout << "\n\nElapsed time in seconds: "
              << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
              << " s"<< std::endl;

    return 0;
}

