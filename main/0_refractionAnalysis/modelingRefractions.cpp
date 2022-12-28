# include <chrono>
# include <iostream>

# include "../../src/acoustic/acoustic.hpp"

int main(int argc, char **argv)
{
    auto acoustic = Acoustic();

    auto start = std::chrono::steady_clock::now();
    
    // Defining velocity model

    acoustic.nx = 1001;
    acoustic.ny = 81;
    acoustic.nz = 101;

    acoustic.dx = 12.5f;
    acoustic.dy = 12.5f;
    acoustic.dz = 12.5f;

    acoustic.nb = 50;
    acoustic.factor = 0.0045;

    acoustic.initialize();
    acoustic.dampingGenerator();

    acoustic.V = new float[acoustic.nPointsB]();

    for (int index = 0; index < acoustic.nPointsB; index++)
    {
        int k = (int) (index / (acoustic.nxx*acoustic.nzz));                  // y direction
        int j = (int) (index - k*acoustic.nxx*acoustic.nzz) / acoustic.nzz;   // x direction
        int i = (int) (index - j*acoustic.nzz - k*acoustic.nxx*acoustic.nzz); // z direction

        if (i < (int)(800/acoustic.dz) + acoustic.nb)
            acoustic.V[i + j*acoustic.nzz + k*acoustic.nxx*acoustic.nzz] = 2000.0f;
        else if ((i >= (int)(800/acoustic.dz) + acoustic.nb) && (i < (int)(1100/acoustic.dz) + acoustic.nb))    
            acoustic.V[i + j*acoustic.nzz + k*acoustic.nxx*acoustic.nzz] = 3000.0f; 
        else 
            acoustic.V[i + j*acoustic.nzz + k*acoustic.nxx*acoustic.nzz] = 4000.0f; 
    }

    // Defining time domain variables

    acoustic.nt = 4151;
    acoustic.dt = 0.001;
    acoustic.nsrc = 330;
    acoustic.fmax = 30.0f;

    // Defining geometry

    acoustic.shots.n_xline = 1;
    acoustic.shots.n_yline = 1;
    acoustic.shots.elevation = 0.0f;

    acoustic.set_SW(500.0f, 500.0f);
    acoustic.set_SE(500.0f, 500.0f);
    acoustic.set_NW(500.0f, 500.0f);

    acoustic.setGridGeometry(acoustic.shots);

    acoustic.nodes.n_xline = 1;
    acoustic.nodes.n_yline = 1;
    acoustic.nodes.elevation = 0.0f;

    acoustic.set_SE(9500.0f, 500.0f);

    acoustic.setGridGeometry(acoustic.nodes);

    acoustic.shotsPath = "shots.txt";
    acoustic.nodesPath = "nodes.txt";
    
    acoustic.exportPositions();

    // Defining zero phase source wavelet

    acoustic.source = acoustic.readBinaryFloat("../../inputs/wavelets/sourceZeroPhase_330_30Hz_1ms.bin", acoustic.nsrc);

    // Modeling acoustic wave

    acoustic.U_pas = new float[acoustic.nPointsB]();
    acoustic.U_pre = new float[acoustic.nPointsB]();
    acoustic.U_fut = new float[acoustic.nPointsB]();

    acoustic.seismogram = new float[acoustic.nt*acoustic.nodes.all]();

    acoustic.seisLabel = "zpSeism";

    acoustic.forwardModeling();

    // Defining min phase source wavelet

    acoustic.source = acoustic.readBinaryFloat("../../inputs/wavelets/sourceMinPhase_330_30Hz_1ms.bin", acoustic.nsrc);

    acoustic.seisLabel = "mpSeism";

    // acoustic.forwardModeling();

    auto end = std::chrono::steady_clock::now();

    std::cout << "\n\nElapsed time in seconds: "
              << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
              << " s"<< std::endl;

    return 0;
}

