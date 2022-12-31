# include <string>
# include <chrono>
# include <iostream>

# include "../../src/essentials/utils.hpp"
# include "../../src/essentials/model.hpp"
# include "../../src/essentials/geometry.hpp"
# include "../../src/eikonal/eikonal.hpp"
# include "../../src/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomo = Tomography();

    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    tomo.setParameters(argv[1]);

    // Observed data generation --------------------------------------------------------------------

    tomo.V = tomo.expand(tomo.readBinaryFloat("outputs/trueModel_" + std::to_string(tomo.nz) + "x" + std::to_string(tomo.nx) + "x" + std::to_string(tomo.ny) + "_" + std::to_string((int) tomo.dx) + "m.bin", tomo.nPoints));
    
    std::cout<<"\nObserved data generation"<<std::endl;

    tomo.arrivalFolder = tomo.dobsPath;

    tomo.T = new float [tomo.nPointsB];

    for (tomo.shotId = 0; tomo.shotId < tomo.shots.all; tomo.shotId++)
    {
        tomo.eikonalComputing();
    }

    delete[] tomo.T;

    // Inversion part ----------------------------------------------------------------------

    tomo.V = tomo.expand(tomo.readBinaryFloat("outputs/initModel_" + std::to_string(tomo.nz) + "x" + std::to_string(tomo.nx) + "x" + std::to_string(tomo.ny) + "_" + std::to_string((int) tomo.dx) + "m.bin", tomo.nPoints));

    tomo.importDobs();
    tomo.setInitialModel();

    tomo.arrivalFolder = tomo.dcalPath;

    ti = std::chrono::system_clock::now();

    while (true)
    {
        tomo.forwardModeling();                
        
        tomo.importDcal();

        if (tomo.converged()) break;

        tomo.optimization();

        tomo.modelUpdate();
    }

    tomo.exportConvergency();

    tf = std::chrono::system_clock::now();

    elapsed_seconds = tf - ti;

    std::cout<<"\nTomography run time: "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}




