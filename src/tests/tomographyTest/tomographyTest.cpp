# include <omp.h>
# include <string>
# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../essentials/model.hpp"
# include "../../essentials/geometry.hpp"
# include "../../eikonal/eikonal.hpp"
# include "../../tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomo = Tomography();

    tomo.setParameters(argv[1]);

    // Observed data generation --------------------------------------------------------------------

    // tomo.V = tomo.expand(tomo.readBinaryFloat("outputs/trueModel_" + std::to_string(tomo.nz) + "x" + std::to_string(tomo.nx) + "x" + std::to_string(tomo.ny) + "_" + std::to_string((int) tomo.dx) + "m.bin", tomo.nPoints));
    
    // std::cout<<"\nObserved data generation"<<std::endl;

    // tomo.arrivalFolder = tomo.dobsPath;

    // tomo.T = new float [tomo.nPointsB];

    // for (tomo.shotId = 0; tomo.shotId < tomo.shots.all; tomo.shotId++)
    // {
    //     tomo.eikonalComputing();
    // }

    // delete[] tomo.T;

    // Inversion part ----------------------------------------------------------------------

    tomo.V = tomo.expand(tomo.readBinaryFloat("outputs/initModel_" + std::to_string(tomo.nz) + "x" + std::to_string(tomo.nx) + "x" + std::to_string(tomo.ny) + "_" + std::to_string((int) tomo.dx) + "m.bin", tomo.nPoints));

    tomo.importDobs();
    tomo.setInitialModel();

    tomo.arrivalFolder = tomo.dcalPath;

    double t0 = omp_get_wtime();

    while (true)
    {
        tomo.forwardModeling();                
        
        tomo.importDcal();

        if (tomo.converged()) break;

        tomo.optimization();

        tomo.modelUpdate();
    }

    tomo.exportConvergency();

    std::cout<<"\nTomography run time: "<<omp_get_wtime() - t0<<" s."<<std::endl;

    return 0;
}




