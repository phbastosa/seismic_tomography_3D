# include <omp.h>
# include <iostream>

# include "../src/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomography = Tomography(argv);

    double start = omp_get_wtime();

    tomography.importDobs();

    tomography.setInitialModel();

    while (true)
    {
        tomography.forwardModeling();
        
        tomography.importDcal();

        tomography.makeGradient();

        if (tomography.converged()) break;

        tomography.cgls_Berriman();

        tomography.modelUpdate();

        tomography.modelSmoothing();
    }

    tomography.exportConvergency();

    std::cout<<"\nRun time "<<omp_get_wtime() - start<<" s."<<std::endl;

    return 0;
}
