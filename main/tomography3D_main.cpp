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

        if (tomography.converged()) break;

        tomography.optimization();

        tomography.modelUpdate();
    }

    tomography.exportConvergency();

    std::cout<<"\nRun time "<<omp_get_wtime() - start<<" s."<<std::endl;

    return 0;
}
