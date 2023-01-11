# include <string>
# include <chrono>
# include <iostream>

# include "../../src/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomo = Tomography();
    
    auto ti = std::chrono::system_clock::now();
    
    tomo.parameters = std::string(argv[1]);

    tomo.setParameters();

    tomo.importDobs();
    tomo.setInitialModel();

    while (true)
    {
        tomo.forwardModeling();                
        
        tomo.importDcal();

        if (tomo.converged()) break;

        tomo.optimization();

        tomo.modelUpdate();
    }

    tomo.exportConvergency();

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds;

    elapsed_seconds = tf - ti;

    std::cout<<"\nTomography run time: "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}
