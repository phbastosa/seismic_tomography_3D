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

    // Observed data generation

    tomo.eikonalType = 1;

    tomo.arrivalFolder = tomo.dobsPath;

    tomo.T = new float[tomo.nPointsB];
    
    std::cout<<"\nObserved data generation"<<std::endl;
    
    for (tomo.shotId = 0; tomo.shotId < tomo.shots.all; tomo.shotId++)
    {
        tomo.eikonalComputing();
    }

    delete[] tomo.T;

    return 0;
}




