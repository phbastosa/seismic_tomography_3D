# include <chrono>
# include <iostream>

# include "../../src/eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    ti = std::chrono::system_clock::now();

    auto eikonal = Eikonal();

    eikonal.setEikonalParameters(argv[1]);

    for (eikonal.shotId = 0; eikonal.shotId < eikonal.shots.all; eikonal.shotId++)
    {
        eikonal.eikonalComputing();
    }

    eikonal.writeIllumination();

    tf = std::chrono::system_clock::now();

    elapsed_seconds = tf - ti;
    std::cout<<"\nRun time = "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}