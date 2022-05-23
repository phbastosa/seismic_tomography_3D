# include "../src/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomography = Tomography3D(argv);

    tomography.importDobs();

    while (true)
    {
        tomography.fwdModeling();
        
        tomography.importDcal();

        tomography.makeGradient();

        if (tomography.converged()) break;

        


        tomography.iteration += 1;
    }


    return 0;
}
