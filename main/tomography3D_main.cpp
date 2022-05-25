# include "../src/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomography = Tomography(argv);

    tomography.importDobs();

    tomography.setInitialModel();

    while (true)
    {
        tomography.fwdModeling();
        
        tomography.importDcal();

        tomography.makeGradient();

        if (tomography.converged()) break;

        tomography.cgls_Berriman();

        // tomography.cgls_zoTikhonov();
        // tomography.cgls_foTikhonov();
        // tomography.cgls_soTikhonov();

        tomography.modelUpdate();

        tomography.modelSmoothing();
    }

    tomography.exportConvergency();

    tomography.deleteVolumes();

    return 0;
}
