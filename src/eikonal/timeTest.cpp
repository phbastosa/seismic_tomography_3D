# include <omp.h>
# include "eikonal.hpp"

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

int main(int argc, char **argv)
{
    auto eikonal = Eikonal();

    eikonal.m3D.nx = 881;
    eikonal.m3D.ny = 881;    
    eikonal.m3D.nz = 45;

    eikonal.m3D.dx = 25;
    eikonal.m3D.dy = 25;
    eikonal.m3D.dz = 25;

    eikonal.m3D.nb = 2;

    eikonal.m3D.init();

    eikonal.m3D.vp = new float[eikonal.m3D.nPointsB];

    int interface = (int)(1000/eikonal.m3D.dz) + eikonal.m3D.nxx*eikonal.m3D.nzz + eikonal.m3D.nyy*eikonal.m3D.nxx*eikonal.m3D.nzz;

    for (int i = 0; i < eikonal.m3D.nPointsB; ++i) 
    {
        eikonal.m3D.vp[i] = 1500.0f;
        
        if (i >= interface)
            eikonal.m3D.vp[i] = 2000.0f;
    }

    // Generate central shot

    eikonal.g3D.SW.x = 11000.0f; eikonal.g3D.SW.y = 11000.0f;    
    eikonal.g3D.NW.x = 11000.0f; eikonal.g3D.NW.y = 11000.0f;    
    eikonal.g3D.SE.x = 11000.0f; eikonal.g3D.SE.y = 11000.0f;    

    eikonal.g3D.shots.nx = 1; 
    eikonal.g3D.shots.ny = 1; 
    eikonal.g3D.shots.elevation = 0.0f;

    eikonal.g3D.setGridShots();

    eikonal.allocateVolumes();    

    eikonal.exportTimesVolume = true;
    eikonal.exportFirstArrivals = false;

    eikonal.eikonalPath = "pod_";
    eikonal.arrivalsPath = "pod_";

    double t0 = omp_get_wtime();
    for (int shot = 0; shot < eikonal.g3D.shots.n; ++shot)
    {
        eikonal.shotId = shot;
        eikonal.podvin();
    }
    std::cout<<"\nPodvin & Lecomte (1991) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;

    eikonal.eikonalPath = "fim_";
    eikonal.arrivalsPath = "fim_";

    t0 = omp_get_wtime();
    for (int shot = 0; shot < eikonal.g3D.shots.n; ++shot)
    {
        eikonal.shotId = shot;
        eikonal.jeongFIM();
    }
    std::cout<<"\nJeong & Whitaker (2008) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;

    eikonal.eikonalPath = "fsm_";
    eikonal.arrivalsPath = "fsm_";

    t0 = omp_get_wtime();
    for (int shot = 0; shot < eikonal.g3D.shots.n; ++shot)
    {
        eikonal.shotId = shot;
        eikonal.nobleFSM();
    }
    std::cout<<"\nNoble, Gesret and Belayouni (2014) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;

    eikonal.deleteVolumes();

    return 0;
}
