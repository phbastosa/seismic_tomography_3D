# include "eikonal.hpp"

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal();

    std::vector<std::string> modelNames {"refractiveModel_23x441x441_50m.bin", "refractiveModel_45x881x881_25m.bin", "refractiveModel_89x1761x1761_12.5m.bin"};

    std::vector<int> nx_all {441, 881, 1761};
    std::vector<int> ny_all {441, 881, 1761};
    std::vector<int> nz_all {23, 45, 89};
    
    std::vector<float> dh_all {50.0f, 25.0f, 12.5f};

    eikonal.m3D.nb = 2;

    // Set fixed circular geometry

    eikonal.g3D.circles.xc = 11000.0f;
    eikonal.g3D.circles.yc = 11000.0f;
    eikonal.g3D.circles.ds = 50.0f;

    eikonal.g3D.circles.offsets = {10000.0f}; // 10 km

    eikonal.g3D.rElev = 0.0f;

    eikonal.g3D.setCircularNodes();

    eikonal.exportFirstArrivals = true;

    // Setting model 

    int n = 2;

    eikonal.m3D.nx = nx_all[n];
    eikonal.m3D.ny = ny_all[n];    
    eikonal.m3D.nz = nz_all[n];

    eikonal.m3D.dx = dh_all[n];
    eikonal.m3D.dy = dh_all[n];
    eikonal.m3D.dz = dh_all[n];

    eikonal.m3D.init();

    eikonal.m3D.vp = new float[eikonal.m3D.nPointsB];

    eikonal.m3D.vpPath = modelNames[n];

    eikonal.m3D.readAndExpandVP();

    // Setting extern shot points  

    eikonal.g3D.SW.x = 1000.0f; eikonal.g3D.SW.y = 1000.0f;    
    eikonal.g3D.NW.x = 1000.0f; eikonal.g3D.NW.y = 21000.0f;    
    eikonal.g3D.SE.x = 21000.0f; eikonal.g3D.SE.y = 1000.0f;    

    eikonal.g3D.nsx = 2; 
    eikonal.g3D.nsy = 2; 
    eikonal.g3D.sElev = 0.0f;

    eikonal.g3D.setGridShots();

    // Shots loop

    eikonal.allocateVolumes();    

    eikonal.arrivalsPath = "podvin_extern_";

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.podvin();
    }

    eikonal.arrivalsPath = "fim_extern_";

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.jeongFIM();
    }

    // Generate central shot

    eikonal.g3D.SW.x = 11000.0f; eikonal.g3D.SW.y = 11000.0f;    
    eikonal.g3D.NW.x = 11000.0f; eikonal.g3D.NW.y = 11000.0f;    
    eikonal.g3D.SE.x = 11000.0f; eikonal.g3D.SE.y = 11000.0f;    

    eikonal.g3D.nsx = 1; 
    eikonal.g3D.nsy = 1; 

    eikonal.g3D.setGridShots();

    eikonal.arrivalsPath = "podvin_central_";

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.podvin();
    }

    eikonal.arrivalsPath = "fim_central_";

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.jeongFIM();
    }

    eikonal.deleteVolumes();

    return 0;
}
