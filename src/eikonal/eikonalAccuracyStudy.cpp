# include "eikonal.hpp"

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal();

    std::vector<std::string> modelNames {"refractiveModel_12x221x221_100m.bin", "refractiveModel_23x441x441_50m.bin", "refractiveModel_45x881x881_25m.bin"};

    std::vector<int> nx_all {221, 441, 881};
    std::vector<int> ny_all {221, 441, 881};
    std::vector<int> nz_all {12, 23, 45};
    
    std::vector<float> dh_all {100.0f, 50.0f, 25.0f};

    eikonal.m3D.nb = 5;

    // Set fixed circular geometry

    eikonal.g3D.nodes.xc = 11000.0f;
    eikonal.g3D.nodes.yc = 11000.0f;
    eikonal.g3D.nodes.ds = 12.5f;
    eikonal.g3D.nodes.elevation = 0.0f;

    eikonal.g3D.nodes.offsets = {10000.0f}; // 10 km

    eikonal.g3D.setCircularNodes();

    eikonal.g3D.nodesPath = "nodesPosition.txt";

    eikonal.g3D.exportPositions();

    eikonal.exportTimesVolume = false;
    eikonal.exportFirstArrivals = true;

    // Setting model 

    std::cout<<"\n\n";
    for (int n = 0; n < dh_all.size(); n++)
    {
        std::cout<<"Generating data with dh = "<<dh_all[n]<<"\n";

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

        eikonal.g3D.shots.nx = 2; 
        eikonal.g3D.shots.ny = 2; 
        eikonal.g3D.shots.elevation = 0.0f;

        eikonal.g3D.setGridShots();

        // Shots loop

        eikonal.allocateVolumes();    

        eikonal.arrivalsPath = "pod_extern_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.podvin();
        }

        eikonal.arrivalsPath = "fim_extern_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.jeongFIM();
        }

        eikonal.arrivalsPath = "fsm_extern_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.jeongFIM();
        }

        // Generate central shot

        eikonal.g3D.SW.x = 11000.0f; eikonal.g3D.SW.y = 11000.0f;    
        eikonal.g3D.NW.x = 11000.0f; eikonal.g3D.NW.y = 11000.0f;    
        eikonal.g3D.SE.x = 11000.0f; eikonal.g3D.SE.y = 11000.0f;    

        eikonal.g3D.shots.nx = 1; 
        eikonal.g3D.shots.ny = 1; 

        eikonal.g3D.setGridShots();
        
        eikonal.arrivalsPath = "pod_central_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.podvin();
        }

        eikonal.arrivalsPath = "fim_central_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.jeongFIM();
        }

        if (n == 2)
        {
            eikonal.exportTimesVolume = true;
            eikonal.eikonalPath = "central_";
        }

        eikonal.arrivalsPath = "fsm_central_"+std::to_string((int) dh_all[n])+"m_";

        for (int shot = 0; shot < eikonal.g3D.shots.n; shot++)
        {
            eikonal.shotId = shot;
            eikonal.nobleFSM();
        }

        eikonal.deleteVolumes();
    }

    return 0;
}
