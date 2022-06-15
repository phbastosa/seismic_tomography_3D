# include <iostream>

# include "eikonal.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal();

    std::vector<std::string> modelNames {"refractiveModel_12x221x221_100m.bin", "refractiveModel_23x441x441_50m.bin", "refractiveModel_45x881x881_25m.bin"};

    std::vector<int> nx_all {221, 441, 881};
    std::vector<int> ny_all {221, 441, 881};
    std::vector<int> nz_all {12, 23, 45};
    
    std::vector<float> dh_all {100.0f, 50.0f, 25.0f};

    eikonal.nb = 5;

    // Set fixed circular geometry

    eikonal.nodes.xcenter = 11000.0f;
    eikonal.nodes.ycenter = 11000.0f;
    eikonal.nodes.offsets = {10000.0f}; // 10 km
    eikonal.nodes.elevation = 0.0f;
    eikonal.nodes.circle_spacing = 12.5f;

    eikonal.setCircularNodes();

    eikonal.nodesPath = "nodesPosition.txt";

    eikonal.exportPositions();

    eikonal.exportTimesVolume = false;
    eikonal.exportFirstArrivals = true;

    // Setting model 

    std::cout<<"\n\n";
    for (int n = 0; n < dh_all.size(); n++)
    {
        std::cout<<"Generating data with dh = "<<dh_all[n]<<"\n";

        eikonal.nx = nx_all[n];
        eikonal.ny = ny_all[n];    
        eikonal.nz = nz_all[n];

        eikonal.dx = dh_all[n];
        eikonal.dy = dh_all[n];
        eikonal.dz = dh_all[n];

        eikonal.initialize();

        float * model = new float[eikonal.nPoints]; 

        eikonal.readBinaryFloat(modelNames[n], model, eikonal.nPoints);

        eikonal.vp = eikonal.expandModel(model);

        delete[] model;

        // Setting extern shot points  

        eikonal.set_SW( 1000.0f,  1000.0f);    
        eikonal.set_NW( 1000.0f, 21000.0f);    
        eikonal.set_SE(21000.0f,  1000.0f);    

        eikonal.shots.n_xline = 2; 
        eikonal.shots.n_yline = 2; 
        eikonal.shots.elevation = 0.0f;

        eikonal.setGridShots();

        // Shots loop

        eikonal.arrivalFolder = "pod_extern_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 0;
        eikonal.forwardModeling();

        eikonal.arrivalFolder = "fim_extern_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 1;
        eikonal.forwardModeling();

        eikonal.arrivalFolder = "fsm_extern_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 2;
        eikonal.forwardModeling();

        // Generate central shot

        eikonal.set_SW(11000.0f, 11000.0f);    
        eikonal.set_NW(11000.0f, 11000.0f);    
        eikonal.set_SE(11000.0f, 11000.0f);    

        eikonal.shots.n_xline = 1; 
        eikonal.shots.n_yline = 1; 

        eikonal.setGridShots();
        
        eikonal.arrivalFolder = "pod_central_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 0;
        eikonal.forwardModeling();

        eikonal.arrivalFolder = "fim_central_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 1;
        eikonal.forwardModeling();

        if (n == 2)
        {
            eikonal.exportTimesVolume = true;
            eikonal.eikonalFolder = "central_";
        }

        eikonal.arrivalFolder = "fsm_central_"+std::to_string((int) dh_all[n])+"m_";
        eikonal.eikonalType = 2;
        eikonal.forwardModeling();
    }

    return 0;
}
