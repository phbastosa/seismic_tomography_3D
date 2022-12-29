# include <iostream>

# include "../../eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal();

    std::vector<std::string> modelNames {"outputs/refractiveModel_12x221x221_100m.bin", 
                                         "outputs/refractiveModel_23x441x441_50m.bin", 
                                         "outputs/refractiveModel_45x881x881_25m.bin"};

    std::vector<int> nx_all {221, 441, 881};
    std::vector<int> ny_all {221, 441, 881};
    std::vector<int> nz_all {12, 23, 45};
    
    std::vector<float> dh_all {100.0f, 50.0f, 25.0f};

    eikonal.nb = 5;

    // Set fixed circular geometry

    eikonal.nodes.elevation = 0.0f;

    eikonal.nodes.xcenter = 11000.0f;
    eikonal.nodes.ycenter = 11000.0f;
    eikonal.nodes.offsets = {10000.0f};   // 10 km
    eikonal.nodes.circle_spacing = 12.5f;

    eikonal.setCircularGeometry(eikonal.nodes);

    eikonal.nodesPath = "outputs/nodesPosition.txt";

    eikonal.exportPositions();

    eikonal.exportTimesVolume = false;
    eikonal.exportFirstArrivals = true;

    // Setting model 

    std::cout<<"\n\n";
    for (int n = 0; n < dh_all.size(); n++)
    {
        std::cout<<"Generating data with dh = "<<dh_all[n]<<" m\n";

        eikonal.nx = nx_all[n];
        eikonal.ny = ny_all[n];    
        eikonal.nz = nz_all[n];

        eikonal.dx = dh_all[n];
        eikonal.dy = dh_all[n];
        eikonal.dz = dh_all[n];

        eikonal.initialize();
        
        eikonal.V = eikonal.expand(eikonal.readBinaryFloat(modelNames[n], eikonal.nPoints));

        // Setting extern shot points  

        eikonal.shots.elevation = 0.0f;

        eikonal.set_SW( 1000.0f,  1000.0f);    
        eikonal.set_NW( 1000.0f, 21000.0f);    
        eikonal.set_SE(21000.0f,  1000.0f);    

        eikonal.shots.n_xline = 2; 
        eikonal.shots.n_yline = 2; 

        eikonal.setGridGeometry(eikonal.shots);

        // Shots loop

        eikonal.T = new float[eikonal.nPointsB]();

        for (eikonal.shotId = 0; eikonal.shotId < eikonal.shots.all; eikonal.shotId++) 
        {   
            eikonal.arrivalFolder = "outputs/pod_extern_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 0;
            eikonal.eikonalComputing();

            eikonal.arrivalFolder = "outputs/fim_extern_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 1;
            eikonal.eikonalComputing();

            eikonal.arrivalFolder = "outputs/fsm_extern_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 2;
            eikonal.eikonalComputing();
        }
        
        delete[] eikonal.T;

        // Generating central shot

        eikonal.shots.elevation = 0.0f;

        eikonal.set_SW(11000.0f, 11000.0f);    
        eikonal.set_NW(11000.0f, 11000.0f);    
        eikonal.set_SE(11000.0f, 11000.0f);    

        eikonal.shots.n_xline = 1; 
        eikonal.shots.n_yline = 1; 

        eikonal.setGridGeometry(eikonal.shots);

        eikonal.T = new float[eikonal.nPointsB]();

        for (eikonal.shotId = 0; eikonal.shotId < eikonal.shots.all; eikonal.shotId++) 
        {   
            eikonal.arrivalFolder = "outputs/pod_central_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 0;
            eikonal.eikonalComputing();

            eikonal.arrivalFolder = "outputs/fim_central_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 1;
            eikonal.eikonalComputing();

            if (n == 2)
            {
                eikonal.exportTimesVolume = true;
                eikonal.eikonalFolder = "outputs/central_";
            }

            eikonal.arrivalFolder = "outputs/fsm_central_"+std::to_string((int) dh_all[n])+"m_";
            eikonal.eikonalType = 2;
            eikonal.eikonalComputing();
        }
        
        delete[] eikonal.T;
    }

    return 0;
}
