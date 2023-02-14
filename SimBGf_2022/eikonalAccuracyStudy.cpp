# include <chrono>
# include <iostream>

# include "../../src/eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal();

    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    std::vector<std::string> modelNames {"outputs/refractiveModel_12x221x221_100m.bin", 
                                         "outputs/refractiveModel_23x441x441_50m.bin", 
                                         "outputs/refractiveModel_45x881x881_25m.bin"};

    std::vector<std::string> labels {"pod_","fim_","fsm_"};

    std::vector<int> nx_all {221, 441, 881};
    std::vector<int> ny_all {221, 441, 881};
    std::vector<int> nz_all {12, 23, 45};
    
    std::vector<float> dh_all {100.0f, 50.0f, 25.0f};

    eikonal.nb = 1;

    // Set fixed circular geometry

    eikonal.nodes.elevation = 0.0f;

    eikonal.nodes.xcenter = 11000.0f;
    eikonal.nodes.ycenter = 11000.0f;
    eikonal.nodes.offsets = {10000.0f};   // 10 km
    eikonal.nodes.circle_spacing = 12.5f;

    eikonal.setCircularGeometry(eikonal.nodes);

    // Setting five shot points  
    
    eikonal.shots.all = 5;

    eikonal.shots.x = new float[eikonal.shots.all];
    eikonal.shots.y = new float[eikonal.shots.all];
    eikonal.shots.z = new float[eikonal.shots.all];

    std::vector<float> x {1000, 21000, 1000, 21000, 11000};
    std::vector<float> y {1000, 1000, 21000, 21000, 11000};

    for (int k = 0; k < eikonal.shots.all; k++)
    {
        eikonal.shots.x[k] = x[k];
        eikonal.shots.y[k] = y[k];
        eikonal.shots.z[k] = 0.0f;
    }

    eikonal.geometryFolder = "outputs/";

    eikonal.exportPositions();

    eikonal.exportTimesVolume = false;
    eikonal.exportFirstArrivals = true;

    // Setting model 

    ti = std::chrono::system_clock::now();

    std::cout<<"\n";
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

        // Shots loop

        eikonal.T = new float[eikonal.nPointsB]();

        for (eikonal.shotId = 0; eikonal.shotId < eikonal.shots.all; eikonal.shotId++) 
        {   
            for (int k = 0; k < labels.size(); k++)
            {
                eikonal.eikonalType = k;

                eikonal.arrivalFolder = "outputs/" + labels[k] + std::to_string((int) dh_all[n]) + "m_";
                eikonal.eikonalFolder = "outputs/" + labels[k] + std::to_string((int) dh_all[n]) + "m_";

                if ((k == 0) && (eikonal.shotId == 4) && (n == 2)) // Podvin; shot 5; 25 m
                    eikonal.exportTimesVolume = true;

                eikonal.eikonalComputing();
            }    
        }
        
        delete[] eikonal.T;
        delete[] eikonal.V;
    }

    tf = std::chrono::system_clock::now();
    
    elapsed_seconds = tf - ti;

    std::cout<<"\n Experiment run time = "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}
