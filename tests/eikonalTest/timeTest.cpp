# include <chrono>
# include <iostream>

# include "../../src/eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    auto eikonal = Eikonal();

    eikonal.nb = 1;

    eikonal.nx = 401;
    eikonal.ny = 401;    
    eikonal.nz = 401;

    eikonal.dx = 25.0f;
    eikonal.dy = 25.0f;
    eikonal.dz = 25.0f;

    eikonal.initialize();

    eikonal.V = new float[eikonal.nPointsB];

    for (int i = 0; i < eikonal.nPointsB; ++i) 
    {
        eikonal.V[i] = 2000.0f;        
    }

    // Generate central shot

    float x = (float)(eikonal.nx/2)*eikonal.dx;
    float y = (float)(eikonal.ny/2)*eikonal.dy;

    eikonal.set_SW(x,y);     
    eikonal.set_NW(x,y);    
    eikonal.set_SE(x,y);    

    eikonal.shots.n_xline = 1; 
    eikonal.shots.n_yline = 1; 
    eikonal.shots.elevation = (float)(eikonal.nz/2)*eikonal.dz;
    
    eikonal.setGridGeometry(eikonal.shots);

    // Generate circular receiver survey
    
    eikonal.nodes.xcenter = x;
    eikonal.nodes.ycenter = y;
    eikonal.nodes.offsets = {(float) eikonal.dx*(eikonal.nx - 5) / 2.0f};
    eikonal.nodes.elevation = (float)(eikonal.nz/2)*eikonal.dz;
    eikonal.nodes.circle_spacing = 25.0f;
    
    eikonal.setCircularGeometry(eikonal.nodes);

    // Compressed loop

    std::cout<<"\nEikonal equation runtime test"<<std::endl;
    std::cout<<"\nModel dimensions:"<<std::endl;
    std::cout<<"Samples in x: "<<eikonal.nx<<" -> "<<(eikonal.nx-1)*eikonal.dx<<" m"<<std::endl;
    std::cout<<"Samples in y: "<<eikonal.ny<<" -> "<<(eikonal.ny-1)*eikonal.dy<<" m"<<std::endl;
    std::cout<<"Samples in z: "<<eikonal.nz<<" -> "<<(eikonal.nz-1)*eikonal.dz<<" m"<<std::endl;

    std::cout<<"\nCentral shot applied in position:"<<std::endl;
    std::cout<<"x = "<<eikonal.shots.x[0]<<" m"<<std::endl;
    std::cout<<"y = "<<eikonal.shots.y[0]<<" m"<<std::endl;
    std::cout<<"elevation = "<<eikonal.shots.z[0]<<" m"<<std::endl;

    std::cout<<"\nCircular geometry applied in configuration:"<<std::endl;
    std::cout<<"x center = "<<eikonal.nodes.xcenter<<" m"<<std::endl;
    std::cout<<"y center = "<<eikonal.nodes.ycenter<<" m"<<std::endl;
    std::cout<<"elevation = "<<eikonal.nodes.elevation<<" m"<<std::endl;
    std::cout<<"offset = "<<eikonal.nodes.offsets[0]<<" m"<<std::endl;

    std::cout<<"\n----------------- Run time ------------------\n"<<std::endl;

    // Comparing eikonal execution time

    eikonal.exportTimesVolume = true;
    eikonal.exportFirstArrivals = true;

    std::vector <std::string> labels {"outputs/pod_", "outputs/fim_", "outputs/fsm_"};

    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    eikonal.T = new float[eikonal.nPointsB]();

    for (int i = 0; i < labels.size(); i++)
    {
        eikonal.eikonalFolder = labels[i];
        eikonal.arrivalFolder = labels[i];

        ti = std::chrono::system_clock::now();

        eikonal.shotId = 0;   
        eikonal.eikonalType = i;
        eikonal.eikonalComputing();
        
        switch (i)
        {
        case 0:
            tf = std::chrono::system_clock::now();
            elapsed_seconds = tf - ti;
            std::cout<<"Podvin & Lecomte (1991) time = "<<elapsed_seconds.count()<<" s."<<std::endl;            
            break;
        case 1:
            tf = std::chrono::system_clock::now();
            elapsed_seconds = tf - ti;
            std::cout<<"Jeong & Whitaker (2008) time = "<<elapsed_seconds.count()<<" s."<<std::endl;
            break;
        case 2:
            tf = std::chrono::system_clock::now();
            elapsed_seconds = tf - ti;
            std::cout<<"Noble, Gesret and Belayouni (2014) time = "<<elapsed_seconds.count()<<" s."<<std::endl;
            break;
        }
    }

    delete[] eikonal.T;

    return 0;
}
