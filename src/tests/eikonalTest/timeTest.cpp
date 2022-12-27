# include <omp.h>
# include <iostream>

# include "../../eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    auto eikonal = Eikonal();

    eikonal.nx = 881;
    eikonal.ny = 881;    
    eikonal.nz = 45;

    eikonal.dx = 25.0f;
    eikonal.dy = 25.0f;
    eikonal.dz = 25.0f;

    eikonal.nb = 2;

    eikonal.initialize();

    eikonal.V = new float[eikonal.nPointsB];

    int interface = (int)(1000.0f/eikonal.dz) + eikonal.nxx*eikonal.nzz + eikonal.nyy*eikonal.nxx*eikonal.nzz;

    for (int i = 0; i < eikonal.nPointsB; ++i) 
    {
        eikonal.V[i] = 1500.0f;
        
        if (i >= interface) eikonal.V[i] = 2000.0f;
    }

    // Generate central shot

    eikonal.set_SW(11000.0f, 11000.0f);     
    eikonal.set_NW(11000.0f, 11000.0f);    
    eikonal.set_SE(11000.0f, 11000.0f);    

    eikonal.shots.n_xline = 1; 
    eikonal.shots.n_yline = 1; 
    eikonal.shots.elevation = 0.0f;
    
    eikonal.setGridGeometry(eikonal.shots);

    // Generate circular receiver survey
    
    eikonal.nodes.xcenter = 11000.0f;
    eikonal.nodes.ycenter = 11000.0f;
    eikonal.nodes.offsets = {10000.0f};
    eikonal.nodes.elevation = 0.0f;
    eikonal.nodes.circle_spacing = 12.5f;
    
    eikonal.setCircularGeometry(eikonal.nodes);

    // Compressed loop

    int n = 3;
    double t0;

    eikonal.exportTimesVolume = false;
    eikonal.exportFirstArrivals = false;

    std::vector <std::string> labels {"pod_", "fim_", "fsm_"};

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

    eikonal.T = new float[eikonal.nPointsB]();

    for (int i = 0; i < n; i++)
    {
        eikonal.eikonalFolder = labels[i];
        eikonal.arrivalFolder = labels[i];

        t0 = omp_get_wtime();

        eikonal.shotId = 0;   
        eikonal.eikonalType = i;
        eikonal.eikonalComputing();
        
        switch (i)
        {
        case 0:
            std::cout<<"Podvin & Lecomte (1991) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;
            break;
        case 1:
            std::cout<<"Jeong & Whitaker (2008) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;
            break;
        case 2:
            std::cout<<"Noble, Gesret and Belayouni (2014) time = "<<omp_get_wtime() - t0<<" s."<<std::endl;
            break;
        }
    }

    delete[] eikonal.T;

    return 0;
}
