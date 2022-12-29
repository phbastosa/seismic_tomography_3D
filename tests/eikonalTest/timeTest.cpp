# include <chrono>
# include <iostream>

# include "../../src/eikonal.hpp"

int main(int argc, char **argv)
{
    int nx = 401;
    int ny = 401;    
    int nz = 401;

    int dx = 25.0f;
    int dy = 25.0f;
    int dz = 25.0f;

    float * V = new float[nx*ny*nz];

    for (int index = 0; index < nx*ny*nz; index++) 
    {
        V[index] = 1500.0f;
    }

    // Generate central shot

    Point SE, SW, NW;

    SW.x = (float)(nx/2)*dx; SW.y = (float)(ny/2)*dy;     
    NW.x = (float)(nx/2)*dx; NW.y = (float)(ny/2)*dy;    
    SE.x = (float)(nx/2)*dx; SE.y = (float)(ny/2)*dy;    

    float elevation = (float)(nz/2)*dz;

    auto shots = setGridGeometry(SE,SW,NW,1,1,elevation);

    // Generate circular receiver survey
    
    float xcenter = (float)(nx/2)*dx;
    float ycenter = (float)(ny/2)*dy;
    float spacing = 12.5f;

    std::vector<float> offsets = {2000.0f};

    auto nodes = setCircularGeometry(offsets,spacing,xcenter,ycenter,elevation);

    // Compressed loop

    std::cout<<"\nEikonal equation runtime test"<<std::endl;
    std::cout<<"\nModel dimensions:"<<std::endl;
    std::cout<<"Samples in x: "<<nx<<" -> "<<(nx-1)*dx<<" m"<<std::endl;
    std::cout<<"Samples in y: "<<ny<<" -> "<<(ny-1)*dy<<" m"<<std::endl;
    std::cout<<"Samples in z: "<<nz<<" -> "<<(nz-1)*dz<<" m"<<std::endl;

    std::cout<<"\nCentral shot applied in position:"<<std::endl;
    std::cout<<"x = "<<shots.x[0]<<" m"<<std::endl;
    std::cout<<"y = "<<shots.y[0]<<" m"<<std::endl;
    std::cout<<"elevation = "<<shots.z[0]<<" m"<<std::endl;

    std::cout<<"\nCircular geometry applied in configuration:"<<std::endl;
    std::cout<<"x center = "<<xcenter<<" m"<<std::endl;
    std::cout<<"y center = "<<ycenter<<" m"<<std::endl;
    std::cout<<"elevation = "<<elevation<<" m"<<std::endl;
    std::cout<<"offset = "<<offsets[0]<<" m"<<std::endl;

    std::cout<<"\n----------------- Run time ------------------\n"<<std::endl;

    // Comparing eikonal execution time

    bool writeTT = true;
    bool writeFA = true;    

    std::string folder;

    float sx = shots.x[0];
    float sy = shots.y[0];
    float sz = shots.z[0];

    std::vector <std::string> labels = {"pod_", "fim_", "fsm_"};
    std::vector <std::string> formulation = {"Podvin & Lecomte (1991) formulation",
                                             "Jeong & Witaker (2008) formulation", 
                                             "Noble, Gesret and Belayouni (2014) formulation"};

    for (int type = 0; type < labels.size(); type++) // eikonal type
    {
        folder = "outputs/" + labels[type]; 

        auto ti = std::chrono::steady_clock::now();

        float * T = eikonalComputing(V,nx,ny,nz,dx,dy,dz,sx,sy,sz,type);
        
        if (writeTT) writeTravelTimes(T,nx,ny,nz,0,folder);    
        if (writeFA) writeFirstArrivals(T,nodes,nx,ny,nz,dx,dy,dz,0,folder); 

        auto tf = std::chrono::steady_clock::now();

        std::cout<<formulation[type]<<std::endl;
        std::cout << "Elapsed time: "
                    << std::chrono::duration_cast<std::chrono::seconds>(tf - ti).count()
                    << " s\n"<< std::endl;
    }

    return 0;
}
