# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

# include "geometry.hpp"

int main(int argc, char**argv)
{
    auto g3D = Geometry();

    std::cout<<"\nPoint acquisition test"<<std::endl;

    g3D.SW.x = 10.8f; g3D.SW.y = 20.4f;    
    g3D.NW.x = 10.8f; g3D.NW.y = 20.4f;    
    g3D.SE.x = 10.8f; g3D.SE.y = 20.4f;    

    g3D.shots.nx = 1; g3D.shots.ny = 1; g3D.shots.elevation = 25.0f;

    g3D.setGridShots();

    std::cout<<g3D.shots.x[0]<<" "<<g3D.shots.y[0]<<" "<<g3D.shots.z[0]<<std::endl;

    g3D.SW.x = 15.8f; g3D.SW.y = 22.4f;    
    g3D.NW.x = 15.8f; g3D.NW.y = 22.4f;    
    g3D.SE.x = 15.8f; g3D.SE.y = 22.4f;    

    g3D.nodes.nx = 1; g3D.nodes.ny = 1; g3D.nodes.elevation = 45.0f;

    g3D.setGridNodes();    

    std::cout<<g3D.nodes.x[0]<<" "<<g3D.nodes.y[0]<<" "<<g3D.nodes.z[0]<<std::endl;

    delete[] g3D.shots.x;
    delete[] g3D.shots.y;
    delete[] g3D.shots.z;

    delete[] g3D.nodes.x;
    delete[] g3D.nodes.y;
    delete[] g3D.nodes.z;

    // -------------------------------------------------------------------------------
    std::cout<<"\nHorizontal line acquisition test (x axis)"<<std::endl;

    g3D.SW.x = 10.8f;  g3D.SW.y = 20.4f;    
    g3D.NW.x = 10.8f;  g3D.NW.y = 20.4f;    
    g3D.SE.x = 350.8f; g3D.SE.y = 20.4f;    

    g3D.shots.nx = 5; g3D.shots.ny = 1; g3D.shots.elevation = 1000.0f;

    g3D.setGridShots();    
    
    for (int i = 0; i < g3D.shots.n; i++)
        std::cout<<" ("<<g3D.shots.x[i]<<", "<<g3D.shots.y[i]<<", "<<g3D.shots.z[i]<<") ";
    
    std::cout<<"\n";

    g3D.SW.x = 15.8f;  g3D.SW.y = 22.4f;    
    g3D.NW.x = 15.8f;  g3D.NW.y = 22.4f;    
    g3D.SE.x = 350.8f; g3D.SE.y = 22.4f;    

    g3D.nodes.nx = 5; g3D.nodes.ny = 1; g3D.nodes.elevation = 2000.0f;

    g3D.setGridNodes();    
    
    for (int i = 0; i < g3D.nodes.n; i++)
        std::cout<<" ("<<g3D.nodes.x[i]<<", "<<g3D.nodes.y[i]<<", "<<g3D.nodes.z[i]<<") ";
    
    std::cout<<"\n";

    delete[] g3D.shots.x;
    delete[] g3D.shots.y;
    delete[] g3D.shots.z;

    delete[] g3D.nodes.x;
    delete[] g3D.nodes.y;
    delete[] g3D.nodes.z;

    // -------------------------------------------------------------------------------
    std::cout<<"\nVertical line acquisition test (y axis)"<<std::endl;

    g3D.SW.x = 10.8f;  g3D.SW.y = 20.4f;    
    g3D.NW.x = 10.8f;  g3D.NW.y = 350.4f;    
    g3D.SE.x = 10.8f;  g3D.SE.y = 20.4f;    

    g3D.shots.nx = 1; g3D.shots.ny = 5; g3D.shots.elevation = 1000.0f;

    g3D.setGridShots();

    for (int i = 0; i < g3D.shots.n; i++)
        std::cout<<" ("<<g3D.shots.x[i]<<", "<<g3D.shots.y[i]<<", "<<g3D.shots.z[i]<<") ";
    
    std::cout<<"\n";

    g3D.SW.x = 15.8f;  g3D.SW.y = 22.4f;    
    g3D.NW.x = 15.8f;  g3D.NW.y = 350.4f;    
    g3D.SE.x = 15.8f;  g3D.SE.y = 22.4f;    

    g3D.nodes.nx = 1; g3D.nodes.ny = 5; g3D.nodes.elevation = 2000.0f;

    g3D.setGridNodes();    
    
    for (int i = 0; i < g3D.nodes.n; i++)
        std::cout<<" ("<<g3D.nodes.x[i]<<", "<<g3D.nodes.y[i]<<", "<<g3D.nodes.z[i]<<") ";
    
    std::cout<<"\n";

    delete[] g3D.shots.x;
    delete[] g3D.shots.y;
    delete[] g3D.shots.z;

    delete[] g3D.nodes.x;
    delete[] g3D.nodes.y;
    delete[] g3D.nodes.z;

    // -------------------------------------------------------------------------------
    std::cout<<"\nFull carpet acquisition test (x and y axis)"<<std::endl;

    g3D.SW.x = 0.0f;  g3D.SW.y = 0.0f;    
    g3D.NW.x = 0.0f;  g3D.NW.y = 50.0f;    
    g3D.SE.x = 50.0f; g3D.SE.y = 0.0f;    

    g3D.shots.nx = 10; g3D.shots.ny = 5; g3D.shots.elevation = 1000.0f;

    g3D.setGridShots();

    for (int i = 0; i < g3D.shots.n; i++)
    {
        std::cout<<" ("<<g3D.shots.x[i]<<", "<<g3D.shots.y[i]<<", "<<g3D.shots.z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    g3D.SW.x = 10.0f;  g3D.SW.y = 10.0f;    
    g3D.NW.x = 10.0f;  g3D.NW.y = 40.0f;    
    g3D.SE.x = 40.0f;  g3D.SE.y = 10.0f;    

    g3D.nodes.nx = 10; g3D.nodes.ny = 5; g3D.nodes.elevation = 1000.0f;

    g3D.setGridNodes();    
    
    for (int i = 0; i < g3D.nodes.n; i++)
    {    
        std::cout<<" ("<<g3D.nodes.x[i]<<", "<<g3D.nodes.y[i]<<", "<<g3D.nodes.z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nReciprocity test"<<std::endl;

    g3D.setReciprocity();

    for (int i = 0; i < g3D.shots.n; i++)
    {
        std::cout<<" ("<<g3D.shots.x[i]<<", "<<g3D.shots.y[i]<<", "<<g3D.shots.z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    for (int i = 0; i < g3D.nodes.n; i++)
    {    
        std::cout<<" ("<<g3D.nodes.x[i]<<", "<<g3D.nodes.y[i]<<", "<<g3D.nodes.z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";
	
    delete[] g3D.shots.x;
    delete[] g3D.shots.y;
    delete[] g3D.shots.z;

    delete[] g3D.nodes.x;
    delete[] g3D.nodes.y;
    delete[] g3D.nodes.z;

    // Testing circular shots

    g3D.shots.xc = 2000.0f;
    g3D.shots.yc = 2000.0f;
    g3D.shots.ds = 50.0f;
    
    g3D.shots.offsets = {1e3f, 1.5e3f, 2e3f};

    g3D.shots.elevation = 10.0f;

    g3D.setCircularShots();

    g3D.SW.x = 500.0f;  g3D.SW.y = 500.0f;    
    g3D.NW.x = 500.0f;  g3D.NW.y = 3500.0f;    
    g3D.SE.x = 3500.0f; g3D.SE.y = 50.0f;    

    g3D.nodes.nx = 11; g3D.nodes.ny = 11; g3D.nodes.elevation = 1000.0f;

    g3D.setGridNodes();    

    g3D.nodesPath = "nodes_n" + std::to_string(g3D.nodes.n) + ".bin";
    g3D.shotsPath = "shots_n" + std::to_string(g3D.shots.n) + ".bin";

    g3D.exportPositions();

    return 0;
}
