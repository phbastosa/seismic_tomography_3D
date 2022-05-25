# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

# include "geometry.hpp"

int main(int argc, char**argv)
{
    auto G = Geometry();

    std::cout<<"\nPoint acquisition test"<<std::endl;

    G.SW.x = 10.8f; G.SW.y = 20.4f;    
    G.NW.x = 10.8f; G.NW.y = 20.4f;    
    G.SE.x = 10.8f; G.SE.y = 20.4f;    

    G.nsx = 1; G.nsy = 1; G.sElev = 25.0f;

    G.setGridShots();

    std::cout<<G.shots->x[0]<<" "<<G.shots->y[0]<<" "<<G.shots->z[0]<<std::endl;

    G.SW.x = 15.8f; G.SW.y = 22.4f;    
    G.NW.x = 15.8f; G.NW.y = 22.4f;    
    G.SE.x = 15.8f; G.SE.y = 22.4f;    

    G.nrx = 1; G.nry = 1; G.rElev = 45.0f;

    G.setGridNodes();    

    std::cout<<G.nodes->x[0]<<" "<<G.nodes->y[0]<<" "<<G.nodes->z[0]<<std::endl;

    // -------------------------------------------------------------------------------
    std::cout<<"\nHorizontal line acquisition test (x axis)"<<std::endl;

    G.SW.x = 10.8f;  G.SW.y = 20.4f;    
    G.NW.x = 10.8f;  G.NW.y = 20.4f;    
    G.SE.x = 350.8f; G.SE.y = 20.4f;    

    G.nsx = 5; G.nsy = 1; G.sElev = 1000.0f;

    G.setGridShots();    
    
    for (int i = 0; i < G.ns; i++)
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
    
    std::cout<<"\n";

    G.SW.x = 15.8f;  G.SW.y = 22.4f;    
    G.NW.x = 15.8f;  G.NW.y = 22.4f;    
    G.SE.x = 350.8f; G.SE.y = 22.4f;    

    G.nrx = 5; G.nry = 1; G.rElev = 2000.0f;

    G.setGridNodes();    
    
    for (int i = 0; i < G.nr; i++)
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
    
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nVertical line acquisition test (y axis)"<<std::endl;

    G.SW.x = 10.8f;  G.SW.y = 20.4f;    
    G.NW.x = 10.8f;  G.NW.y = 350.4f;    
    G.SE.x = 10.8f;  G.SE.y = 20.4f;    

    G.nsx = 1; G.nsy = 5; G.sElev = 1000.0f;

    G.setGridShots();

    for (int i = 0; i < G.ns; i++)
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
    
    std::cout<<"\n";

    G.SW.x = 15.8f;  G.SW.y = 22.4f;    
    G.NW.x = 15.8f;  G.NW.y = 350.4f;    
    G.SE.x = 15.8f;  G.SE.y = 22.4f;    

    G.nrx = 1; G.nry = 5; G.rElev = 2000.0f;

    G.setGridNodes();    
    
    for (int i = 0; i < G.nr; i++)
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
    
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nFull carpet acquisition test (x and y axis)"<<std::endl;

    G.SW.x = 0.0f;  G.SW.y = 0.0f;    
    G.NW.x = 0.0f;  G.NW.y = 50.0f;    
    G.SE.x = 50.0f; G.SE.y = 0.0f;    

    G.nsx = 10; G.nsy = 5; G.sElev = 1000.0f;

    G.setGridShots();

    for (int i = 0; i < G.ns; i++)
    {
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    G.SW.x = 10.0f;  G.SW.y = 10.0f;    
    G.NW.x = 10.0f;  G.NW.y = 40.0f;    
    G.SE.x = 40.0f;  G.SE.y = 10.0f;    

    G.nrx = 10; G.nry = 5; G.rElev = 1000.0f;

    G.setGridNodes();    
    
    for (int i = 0; i < G.nr; i++)
    {    
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nReciprocity test"<<std::endl;

    G.setReciprocity();

    for (int i = 0; i < G.ns; i++)
    {
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    for (int i = 0; i < G.nr; i++)
    {    
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";
	
    G.circles.xc = 2000.0f;
    G.circles.yc = 2000.0f;
    G.circles.ds = 50.0f;
    
    G.circles.offsets = {1e3f, 1.5e3f, 2e3f};

    G.sElev = 10.0f;

    G.setCircularShots();

    G.SW.x = 500.0f;  G.SW.y = 500.0f;    
    G.NW.x = 500.0f;  G.NW.y = 3500.0f;    
    G.SE.x = 3500.0f; G.SE.y = 50.0f;    

    G.nrx = 11; G.nry = 11; G.rElev = 1000.0f;

    G.setGridNodes();    

    G.nodesPath = "nodes_n" + std::to_string(G.nr) + ".bin";
    G.shotsPath = "shots_n" + std::to_string(G.ns) + ".bin";

    G.exportPositions();

    return 0;
}
