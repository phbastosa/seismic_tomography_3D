# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

# include "geometry.hpp"

int main(int argc, char**argv)
{
    auto G = Geometry3D();

    std::cout<<"\nPoint acquisition test"<<std::endl;

    Utils::point2D SW, NW, SE;

    SW.x = 10.8f; SW.y = 20.4f;    
    NW.x = 10.8f; NW.y = 20.4f;    
    SE.x = 10.8f; SE.y = 20.4f;    

    G.setOBNS(SW,NW,SE,1,1,1000.0f);

    std::cout<<G.shots->x[0]<<" "<<G.shots->y[0]<<" "<<G.shots->z[0]<<std::endl;

    SW.x = 15.8f; SW.y = 22.4f;    
    NW.x = 15.8f; NW.y = 22.4f;    
    SE.x = 15.8f; SE.y = 22.4f;    

    G.setOBNR(SW,NW,SE,1,1,2000.0f);    

    std::cout<<G.nodes->x[0]<<" "<<G.nodes->y[0]<<" "<<G.nodes->z[0]<<std::endl;

    // -------------------------------------------------------------------------------
    std::cout<<"\nHorizontal line acquisition test (x axis)"<<std::endl;

    SW.x = 10.8f;  SW.y = 20.4f;    
    NW.x = 10.8f;  NW.y = 20.4f;    
    SE.x = 350.8f; SE.y = 20.4f;    

    G.setOBNS(SW,NW,SE,5,1,1000.0f);

    for (int i = 0; i < G.ns; i++)
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
    
    std::cout<<"\n";

    SW.x = 15.8f;  SW.y = 22.4f;    
    NW.x = 15.8f;  NW.y = 22.4f;    
    SE.x = 350.8f; SE.y = 22.4f;    

    G.setOBNR(SW,NW,SE,5,1,2000.0f);    
    
    for (int i = 0; i < G.nr; i++)
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
    
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nVertical line acquisition test (y axis)"<<std::endl;

    SW.x = 10.8f;  SW.y = 20.4f;    
    NW.x = 10.8f;  NW.y = 350.4f;    
    SE.x = 10.8f;  SE.y = 20.4f;    

    G.setOBNS(SW,NW,SE,1,5,1000.0f);

    for (int i = 0; i < G.ns; i++)
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
    
    std::cout<<"\n";

    SW.x = 15.8f;  SW.y = 22.4f;    
    NW.x = 15.8f;  NW.y = 350.4f;    
    SE.x = 15.8f;  SE.y = 22.4f;    

    G.setOBNR(SW,NW,SE,1,5,2000.0f);    
    
    for (int i = 0; i < G.nr; i++)
        std::cout<<" ("<<G.nodes->x[i]<<", "<<G.nodes->y[i]<<", "<<G.nodes->z[i]<<") ";
    
    std::cout<<"\n";

    // -------------------------------------------------------------------------------
    std::cout<<"\nFull carpet acquisition test (x and y axis)"<<std::endl;

    SW.x = 0.0f;  SW.y = 0.0f;    
    NW.x = 0.0f;  NW.y = 50.0f;    
    SE.x = 50.0f; SE.y = 0.0f;    

    G.setOBNS(SW,NW,SE,10,5,1000.0f);

    for (int i = 0; i < G.ns; i++)
    {
        std::cout<<" ("<<G.shots->x[i]<<", "<<G.shots->y[i]<<", "<<G.shots->z[i]<<") ";
        if ((i+1) % 5 == 0) std::cout<<"\n";
    }
    std::cout<<"\n";

    SW.x = 10.0f;  SW.y = 10.0f;    
    NW.x = 10.0f;  NW.y = 40.0f;    
    SE.x = 40.0f;  SE.y = 10.0f;    

    G.setOBNR(SW,NW,SE,5,10,2000.0f);    
    
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
	
    return 0;
}
