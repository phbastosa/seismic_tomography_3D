# include "eikonal.hpp"

# include "../../essentials/inout/inout.hpp"
# include "../../essentials/utils/utils.hpp"
# include "../../essentials/model/model.hpp"
# include "../../essentials/geometry/geometry.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal3D();
     
    eikonal.eikonalPath = "";
    eikonal.arrivalsPath = "";

    eikonal.exportTimesVolume = true;
    eikonal.exportFirstArrivals = true;

    // Setting model 

    auto m3D = Model3D();

    m3D.nx = 201;
    m3D.ny = 201;
    m3D.nz = 45;

    m3D.nb = 2;

    m3D.dx = 25.0f;
    m3D.dy = 25.0f;
    m3D.dz = 25.0f;

    m3D.init();

    m3D.vp = new float[m3D.nPointsB];

    m3D.readAndExpandVP("refractiveModel_45x201x201_25m.bin");

    // Setting geometry

    auto geom = Geometry3D();

    Utils::point2D SW, NW, SE;

    SW.x = 100.0f; SW.y = 2500.0f;    
    NW.x = 100.0f; NW.y = 2500.0f;    
    SE.x = 100.0f; SE.y = 2500.0f;    

    geom.setOBNS(SW,NW,SE,1,1,25.0f);

    SW.x = 4900.0f; SW.y = 100.0f;    
    NW.x = 4900.0f; NW.y = 4900.0f;    
    SE.x = 4900.0f; SE.y = 100.0f;    

    geom.setOBNR(SW,NW,SE,1,481,25.0f);

    geom.exportPositions("shots.txt","nodes.txt");

    // Shots loop

    eikonal.allocateVolumes(m3D);    

    for (int shot = 0; shot < geom.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.podvin3D(m3D, geom);
    }

    eikonal.deleteVolumes();

    return 0;
}
