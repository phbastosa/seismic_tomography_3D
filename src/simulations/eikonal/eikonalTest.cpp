# include "eikonal.hpp"

# include "../../essentials/inout/inout.hpp"
# include "../../essentials/utils/utils.hpp"
# include "../../essentials/model/model.hpp"
# include "../../essentials/geometry/geometry.hpp"

int main(int argc, char **argv)
{
    // Setting eikonal

    auto eikonal = Eikonal3D();
     
    // Setting model 

    eikonal.m3D.nx = 441;
    eikonal.m3D.ny = 441;
    eikonal.m3D.nz = 23;

    eikonal.m3D.nb = 2;

    eikonal.m3D.dx = 50.0f;
    eikonal.m3D.dy = 50.0f;
    eikonal.m3D.dz = 50.0f;

    eikonal.m3D.init();

    eikonal.m3D.vp = new float[eikonal.m3D.nPointsB];

    eikonal.m3D.vpPath = "refractiveModel_23x441x441_50m.bin";

    eikonal.m3D.readAndExpandVP();

    // Setting geometry

    eikonal.g3D.circles.xc = 11000.0f;
    eikonal.g3D.circles.yc = 11000.0f;
    eikonal.g3D.circles.ds = 50.0f;
    
    eikonal.g3D.circles.offsets = {1e4f};

    eikonal.g3D.sElev = 10.0f;

    eikonal.g3D.setCircularShots();

    eikonal.g3D.SW.x = 11000.0f; eikonal.g3D.SW.y = 11000.0f;    
    eikonal.g3D.NW.x = 11000.0f; eikonal.g3D.NW.y = 11000.0f;    
    eikonal.g3D.SE.x = 11000.0f; eikonal.g3D.SE.y = 11000.0f;    

    eikonal.g3D.nrx = 1; 
    eikonal.g3D.nry = 1; 
    eikonal.g3D.rElev = 10.0f;

    eikonal.g3D.setGridNodes();

    eikonal.g3D.setReciprocity();

    eikonal.g3D.nodesPath = "nodes_n" + InOut::toString(eikonal.g3D.nr) + ".txt";
    eikonal.g3D.shotsPath = "shots_n" + InOut::toString(eikonal.g3D.ns) + ".txt";

    eikonal.g3D.exportPositions();

    // Shots loop

    eikonal.exportTimesVolume = true;
    eikonal.exportFirstArrivals = true;

    eikonal.eikonalPath = "podvin_";
    eikonal.arrivalsPath = "podvin_";

    eikonal.allocateVolumes();    

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.podvin3D();
    }

    eikonal.eikonalPath = "fim_";
    eikonal.arrivalsPath = "fim_";

    for (int shot = 0; shot < eikonal.g3D.ns; shot++)
    {
        eikonal.shotId = shot;
        eikonal.fim3D();
    }

    eikonal.deleteVolumes();

    return 0;
}
