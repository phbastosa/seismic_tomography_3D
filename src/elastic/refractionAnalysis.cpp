# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

# include "elastic.hpp"

int main(int argc, char **argv)
{
    auto elastic = Elastic();

    elastic.nt = 7001;
    elastic.dt = 1e-3f;    

    elastic.m3D.nx = 2001;
    elastic.m3D.ny = 21;
    elastic.m3D.nz = 451;
    elastic.m3D.nb = 50;

    elastic.m3D.dx = 10.0f;
    elastic.m3D.dy = 10.0f;
    elastic.m3D.dz = 10.0f;

    elastic.m3D.initialize();

    elastic.nsrc = 201;
    elastic.source = new float[elastic.nsrc];
    elastic.io.readBinaryFloat("../../inputs/wavelets/sourceZeroPhase_201_1ms.bin", elastic.source, elastic.nsrc);
    // elastic.io.readBinaryFloat("../../inputs/wavelets/sourceMinPhase_201_1ms.bin", elastic.source, elastic.nsrc);

    elastic.g3D.SW.x = 0.0f; elastic.g3D.SW.y = 100.0f;    
    elastic.g3D.NW.x = 0.0f; elastic.g3D.NW.y = 100.0f;    
    elastic.g3D.SE.x = 0.0f; elastic.g3D.SE.y = 100.0f;   

    elastic.g3D.shots.nx = 1;
    elastic.g3D.shots.ny = 1;
    elastic.g3D.shots.elevation = 0.0f;

    elastic.g3D.setGridShots();

    elastic.g3D.SW.x =     0.0f; elastic.g3D.SW.y = 100.0f;    
    elastic.g3D.NW.x =     0.0f; elastic.g3D.NW.y = 100.0f;    
    elastic.g3D.SE.x = 20000.0f; elastic.g3D.SE.y = 100.0f;   

    elastic.g3D.nodes.nx = 1001;
    elastic.g3D.nodes.ny = 1;
    elastic.g3D.nodes.elevation = 0.0f;

    elastic.g3D.setGridNodes();

    elastic.cerjan = new float[elastic.m3D.nPointsB];
    elastic.m3D.vp = new float[elastic.m3D.nPointsB];
    elastic.m3D.vs = new float[elastic.m3D.nPointsB];
    elastic.m3D.rho = new float[elastic.m3D.nPointsB];

    elastic.cerjan = elastic.m3D.readAndExpandModel("../../inputs/models/cerjanVolumeABC.bin");
    elastic.m3D.vp = elastic.m3D.readAndExpandModel("../../inputs/models/vpModel_451x2001x21_10m.bin");
    elastic.m3D.vs = elastic.m3D.readAndExpandModel("../../inputs/models/vsModel_451x2001x21_10m.bin");
    elastic.m3D.rho = elastic.m3D.readAndExpandModel("../../inputs/models/rhoModel_451x2001x21_10m.bin");

    elastic.forwardModeling();

    elastic.io.writeBinaryFloat("seismogram.bin",elastic.seismogram,elastic.nt*elastic.g3D.nodes.n);

    return 0;
}
