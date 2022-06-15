# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

# include "elastic.hpp"

int main(int argc, char **argv)
{
    auto elastic = Elastic();

    elastic.nt = 1001;
    elastic.dt = 1e-3f;    

    elastic.nx = 101;
    elastic.ny = 101;
    elastic.nz = 101;
    elastic.nb = 50;

    elastic.dx = 10.0f;
    elastic.dy = 10.0f;
    elastic.dz = 10.0f;

    elastic.initialize();

    elastic.nsrc = 201;
    // elastic.source = elastic.readBinaryFloat("../../inputs/wavelets/sourceMinPhase_201_1ms.bin", elastic.nsrc);
    elastic.source = elastic.readBinaryFloat("../../inputs/wavelets/sourceZeroPhase_201_1ms.bin", elastic.nsrc);

    elastic.set_SW(500.0f, 500.0f);    
    elastic.set_NW(500.0f, 500.0f);    
    elastic.set_SE(500.0f, 500.0f);   

    elastic.shots.n_xline = 1;
    elastic.shots.n_yline = 1;
    elastic.shots.elevation = 0.0f;

    elastic.setGridShots();

    elastic.set_SW(0.0f, 500.0f);    
    elastic.set_NW(0.0f, 500.0f);    
    elastic.set_SE(1000.0f, 500.0f);   

    elastic.nodes.n_xline = 21;
    elastic.nodes.n_yline = 1;
    elastic.nodes.elevation = 0.0f;

    elastic.setGridNodes();
    
    elastic.vp = elastic.expandModel(elastic.readBinaryFloat("../../inputs/models/vpModel_101x101x101_10m.bin", elastic.nPoints));
    elastic.vs = elastic.expandModel(elastic.readBinaryFloat("../../inputs/models/vsModel_101x101x101_10m.bin", elastic.nPoints));
    elastic.rho = elastic.expandModel(elastic.readBinaryFloat("../../inputs/models/rhoModel_101x101x101_10m.bin", elastic.nPoints));
        
    elastic.cerjan = elastic.readBinaryFloat("../../inputs/models/cerjanVolumeABC.bin", elastic.nPointsB);
    
    elastic.forwardModeling();

    elastic.writeBinaryFloat("seismogram.bin", elastic.seismogram, elastic.nt*elastic.nodes.all);

    return 0;
}
