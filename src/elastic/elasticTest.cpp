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

    elastic.setWavelet();

    elastic.set_SW(500.0f, 500.0f);    
    elastic.set_NW(500.0f, 500.0f);    
    elastic.set_SE(500.0f, 500.0f);   

    elastic.shots.n_xline = 1;
    elastic.shots.n_yline = 1;
    elastic.shots.elevation = 500.0f;

    elastic.setGridShots();

    elastic.set_SW(0.0f, 0.0f);    
    elastic.set_NW(0.0f, 1000.0f);    
    elastic.set_SE(1000.0f, 1000.0f);   

    elastic.nodes.n_xline = 21;
    elastic.nodes.n_yline = 21;
    elastic.nodes.elevation = 500.0f;

    elastic.setGridNodes();
    
    elastic.vp = new float[elastic.nPointsB];
    elastic.vs = new float[elastic.nPointsB];
    elastic.rho = new float[elastic.nPointsB];
    
    for (int i = 0; i < elastic.nPointsB; i++)
    {
        elastic.vp[i] = 1500.0f;
        elastic.vs[i] = 0.0f;
        elastic.rho[i] = 1000.0f;
    }

    elastic.forwardModeling();

    elastic.writeBinaryFloat("seismogram.bin", elastic.seismogram, elastic.nt*elastic.nodes.all);

    return 0;
}
