# include "model.hpp"

# include "../inout/inout.hpp"

int main(int argc, char **argv)
{
    auto m3D = Model();

    m3D.nx = 338;
    m3D.ny = 338;
    m3D.nz = 88;

    m3D.nb = 20;

    m3D.dx = 50.0f;
    m3D.dy = 50.0f;
    m3D.dz = 50.0f;

    m3D.initialize();

    m3D.vpPath = "saltDome3D_z88_x338_y338.bin";

    m3D.vp = m3D.readAndExpandModel(m3D.vpPath);

    std::string newPath = "saltDome3D_expand_test_z" + std::to_string(m3D.nzz) + "_x" + std::to_string(m3D.nxx) + "_y" + std::to_string(m3D.nyy) + ".bin";
 
    InOut::writeBinaryFloat(newPath,m3D.vp,m3D.nPointsB);

    std::cout<<"\n Expanded model written..."<<std::endl;

    return 0;
}  

