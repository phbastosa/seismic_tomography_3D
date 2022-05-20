# include "model.hpp"

# include "../inout/inout.hpp"

int main(int argc, char **argv)
{
    auto m3D = Model3D();

    m3D.nx = 338;
    m3D.ny = 338;
    m3D.nz = 88;

    m3D.nb = 20;

    m3D.dx = 50.0f;
    m3D.dy = 50.0f;
    m3D.dz = 50.0f;

    m3D.vpPath = "saltDome3D_z88_x338_y338.bin";

    m3D.readAndExpandVP();

    std::string newPath = "saltDome3D_expand_test_z"+ InOut::toString(m3D.nzz) +"_x"+ InOut::toString(m3D.nxx) +"_y"+ InOut::toString(m3D.nyy) +".bin";

    InOut::writeBinaryFloat(newPath,m3D.vp,m3D.nPointsB);

    std::cout<<"\n Expanded model written..."<<std::endl;

    return 0;
}  
