# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../acoustic/acoustic.hpp"

int main(int argc, char**argv)
{
    auto acoustic = Acoustic();

    // Testing source 

    acoustic.nt = 1001;
    acoustic.dt = 0.001f;
    acoustic.fcut = 30.0f;
    acoustic.tlag = 0.15f;

    acoustic.sourceGenerator();

    acoustic.writeBinaryFloat("outputs/source_Nt"+ std::to_string(acoustic.nt) + ".bin", acoustic.source, acoustic.nt);

    // Testing Cerjan factor attenuation

    acoustic.nb = 50;
    acoustic.factor = 0.0045;

    acoustic.dampingGenerator();

    acoustic.writeBinaryFloat("outputs/damp1D_" + std::to_string(acoustic.nb) + "_samples.bin", acoustic.damp1D, acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp2D_" + std::to_string(acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp2D, acoustic.nb*acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp3D_" + std::to_string(acoustic.nb*acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp3D, acoustic.nb*acoustic.nb*acoustic.nb);




    return 0;
}
