# include <chrono>
# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../acoustic/acoustic.hpp"

int main(int argc, char **argv)
{
    auto acoustic = Acoustic();

    acoustic.paramFile = argv[1];

    auto start = std::chrono::steady_clock::now();
    
    acoustic.setParameters();

    // Testing source wavelet

    acoustic.writeBinaryFloat("outputs/source_Nt"+ std::to_string(acoustic.nsrc) + ".bin", acoustic.source, acoustic.nsrc);

    // Testing Cerjan factor attenuation

    acoustic.writeBinaryFloat("outputs/damp1D_" + std::to_string(acoustic.nb) + "_samples.bin", acoustic.damp1D, acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp2D_" + std::to_string(acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp2D, acoustic.nb*acoustic.nb);
    acoustic.writeBinaryFloat("outputs/damp3D_" + std::to_string(acoustic.nb*acoustic.nb*acoustic.nb) + "_samples.bin", acoustic.damp3D, acoustic.nb*acoustic.nb*acoustic.nb);
        
    // Preparing shots loop

    // acoustic.forwardModeling();

    auto end = std::chrono::steady_clock::now();

    std::cout << "\n\nElapsed time in seconds: "
              << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
              << " s"<< std::endl;

    return 0;
}

