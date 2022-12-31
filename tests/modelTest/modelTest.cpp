# include <string>
# include <iostream>

# include "../../src/essentials/utils.hpp"
# include "../../src/essentials/model.hpp"

int main(int argc, char**argv)
{
    auto model = Model();
    auto utils = Utils();

    model.nb = 10;

    model.nx = 201;
    model.ny = 201;
    model.nz = 51;

    model.initialize();

    model.dx = 25.0f;
    model.dy = 25.0f;
    model.dz = 25.0f;

    float * Vp = new float[model.nPoints]();

    for (int i = 0; i < model.nz; i++)
    {
        for (int j = 0; j < model.nx; j++)
        {
            for (int k = 0; k < model.ny; k++)
            {
                Vp[i + j*model.nz + k*model.nx*model.nz] = 1500 + 0.6*i*model.dz + 0.2*j*model.dx + 0.1*k*model.dy; 
            }
        }
    }

    std::cout<<"Expanding velocity model dimensions to (z = "<<model.nzz<<", x = "<<model.nxx<<", y = "<<model.nyy<<") samples"<<std::endl;

    float * fullVp = model.expand(Vp);

    std::cout<<"Writing the expanded velocity model"<<std::endl;

    utils.writeBinaryFloat("outputs/expandedModel.bin", fullVp, model.nPointsB);

    std::cout<<"Reducing expanded velocity model back to (z = "<<model.nz<<", x = "<<model.nx<<", y = "<<model.ny<<") samples"<<std::endl;

    float * innerVp = model.reduce(fullVp);

    std::cout<<"Writing reduced velocity model"<<std::endl;

    utils.writeBinaryFloat("outputs/innerModel.bin", innerVp, model.nPoints);

    delete[] Vp;
    delete[] fullVp;
    delete[] innerVp;

    return 0;
}
