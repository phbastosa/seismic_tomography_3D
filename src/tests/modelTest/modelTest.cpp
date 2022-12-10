# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../essentials/model.hpp"

int main(int argc, char**argv)
{
    auto model = Model();
    auto utils = Utils();

    model.nb = 50;

    model.nx = 801;
    model.ny = 801;
    model.nz = 187;

    model.initialize();

    model.dx = 25.0f;
    model.dy = 25.0f;
    model.dz = 25.0f;

    std::cout<<"Reading velocity model with (z = "<<model.nz<<", x = "<<model.nx<<", y = "<<model.ny<<") samples"<<std::endl;

    model.fileName = "../../../inputs/models/overthrust_z187_x801_y801_h25.bin";

    float * Vp = utils.readBinaryFloat(model.fileName, model.nPoints);

    std::cout<<"Expanding velocity model dimensions to (z = "<<model.nzz<<", x = "<<model.nxx<<", y = "<<model.nyy<<") samples"<<std::endl;

    float * fullVp = model.expand(Vp);

    std::cout<<"Writing the expanded velocity model"<<std::endl;

    utils.writeBinaryFloat("outputs/expandedModel.bin", fullVp, model.nPointsB);

    std::cout<<"Reducing expanded velocity model back to (z = "<<model.nz<<", x = "<<model.nx<<", y = "<<model.ny<<") samples"<<std::endl;

    float * innerVp = model.reduce(fullVp);

    std::cout<<"Writing reduced velocity model"<<std::endl;

    utils.writeBinaryFloat("outputs/innerModel.bin", innerVp, model.nPoints);

    return 0;
}
