# include <string>
# include <iostream>

# include "../../src/utils.hpp"
# include "../../src/model.hpp"

int main(int argc, char**argv)
{
    int nb = 10;

    int nx = 201;
    int ny = 201;
    int nz = 51;
    
    float dx = 25.0f;
    float dy = 25.0f;
    float dz = 25.0f;

    float * Vp = new float[nx*ny*nz]();

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            for (int k = 0; k < ny; k++)
            {
                Vp[i + j*nz + k*nx*nz] = 1500 + 0.6*i*dz + 0.2*j*dx + 0.1*k*dy; 
            }
        }
    }

    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    std::cout<<"Expanding velocity model dimensions to (z = "<<nzz<<", x = "<<nxx<<", y = "<<nyy<<") samples"<<std::endl;

    float * fullVp = expandModel(Vp, nx, ny, nz, nb);

    std::cout<<"Writing the expanded velocity model"<<std::endl;

    writeBinaryFloat("outputs/expandedModel.bin", fullVp, nxx*nyy*nzz);

    std::cout<<"Reducing expanded velocity model back to (z = "<<nz<<", x = "<<nx<<", y = "<<ny<<") samples"<<std::endl;

    float * innerVp = reduceModel(fullVp, nx, ny, nz, nb);

    std::cout<<"Writing reduced velocity model"<<std::endl;

    writeBinaryFloat("outputs/innerModel.bin", innerVp, nx*ny*nz);

    delete[] Vp;
    delete[] fullVp;
    delete[] innerVp;

    return 0;
}
