# include "model.hpp"

void Model::initialize()
{
    nxx = nx + 2 * nb;
    nyy = ny + 2 * nb;
    nzz = nz + 2 * nb;

    nPoints = nx * ny * nz;
    nPointsB = nxx * nyy * nzz;
}

float * Model::expandModel()
{
    float * expModel = new float[nPointsB];

    // Centering
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                expModel[z + x*nzz + y*nxx*nzz] = model[(z - nb) + (x - nb)*nz + (y - nb)*nx*nz];
            }
        }
    }

    // Z direction
    for (int z = 0; z < nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                expModel[z + x*nzz + y*nxx*nzz] = model[0 + (x - nb)*nz + (y - nb)*nx*nz];
                expModel[(nzz - z - 1) + x*nzz + y*nxx*nzz] = model[(nz - 1) + (x - nb)*nz + (y - nb)*nx*nz];
            }
        }
    }

    // X direction
    for (int x = 0; x < nb; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = nb; y < nyy - nb; y++)
            {
                expModel[z + x*nzz + y*nxx*nzz] = expModel[z + nb*nzz + y*nxx*nzz];
                expModel[z + (nxx - x - 1)*nzz + y*nxx*nzz] = expModel[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < nb; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                expModel[z + x*nzz + y*nxx*nzz] = expModel[z + x*nzz + nb*nxx*nzz];
                expModel[z + x*nzz + (nyy - y - 1)*nxx*nzz] = expModel[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }

    delete[] model;

    return expModel;
}
