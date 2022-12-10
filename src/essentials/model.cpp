# include "model.hpp"

void Model::initialize()
{
    nxx = nx + 2 * nb;
    nyy = ny + 2 * nb;
    nzz = nz + 2 * nb;

    nPoints = nx * ny * nz;
    nPointsB = nxx * nyy * nzz;
}

float * Model::expand(float * volume)
{
    float * expVolume = new float[nPointsB]();

    // Centering
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                expVolume[z + x*nzz + y*nxx*nzz] = volume[(z - nb) + (x - nb)*nz + (y - nb)*nx*nz];
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
                expVolume[z + x*nzz + y*nxx*nzz] = volume[0 + (x - nb)*nz + (y - nb)*nx*nz];
                expVolume[(nzz - z - 1) + x*nzz + y*nxx*nzz] = volume[(nz - 1) + (x - nb)*nz + (y - nb)*nx*nz];
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
                expVolume[z + x*nzz + y*nxx*nzz] = expVolume[z + nb*nzz + y*nxx*nzz];
                expVolume[z + (nxx - x - 1)*nzz + y*nxx*nzz] = expVolume[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
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
                expVolume[z + x*nzz + y*nxx*nzz] = expVolume[z + x*nzz + nb*nxx*nzz];
                expVolume[z + x*nzz + (nyy - y - 1)*nxx*nzz] = expVolume[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }

    delete[] volume;

    return expVolume;
}

float * Model::reduce(float * expVolume)
{
    float * volume = new float[nPoints];

    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         // y direction
        int x = (int) (index - y*nx*nz) / nz;    // x direction
        int z = (int) (index - x*nz - y*nx*nz);  // z direction

        volume[z + x*nz + y*nx*nz] = expVolume[(z + nb) + (x + nb)*nzz + (y + nb)*nxx*nzz];
    }

    delete[] expVolume;

    return volume;
}
