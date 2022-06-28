# include "model.hpp"

void Model::initialize(int nx, int ny, int nz, int nb)
{
    nxx = nx + 2 * nb;
    nyy = ny + 2 * nb;
    nzz = nz + 2 * nb;

    nPoints = nx * ny * nz;
    nPointsB = nxx * nyy * nzz;
}

float * Model::expandVolume(float * volume, int nx, int ny, int nz, int nb)
{
    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    float * expVolume = new float[nxx * nyy * nzz]();

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

float * Model::reduceVolume(float * expVolume, int nx, int ny, int nz, int nb)
{
    int nxb = nx + 2*nb;
    int nzb = nz + 2*nb;

    float * volume = new float[nx*ny*nz];

    for (int index = 0; index < nx*ny*nz; index++)
    {
        int y = (int) (index / (nx*nz));         // y direction
        int x = (int) (index - y*nx*nz) / nz;    // x direction
        int z = (int) (index - x*nz - y*nx*nz);  // z direction

        volume[z + x*nz + y*nx*nz] = expVolume[(z + nb) + (x + nb)*nzb + (y + nb)*nxb*nzb];
    }

    delete[] expVolume;

    return volume;
}
