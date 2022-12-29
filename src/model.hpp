# ifndef MODEL_HPP
# define MODEL_HPP

/* Returns the expanded volume based on input model properties */
float * expandModel(float * inputModel, int nx, int ny, int nz, int nb)
{
    int nxx = nx + 2*nb; 
    int nyy = ny + 2*nb; 
    int nzz = nz + 2*nb; 

    int nPoints = nxx*nyy*nzz;

    float * outputModel = new float[nPoints]();

    // Centering
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                outputModel[z + x*nzz + y*nxx*nzz] = inputModel[(z - nb) + (x - nb)*nz + (y - nb)*nx*nz];
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
                outputModel[z + x*nzz + y*nxx*nzz] = inputModel[0 + (x - nb)*nz + (y - nb)*nx*nz];
                outputModel[(nzz - z - 1) + x*nzz + y*nxx*nzz] = inputModel[(nz - 1) + (x - nb)*nz + (y - nb)*nx*nz];
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
                outputModel[z + x*nzz + y*nxx*nzz] = outputModel[z + nb*nzz + y*nxx*nzz];
                outputModel[z + (nxx - x - 1)*nzz + y*nxx*nzz] = outputModel[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
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
                outputModel[z + x*nzz + y*nxx*nzz] = outputModel[z + x*nzz + nb*nxx*nzz];
                outputModel[z + x*nzz + (nyy - y - 1)*nxx*nzz] = outputModel[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }

    return outputModel;
}

/* Returns the volume in original size without boundaries */
float * reduceModel(float * inputModel, int nx, int ny, int nz, int nb)
{
    int nxx = nx + 2*nb; 
    int nzz = nz + 2*nb; 

    int nPoints = nx*ny*nz;

    float * outputModel = new float[nPoints];

    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         // y direction
        int x = (int) (index - y*nx*nz) / nz;    // x direction
        int z = (int) (index - x*nz - y*nx*nz);  // z direction

        outputModel[z + x*nz + y*nx*nz] = inputModel[(z + nb) + (x + nb)*nzz + (y + nb)*nxx*nzz];
    }

    return outputModel;    
} 

# endif
