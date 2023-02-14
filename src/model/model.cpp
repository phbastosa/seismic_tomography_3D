# include <cmath>

# include "model.hpp"

# include "../file_manager/file_manager.hpp"
# include "../utils/smoothing/gaussian.hpp"
# include "../utils/interpolation/trilinear.hpp"

void Model::set_parameters(std::string file)
{
    auto fm = File_manager();

    fm.parameter_file = file;

    nx = std::stoi(fm.catch_parameter("nx"));
    ny = std::stoi(fm.catch_parameter("ny"));
    nz = std::stoi(fm.catch_parameter("nz"));
    nb = 1;

    nPoints = nx * ny * nz;

    dx = std::stof(fm.catch_parameter("dx"));    
    dy = std::stof(fm.catch_parameter("dy"));    
    dz = std::stof(fm.catch_parameter("dz"));    

    value = new float[nPoints];

    std::string model_file = fm.catch_parameter("model_file");

    fm.read_binary_float(model_file, value, nPoints);
}

void Model::expand()
{
    int nzz = nz + 2*nb;
    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;

    float * expandedModel = new float[nzz*nxx*nyy]();

    // Centering
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                expandedModel[z + x*nzz + y*nxx*nzz] = value[(z - nb) + (x - nb)*nz + (y - nb)*nx*nz];
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
                expandedModel[z + x*nzz + y*nxx*nzz] = value[0 + (x - nb)*nz + (y - nb)*nx*nz];
                expandedModel[(nzz - z - 1) + x*nzz + y*nxx*nzz] = value[(nz - 1) + (x - nb)*nz + (y - nb)*nx*nz];
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
                expandedModel[z + x*nzz + y*nxx*nzz] = expandedModel[z + nb*nzz + y*nxx*nzz];
                expandedModel[z + (nxx - x - 1)*nzz + y*nxx*nzz] = expandedModel[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
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
                expandedModel[z + x*nzz + y*nxx*nzz] = expandedModel[z + x*nzz + nb*nxx*nzz];
                expandedModel[z + x*nzz + (nyy - y - 1)*nxx*nzz] = expandedModel[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }

    delete[] value;

    nz = nzz;
    nx = nxx;
    ny = nyy;

    nPoints = nxx*nyy*nzz;

    value = new float[nPoints]();

    std::swap(value, expandedModel);

    delete[] expandedModel;
}

void Model::reduce()
{
    int nzz = nz;
    int nxx = nx;
    
    nz -= 2*nb;
    nx -= 2*nb;
    ny -= 2*nb;
    
    nPoints = nx*ny*nz;

    float * reducedModel = new float[nPoints];

    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         // y direction
        int x = (int) (index - y*nx*nz) / nz;    // x direction
        int z = (int) (index - x*nz - y*nx*nz);  // z direction

        reducedModel[z + x*nz + y*nx*nz] = value[(z + nb) + (x + nb)*nzz + (y + nb)*nxx*nzz];
    }

    delete[] value;

    value = new float[nPoints]();
    
    std::swap(value, reducedModel);

    delete[] reducedModel;
}

void Model::smooth(float stdv, int samples)
{
    auto smooth = Gaussian();

    smooth.zdim = nz; 
    smooth.xdim = nx; 
    smooth.ydim = ny; 

    smooth.stdv = stdv;
    smooth.samples = samples;

    smooth.volume = new float[nPoints]();

    for (int index = 0; index < nPoints; index++)
        smooth.volume[index] = 1.0f / value[index];

    smooth.gaussian();

    for (int index = 0; index < nPoints; index++)
        value[index] = 1.0f / smooth.volume[index];

    delete[] smooth.volume;
}

void Model::resize(float new_dz, float new_dx, float new_dy)
{
    auto interp = Trilinear();

    interp.new_dz = new_dz;
    interp.new_dx = new_dx;
    interp.new_dy = new_dy;

    interp.new_nz = (int)(static_cast<float>(nz-1) * dz / interp.new_dz) + 1;    
    interp.new_nx = (int)(static_cast<float>(nx-1) * dx / interp.new_dx) + 1;    
    interp.new_ny = (int)(static_cast<float>(ny-1) * dy / interp.new_dy) + 1;    

    float * z = new float[interp.new_nz]();
    float * x = new float[interp.new_nx]();
    float * y = new float[interp.new_ny]();

    for (int i = 0; i < interp.new_nz; i++) 
        z[i] = static_cast<float>(i) * interp.new_dz;
    
    for (int j = 0; j < interp.new_nx; j++) 
        x[j] = static_cast<float>(j) * interp.new_dx;
    
    for (int k = 0; k < interp.new_ny; k++) 
        y[k] = static_cast<float>(k) * interp.new_dy;

    float * interpolatedModel = new float[interp.new_nz * interp.new_nx * interp.new_ny]();

    for (int k = 0; k < interp.new_ny; k++)
    {    
        for (int j = 0; j < interp.new_nx; j++)
        {
            for (int i = 0; i < interp.new_nz; i++)
            {
                if ((z[i] > 0.0f) && (z[i] < (nz-1)*dz) && 
                    (x[j] > 0.0f) && (x[j] < (nx-1)*dx) && 
                    (y[k] > 0.0f) && (y[k] < (ny-1)*dy))
                {
                    interp.z = z[i];
                    interp.x = x[j];
                    interp.y = y[k];

                    interp.z0 = floorf(z[i] / dz) * dz;
                    interp.x0 = floorf(x[j] / dx) * dx;
                    interp.y0 = floorf(y[k] / dy) * dy;

                    interp.z1 = floorf(z[i] / dz) * dz + dz;
                    interp.x1 = floorf(x[j] / dx) * dx + dx;
                    interp.y1 = floorf(y[k] / dy) * dy + dy;

                    int id = static_cast<int>(z[i]/dz) + static_cast<int>(x[j]/dx)*nz + static_cast<int>(y[k]/dy)*nx*nz;

                    interp.c000 = value[id];
                    interp.c001 = value[id + 1];
                    interp.c100 = value[id + nz]; 
                    interp.c101 = value[id + 1 + nz]; 
                    interp.c010 = value[id + nx*nz]; 
                    interp.c011 = value[id + 1 + nx*nz]; 
                    interp.c110 = value[id + nz + nx*nz]; 
                    interp.c111 = value[id + 1 + nz + nx*nz];

                    int interp_id = static_cast<int>(z[i]/interp.new_dz) + static_cast<int>(x[j]/interp.new_dx)*interp.new_nz + static_cast<int>(y[k]/interp.new_dy)*interp.new_nx*interp.new_nz;

                    interpolatedModel[interp_id] = interp.trilinear();
                }
            }
        }
    }

    int zfactor = static_cast<int>(floorf(dz / interp.new_dz)) + 1;
    int xfactor = static_cast<int>(floorf(dx / interp.new_dx)) + 1;
    int yfactor = static_cast<int>(floorf(dy / interp.new_dy)) + 1;

    for (int idz = 0; idz < zfactor; idz++)
    {
        for (int idy = yfactor; idy < interp.new_ny - yfactor; idy++)
        {
            for (int idx = xfactor; idx < interp.new_nx - xfactor; idx++)
            {
                interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[zfactor + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz];    
                interpolatedModel[(interp.new_nz - idz - 1) + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[(interp.new_nz - zfactor - 1) + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz];          
            }
        }
    }

    for (int idx = 0; idx < xfactor; idx++)
    {
        for (int idz = 0; idz < interp.new_nz; idz++)
        {
            for (int idy = yfactor; idy < interp.new_ny - yfactor; idy++)
            {
                interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + xfactor*interp.new_nz + idy*interp.new_nx*interp.new_nz];    
                interpolatedModel[idz + (interp.new_nx - idx - 1)*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + (interp.new_nx - xfactor - 1)*interp.new_nz + idy*interp.new_nx*interp.new_nz];          
            }
        }
    }

    for (int idy = 0; idy < yfactor; idy++)
    {
        for (int idz = 0; idz < interp.new_nz; idz++)
        {
            for (int idx = 0; idx < interp.new_nx; idx++)
            {
                interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + idx*interp.new_nz + yfactor*interp.new_nx*interp.new_nz];    
                interpolatedModel[idz + idx*interp.new_nz + (interp.new_ny - idy - 1)*interp.new_nx*interp.new_nz] = interpolatedModel[idz + idx*interp.new_nz + (interp.new_ny - yfactor - 1)*interp.new_nx*interp.new_nz];          
            }
        }
    }

    delete[] value;

    nx = interp.new_nx;
    ny = interp.new_ny;
    nz = interp.new_nz;

    dx = interp.new_dx;
    dy = interp.new_dy;
    dz = interp.new_dz;

    nPoints = nx * ny * nz;

    value = new float[nPoints];

    std::swap(value, interpolatedModel);

    delete[] interpolatedModel;
    
    interp.~Trilinear();
}


