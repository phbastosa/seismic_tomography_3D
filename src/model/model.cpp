# include "model.hpp"

# include "../utils/file_manager/file_manager.hpp"

void Model::set_parameters(std::string file)
{
    auto fm = File_manager();

    fm.parameter_file = file;

    x_samples = std::stoi(fm.catch_parameter("x_samples"));
    y_samples = std::stoi(fm.catch_parameter("y_samples"));
    z_samples = std::stoi(fm.catch_parameter("z_samples"));
    
    total_samples = x_samples * y_samples * z_samples;

    x_spacing = std::stof(fm.catch_parameter("x_spacing"));    
    y_spacing = std::stof(fm.catch_parameter("y_spacing"));    
    z_spacing = std::stof(fm.catch_parameter("z_spacing"));    

    property = new float[total_samples];
}

void Model::expand()
{
    int nzz = z_samples + 2*boundary_samples;
    int nxx = x_samples + 2*boundary_samples;
    int nyy = y_samples + 2*boundary_samples;

    float * expandedModel = new float[nzz*nxx*nyy]();

    // Centering
    for (int z = boundary_samples; z < nzz - boundary_samples; z++)
    {
        for (int y = boundary_samples; y < nyy - boundary_samples; y++)
        {
            for (int x = boundary_samples; x < nxx - boundary_samples; x++)
            {
                expandedModel[z + x*nzz + y*nxx*nzz] = property[(z - boundary_samples) + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
            }
        }
    }

    // Z direction
    for (int z = 0; z < boundary_samples; z++)
    {
        for (int y = boundary_samples; y < nyy - boundary_samples; y++)
        {
            for (int x = boundary_samples; x < nxx - boundary_samples; x++)
            {
                expandedModel[z + x*nzz + y*nxx*nzz] = property[0 + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
                expandedModel[(nzz - z - 1) + x*nzz + y*nxx*nzz] = property[(z_samples - 1) + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
            }
        }
    }

    // X direction
    for (int x = 0; x < boundary_samples; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = boundary_samples; y < nyy - boundary_samples; y++)
            {
                expandedModel[z + x*nzz + y*nxx*nzz] = expandedModel[z + boundary_samples*nzz + y*nxx*nzz];
                expandedModel[z + (nxx - x - 1)*nzz + y*nxx*nzz] = expandedModel[z + (nxx - boundary_samples - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < boundary_samples; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                expandedModel[z + x*nzz + y*nxx*nzz] = expandedModel[z + x*nzz + boundary_samples*nxx*nzz];
                expandedModel[z + x*nzz + (nyy - y - 1)*nxx*nzz] = expandedModel[z + x*nzz + (nyy - boundary_samples - 1)*nxx*nzz];
            }
        }
    }

    delete[] property;

    z_samples = nzz;
    x_samples = nxx;
    y_samples = nyy;

    total_samples = nxx*nyy*nzz;

    property = new float[total_samples]();

    std::swap(property, expandedModel);

    delete[] expandedModel;
}

void Model::reduce()
{
    int nzz = z_samples;
    int nxx = x_samples;
    
    z_samples -= 2*boundary_samples;
    x_samples -= 2*boundary_samples;
    y_samples -= 2*boundary_samples;
    
    total_samples = x_samples*y_samples*z_samples;

    float * reducedModel = new float[total_samples];

    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        reducedModel[z + x*z_samples + y*x_samples*z_samples] = property[(z + boundary_samples) + (x + boundary_samples)*nzz + (y + boundary_samples)*nxx*nzz];
    }

    delete[] property;

    property = new float[total_samples]();
    
    std::swap(property, reducedModel);

    delete[] reducedModel;
}

// void Model::smooth(float stdv, int samples)
// {
//     auto smooth = Gaussian();

//     smooth.zdim = z_samples; 
//     smooth.xdim = x_samples; 
//     smooth.ydim = y_samples; 

//     smooth.stdv = stdv;
//     smooth.samples = samples;

//     smooth.volume = new float[total_samples]();

//     for (int index = 0; index < total_samples; index++)
//         smooth.volume[index] = 1.0f / property[index];

//     smooth.gaussian();

//     for (int index = 0; index < total_samples; index++)
//         property[index] = 1.0f / smooth.volume[index];

//     delete[] smooth.volume;
// }

// void Model::resize(float new_dz, float new_dx, float new_dy)
// {
//     auto interp = Trilinear();

//     interp.new_dz = new_dz;
//     interp.new_dx = new_dx;
//     interp.new_dy = new_dy;

//     interp.new_nz = (int)(static_cast<float>(z_samples-1) * dz / interp.new_dz) + 1;    
//     interp.new_nx = (int)(static_cast<float>(x_samples-1) * dx / interp.new_dx) + 1;    
//     interp.new_ny = (int)(static_cast<float>(y_samples-1) * dy / interp.new_dy) + 1;    

//     float * z = new float[interp.new_nz]();
//     float * x = new float[interp.new_nx]();
//     float * y = new float[interp.new_ny]();

//     for (int i = 0; i < interp.new_nz; i++) 
//         z[i] = static_cast<float>(i) * interp.new_dz;
    
//     for (int j = 0; j < interp.new_nx; j++) 
//         x[j] = static_cast<float>(j) * interp.new_dx;
    
//     for (int k = 0; k < interp.new_ny; k++) 
//         y[k] = static_cast<float>(k) * interp.new_dy;

//     float * interpolatedModel = new float[interp.new_nz * interp.new_nx * interp.new_ny]();

//     for (int k = 0; k < interp.new_ny; k++)
//     {    
//         for (int j = 0; j < interp.new_nx; j++)
//         {
//             for (int i = 0; i < interp.new_nz; i++)
//             {
//                 if ((z[i] > 0.0f) && (z[i] < (z_samples-1)*dz) && 
//                     (x[j] > 0.0f) && (x[j] < (x_samples-1)*dx) && 
//                     (y[k] > 0.0f) && (y[k] < (y_samples-1)*dy))
//                 {
//                     interp.z = z[i];
//                     interp.x = x[j];
//                     interp.y = y[k];

//                     interp.z0 = floorf(z[i] / dz) * dz;
//                     interp.x0 = floorf(x[j] / dx) * dx;
//                     interp.y0 = floorf(y[k] / dy) * dy;

//                     interp.z1 = floorf(z[i] / dz) * dz + dz;
//                     interp.x1 = floorf(x[j] / dx) * dx + dx;
//                     interp.y1 = floorf(y[k] / dy) * dy + dy;

//                     int id = static_cast<int>(z[i]/dz) + 
//                              static_cast<int>(x[j]/dx)*z_samples + 
//                              static_cast<int>(y[k]/dy)*x_samples*z_samples;

//                     interp.c000 = property[id];
//                     interp.c001 = property[id + 1];
//                     interp.c100 = property[id + z_samples]; 
//                     interp.c101 = property[id + 1 + z_samples]; 
//                     interp.c010 = property[id + x_samples*z_samples]; 
//                     interp.c011 = property[id + 1 + x_samples*z_samples]; 
//                     interp.c110 = property[id + z_samples + x_samples*z_samples]; 
//                     interp.c111 = property[id + 1 + z_samples + x_samples*z_samples];

//                     int interp_id = static_cast<int>(z[i]/interp.new_dz) + 
//                                     static_cast<int>(x[j]/interp.new_dx)*interp.new_nz + 
//                                     static_cast<int>(y[k]/interp.new_dy)*interp.new_nx*interp.new_nz;

//                     interpolatedModel[interp_id] = interp.trilinear();
//                 }
//             }
//         }
//     }

//     int zfactor = static_cast<int>(floorf(dz / interp.new_dz)) + 1;
//     int xfactor = static_cast<int>(floorf(dx / interp.new_dx)) + 1;
//     int yfactor = static_cast<int>(floorf(dy / interp.new_dy)) + 1;

//     for (int idz = 0; idz < zfactor; idz++)
//     {
//         for (int idy = yfactor; idy < interp.new_ny - yfactor; idy++)
//         {
//             for (int idx = xfactor; idx < interp.new_nx - xfactor; idx++)
//             {
//                 interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[zfactor + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz];    
//                 interpolatedModel[(interp.new_nz - idz - 1) + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[(interp.new_nz - zfactor - 1) + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz];          
//             }
//         }
//     }

//     for (int idx = 0; idx < xfactor; idx++)
//     {
//         for (int idz = 0; idz < interp.new_nz; idz++)
//         {
//             for (int idy = yfactor; idy < interp.new_ny - yfactor; idy++)
//             {
//                 interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + xfactor*interp.new_nz + idy*interp.new_nx*interp.new_nz];    
//                 interpolatedModel[idz + (interp.new_nx - idx - 1)*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + (interp.new_nx - xfactor - 1)*interp.new_nz + idy*interp.new_nx*interp.new_nz];          
//             }
//         }
//     }

//     for (int idy = 0; idy < yfactor; idy++)
//     {
//         for (int idz = 0; idz < interp.new_nz; idz++)
//         {
//             for (int idx = 0; idx < interp.new_nx; idx++)
//             {
//                 interpolatedModel[idz + idx*interp.new_nz + idy*interp.new_nx*interp.new_nz] = interpolatedModel[idz + idx*interp.new_nz + yfactor*interp.new_nx*interp.new_nz];    
//                 interpolatedModel[idz + idx*interp.new_nz + (interp.new_ny - idy - 1)*interp.new_nx*interp.new_nz] = interpolatedModel[idz + idx*interp.new_nz + (interp.new_ny - yfactor - 1)*interp.new_nx*interp.new_nz];          
//             }
//         }
//     }

//     delete[] property;

//     x_samples = interp.new_nx;
//     y_samples = interp.new_ny;
//     z_samples = interp.new_nz;

//     dx = interp.new_dx;
//     dy = interp.new_dy;
//     dz = interp.new_dz;

//     total_samples = x_samples * y_samples * z_samples;

//     property = new float[total_samples];

//     std::swap(property, interpolatedModel);

//     delete[] interpolatedModel;
    
//     interp.~Trilinear();
// }


