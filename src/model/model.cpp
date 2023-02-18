# include <cmath> 

# include "model.hpp"

# include "../utils/smoothing/gaussian.hpp"
# include "../utils/interpolation/trilinear.hpp"

float * Model::expand(float * volume)
{
    float * new_volume = new float[total_samples_b]();

    // Centering
    for (int z = boundary_samples; z < z_samples_b - boundary_samples; z++)
    {
        for (int y = boundary_samples; y < y_samples_b - boundary_samples; y++)
        {
            for (int x = boundary_samples; x < x_samples_b - boundary_samples; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[(z - boundary_samples) + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
            }
        }
    }

    // Z direction
    for (int z = 0; z < boundary_samples; z++)
    {
        for (int y = boundary_samples; y < y_samples_b - boundary_samples; y++)
        {
            for (int x = boundary_samples; x < x_samples_b - boundary_samples; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[0 + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
                new_volume[(z_samples_b - z - 1) + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[(z_samples - 1) + (x - boundary_samples)*z_samples + (y - boundary_samples)*x_samples*z_samples];
            }
        }
    }

    // X direction
    for (int x = 0; x < boundary_samples; x++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int y = boundary_samples; y < y_samples_b - boundary_samples; y++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + boundary_samples*z_samples_b + y*x_samples_b*z_samples_b];
                new_volume[z + (x_samples_b - x - 1)*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + (x_samples_b - boundary_samples - 1)*z_samples_b + y*x_samples_b*z_samples_b];
            }
        }
    }

    // Y direction
    for (int y = 0; y < boundary_samples; y++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int x = 0; x < x_samples_b; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + x*z_samples_b + boundary_samples*x_samples_b*z_samples_b];
                new_volume[z + x*z_samples_b + (y_samples_b - y - 1)*x_samples_b*z_samples_b] = new_volume[z + x*z_samples_b + (y_samples_b - boundary_samples - 1)*x_samples_b*z_samples_b];
            }
        }
    }

    return new_volume;
}

float * Model::reduce(float * volume)
{    
    float * new_volume = new float[total_samples];

    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        new_volume[z + x*z_samples + y*x_samples*z_samples] = volume[(z + boundary_samples) + (x + boundary_samples)*z_samples_b + (y + boundary_samples)*x_samples_b*z_samples_b];
    }
    
    return new_volume;
}

float * Model::smooth(float * volume)
{
    auto smoother = Gaussian();

    smoother.zdim = z_samples; 
    smoother.xdim = x_samples; 
    smoother.ydim = y_samples; 

    smoother.stdv = standard_deviation;
    smoother.samples = smoothing_samples;

    smoother.volume = new float[total_samples]();

    for (int index = 0; index < total_samples; index++)
        smoother.volume[index] = 1.0f / volume[index];

    smoother.gaussian();

    for (int index = 0; index < total_samples; index++)
        volume[index] = 1.0f / smoother.volume[index];

    return smoother.volume;
}

float * Model::resize(float * volume)
{
    auto interpolate = Trilinear();

    new_z_samples = (int)(static_cast<float>(z_samples-1) * z_spacing / new_z_spacing) + 1;    
    new_x_samples = (int)(static_cast<float>(x_samples-1) * x_spacing / new_x_spacing) + 1;    
    new_y_samples = (int)(static_cast<float>(y_samples-1) * y_spacing / new_y_spacing) + 1;    

    float * z = new float[new_z_samples]();
    float * x = new float[new_x_samples]();
    float * y = new float[new_y_samples]();

    for (int i = 0; i < new_z_samples; i++) 
        z[i] = static_cast<float>(i) * new_z_spacing;
    
    for (int j = 0; j < new_x_samples; j++) 
        x[j] = static_cast<float>(j) * new_x_spacing;
    
    for (int k = 0; k < new_y_samples; k++) 
        y[k] = static_cast<float>(k) * new_y_spacing;

    float * new_volume = new float[new_z_samples * new_x_samples * new_y_samples]();

    for (int k = 0; k < new_y_samples; k++)
    {    
        for (int j = 0; j < new_x_samples; j++)
        {
            for (int i = 0; i < new_z_samples; i++)
            {
                if ((z[i] > 0.0f) && (z[i] < (z_samples-1)*z_spacing) && 
                    (x[j] > 0.0f) && (x[j] < (x_samples-1)*x_spacing) && 
                    (y[k] > 0.0f) && (y[k] < (y_samples-1)*y_spacing))
                {
                    interpolate.z = z[i];
                    interpolate.x = x[j];
                    interpolate.y = y[k];

                    interpolate.z0 = floorf(z[i] / z_spacing) * z_spacing;
                    interpolate.x0 = floorf(x[j] / x_spacing) * x_spacing;
                    interpolate.y0 = floorf(y[k] / y_spacing) * y_spacing;

                    interpolate.z1 = floorf(z[i] / z_spacing) * z_spacing + z_spacing;
                    interpolate.x1 = floorf(x[j] / x_spacing) * x_spacing + x_spacing;
                    interpolate.y1 = floorf(y[k] / y_spacing) * y_spacing + y_spacing;

                    int id = static_cast<int>(z[i]/z_spacing) + 
                             static_cast<int>(x[j]/x_spacing)*z_samples + 
                             static_cast<int>(y[k]/y_spacing)*x_samples*z_samples;

                    interpolate.c000 = volume[id];
                    interpolate.c001 = volume[id + 1];
                    interpolate.c100 = volume[id + z_samples]; 
                    interpolate.c101 = volume[id + 1 + z_samples]; 
                    interpolate.c010 = volume[id + x_samples*z_samples]; 
                    interpolate.c011 = volume[id + 1 + x_samples*z_samples]; 
                    interpolate.c110 = volume[id + z_samples + x_samples*z_samples]; 
                    interpolate.c111 = volume[id + 1 + z_samples + x_samples*z_samples];

                    int interp_id = static_cast<int>(z[i]/new_z_spacing) + 
                                    static_cast<int>(x[j]/new_x_spacing)*new_z_samples + 
                                    static_cast<int>(y[k]/new_y_spacing)*new_x_samples*new_z_samples;

                    new_volume[interp_id] = interpolate.trilinear();
                }
            }
        }
    }

    int zfactor = static_cast<int>(floorf(z_spacing / new_z_spacing)) + 1;
    int xfactor = static_cast<int>(floorf(x_spacing / new_x_spacing)) + 1;
    int yfactor = static_cast<int>(floorf(y_spacing / new_y_spacing)) + 1;

    for (int idz = 0; idz < zfactor; idz++)
    {
        for (int idy = yfactor; idy < new_y_samples - yfactor; idy++)
        {
            for (int idx = xfactor; idx < new_x_samples - xfactor; idx++)
            {
                new_volume[idz + idx*new_z_samples + idy*new_x_samples*new_z_samples] = new_volume[zfactor + idx*new_z_samples + idy*new_x_samples*new_z_samples];    
                new_volume[(new_z_samples - idz - 1) + idx*new_z_samples + idy*new_x_samples*new_z_samples] = new_volume[(new_z_samples - zfactor - 1) + idx*new_z_samples + idy*new_x_samples*new_z_samples];          
            }
        }
    }

    for (int idx = 0; idx < xfactor; idx++)
    {
        for (int idz = 0; idz < new_z_samples; idz++)
        {
            for (int idy = yfactor; idy < new_y_samples - yfactor; idy++)
            {
                new_volume[idz + idx*new_z_samples + idy*new_x_samples*new_z_samples] = new_volume[idz + xfactor*new_z_samples + idy*new_x_samples*new_z_samples];    
                new_volume[idz + (new_x_samples - idx - 1)*new_z_samples + idy*new_x_samples*new_z_samples] = new_volume[idz + (new_x_samples - xfactor - 1)*new_z_samples + idy*new_x_samples*new_z_samples];          
            }
        }
    }

    for (int idy = 0; idy < yfactor; idy++)
    {
        for (int idz = 0; idz < new_z_samples; idz++)
        {
            for (int idx = 0; idx < new_x_samples; idx++)
            {
                new_volume[idz + idx*new_z_samples + idy*new_x_samples*new_z_samples] = new_volume[idz + idx*new_z_samples + yfactor*new_x_samples*new_z_samples];    
                new_volume[idz + idx*new_z_samples + (new_y_samples - idy - 1)*new_x_samples*new_z_samples] = new_volume[idz + idx*new_z_samples + (new_y_samples - yfactor - 1)*new_x_samples*new_z_samples];          
            }
        }
    }
    
    interpolate.~Trilinear();

    return new_volume;
}


