# include "eikonal_model.hpp"

float * Eikonal_model::expand_fdm(float * volume)
{
    x_samples_b = x_samples + 2;    
    y_samples_b = y_samples + 2;    
    z_samples_b = z_samples + 2;    
    
    total_samples_b = x_samples_b * y_samples_b * z_samples_b;

    float * new_volume = new float[total_samples_b]();

    // Centering
    for (int z = 1; z < z_samples_b - 1; z++)
    {
        for (int y = 1; y < y_samples_b - 1; y++)
        {
            for (int x = 1; x < x_samples_b - 1; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[(z - 1) + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
            }
        }
    }

    // Z direction
    for (int z = 0; z < 1; z++)
    {
        for (int y = 1; y < y_samples_b - 1; y++)
        {
            for (int x = 1; x < x_samples_b - 1; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[0 + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
                new_volume[(z_samples_b - z - 1) + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[(z_samples - 1) + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
            }
        }
    }

    // X direction
    for (int x = 0; x < 1; x++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int y = 1; y < y_samples_b - 1; y++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + 1*z_samples_b + y*x_samples_b*z_samples_b];
                new_volume[z + (x_samples_b - x - 1)*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + (x_samples_b - 1 - 1)*z_samples_b + y*x_samples_b*z_samples_b];
            }
        }
    }

    // Y direction
    for (int y = 0; y < 1; y++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int x = 0; x < x_samples_b; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + x*z_samples_b + 1*x_samples_b*z_samples_b];
                new_volume[z + x*z_samples_b + (y_samples_b - y - 1)*x_samples_b*z_samples_b] = new_volume[z + x*z_samples_b + (y_samples_b - 1 - 1)*x_samples_b*z_samples_b];
            }
        }
    }

    return new_volume;
}

float * Eikonal_model::reduce_fdm(float * volume)
{    
    float * new_volume = new float[total_samples];

    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        new_volume[z + x*z_samples + y*x_samples*z_samples] = volume[(z + 1) + (x + 1)*z_samples_b + (y + 1)*x_samples_b*z_samples_b];
    }
    
    return new_volume;
}

float * Eikonal_model::expand_pad(float * volume, int padx, int pady, int padz)
{
    x_samples_b = x_samples + padx;    
    y_samples_b = y_samples + pady;    
    z_samples_b = z_samples + padz;    
    
    total_samples_b = x_samples_b * y_samples_b * z_samples_b;

    float * new_volume = new float[total_samples_b]();

    // Centering
    for (int z = 0; z < z_samples; z++)
    {
        for (int y = 0; y < y_samples; y++)
        {
            for (int x = 0; x < x_samples; x++)
            {
                new_volume[z + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[z + x*z_samples + y*x_samples*z_samples];
            }
        }
    }

    // Z direction
    for (int z = 0; z < padz; z++)
    {
        for (int y = 0; y < y_samples_b - pady; y++)
        {
            for (int x = 0; x < x_samples_b - padx; x++)
            {
                new_volume[(z_samples_b - z - 1) + x*z_samples_b + y*x_samples_b*z_samples_b] = volume[(z_samples - 1) + x*z_samples + y*x_samples*z_samples];
            }
        }
    }

    // X direction
    for (int x = 0; x < padx; x++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int y = 0; y < y_samples_b - pady; y++)
            {
                new_volume[z + (x_samples_b - x - 1)*z_samples_b + y*x_samples_b*z_samples_b] = new_volume[z + (x_samples_b - padx - 1)*z_samples_b + y*x_samples_b*z_samples_b];
            }
        }
    }

    // Y direction
    for (int y = 0; y < pady; y++)
    {
        for (int z = 0; z < z_samples_b; z++)
        {
            for (int x = 0; x < x_samples_b; x++)
            {
                new_volume[z + x*z_samples_b + (y_samples_b - y - 1)*x_samples_b*z_samples_b] = new_volume[z + x*z_samples_b + (y_samples_b - pady - 1)*x_samples_b*z_samples_b];
            }
        }
    }

    return new_volume;
}

float * Eikonal_model::reduce_pad(float * volume, int padx, int pady, int padz)
{
    float * new_volume = new float[total_samples];

    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        new_volume[z + x*z_samples + y*x_samples*z_samples] = volume[z + x*z_samples_b + y*x_samples_b*z_samples_b];
    }
    
    return new_volume;
}