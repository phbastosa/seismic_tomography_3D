# include "eikonal_model.hpp"

void Eikonal_model::expand_fdm(float * input, float * output)
{
    // Centering
    for (int z = 1; z < z_samples_b - 1; z++)
    {
        for (int y = 1; y < y_samples_b - 1; y++)
        {
            for (int x = 1; x < x_samples_b - 1; x++)
            {
                output[z + x*z_samples_b + y*x_samples_b*z_samples_b] = input[(z - 1) + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
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
                output[z + x*z_samples_b + y*x_samples_b*z_samples_b] = input[0 + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
                output[(z_samples_b - z - 1) + x*z_samples_b + y*x_samples_b*z_samples_b] = input[(z_samples - 1) + (x - 1)*z_samples + (y - 1)*x_samples*z_samples];
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
                output[z + x*z_samples_b + y*x_samples_b*z_samples_b] = output[z + 1*z_samples_b + y*x_samples_b*z_samples_b];
                output[z + (x_samples_b - x - 1)*z_samples_b + y*x_samples_b*z_samples_b] = output[z + (x_samples_b - 1 - 1)*z_samples_b + y*x_samples_b*z_samples_b];
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
                output[z + x*z_samples_b + y*x_samples_b*z_samples_b] = output[z + x*z_samples_b + 1*x_samples_b*z_samples_b];
                output[z + x*z_samples_b + (y_samples_b - y - 1)*x_samples_b*z_samples_b] = output[z + x*z_samples_b + (y_samples_b - 1 - 1)*x_samples_b*z_samples_b];
            }
        }
    }
}

void Eikonal_model::reduce_fdm(float * input, float * output)
{    
    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        output[z + x*z_samples + y*x_samples*z_samples] = input[(z + 1) + (x + 1)*z_samples_b + (y + 1)*x_samples_b*z_samples_b];
    }
}

void Eikonal_model::expand_pad(float * input, float * output, int padx, int pady, int padz)
{
    // Centering
    for (int z = 0; z < z_samples; z++)
    {
        for (int y = 0; y < y_samples; y++)
        {
            for (int x = 0; x < x_samples; x++)
            {
                output[z + x*z_samples_b + y*x_samples_b*z_samples_b] = input[z + x*z_samples + y*x_samples*z_samples];
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
                output[(z_samples_b - z - 1) + x*z_samples_b + y*x_samples_b*z_samples_b] = input[(z_samples - 1) + x*z_samples + y*x_samples*z_samples];
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
                output[z + (x_samples_b - x - 1)*z_samples_b + y*x_samples_b*z_samples_b] = output[z + (x_samples_b - padx - 1)*z_samples_b + y*x_samples_b*z_samples_b];
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
                output[z + x*z_samples_b + (y_samples_b - y - 1)*x_samples_b*z_samples_b] = output[z + x*z_samples_b + (y_samples_b - pady - 1)*x_samples_b*z_samples_b];
            }
        }
    }
}

void Eikonal_model::reduce_pad(float * input, float * output, int padx, int pady, int padz)
{
    for (int index = 0; index < total_samples; index++)
    {
        int y = (int) (index / (x_samples*z_samples));         
        int x = (int) (index - y*x_samples*z_samples) / z_samples;    
        int z = (int) (index - x*z_samples - y*x_samples*z_samples);  

        output[z + x*z_samples + y*x_samples*z_samples] = input[z + x*z_samples_b + y*x_samples_b*z_samples_b];
    }
}