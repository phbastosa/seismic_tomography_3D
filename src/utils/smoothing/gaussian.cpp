# include <cmath>

# include "gaussian.hpp"

void Gaussian::gaussian()
{
    int init = samples / 2;
    int nPoints = xdim * ydim * zdim;
    int nKernel = samples * samples * samples;

    float pi = 4.0f * atanf(1.0f); 

    float * kernel = new float[nKernel]();
    float * smoothed = new float[nPoints]();

    for (int i = 0; i < nPoints; i++) 
        smoothed[i] = volume[i];

    int mid = (int)(samples / 2); 

    kernel[mid + mid*samples + mid*samples*samples] = 1.0f;

    if (stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int y = -init; y <= init; y++)
        {
            for (int x = -init; x <= init; x++)
            {
                for (int z = -init; z <= init; z++)
                {          
                    int index = (z+init) + (x+init)*samples + (y+init)*samples*samples; 

                    kernel[index] = 1.0f / (2.0f*pi*stdv*stdv) * expf(-((x*x + y*y + z*z) / (2.0f*stdv*stdv)));
        
                    sum += kernel[index]; 
                }
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
        
    for (int k = init; k < ydim - init; k++)
    {   
        for (int j = init; j < xdim - init; j++)
        {
            for (int i = init; i < zdim - init; i++)
            {       
                float accum = 0.0f;
                
                for (int yk = 0; yk < samples; yk++)
                {      
                    for (int xk = 0; xk < samples; xk++)
                    {      
                        for (int zk = 0; zk < samples; zk++)
                        {   
                            int index = zk + xk*samples + yk*samples*samples;   
                            int partial = (i-init+zk) + (j-init+xk)*zdim + (k-init+yk)*xdim*zdim; 

                            accum += volume[partial] * kernel[index];
                        }        
                    }
                }
                
                smoothed[i + j*zdim + k*xdim*zdim] = accum;
            }
        }   
    }

    std::swap(volume, smoothed);

    delete[] kernel;
    delete[] smoothed;
}