# include <cmath>
# include <iostream>

# include "acoustic.hpp"

void Acoustic::sourceGenerator()
{
    source = new float[nt]();

    float pi = 4.0f * atanf(1.0f);
    float fc = fcut / (3.0f * sqrtf(pi));

    for (int t = 0; t < nt; t++)
    {
        float aux1 = 1.0f - 2.0f * pi * powf(t*dt - tlag, 2.0f) * powf(fc, 2.0f) * pow(pi, 2.0f);
        float aux2 = expf(-pi * powf(t*dt - tlag, 2.0f)*powf(fc, 2.0f)*powf(pi, 2.0f));

        source[t] = aux1 * aux2;
    }
}

void Acoustic::dampingGenerator()
{
    damp1D = new float[nb]();
    damp2D = new float[nb * nb]();
    damp3D = new float[nb * nb * nb]();

    /* 1D damp construction */
    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = expf(-powf(factor * (nb - i), 2.0f)); 
    }

    /* 2D damp construction */
    for(int i = 0; i < nb; i++) 
    {
        damp2D[i + i*nb] = damp1D[i];            

        for (int j = i; j < nb; j++)
        {   
            damp2D[j + i*nb] = damp1D[i];
            damp2D[i + j*nb] = damp1D[i];
        }    
    }

    /* 3D damp construction */
    for (int i  = 0; i < nb; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            damp3D[i + j*nb + j*nb*nb] = damp2D[i + j*nb];
        }

        for(int j = 0; j < nb-1; j++)
        {
            for(int k = j+1; k < nb; k++)
            {
                damp3D[i + j*nb + k*nb*nb] = damp3D[i + j*nb + j*nb*nb];
                damp3D[i + k*nb + j*nb*nb] = damp3D[i + j*nb + j*nb*nb];
            }
        }
    }        
}

void Acoustic::acousticWaveSolver()
{


}

void Acoustic::dampingApplicator()
{


}

void Acoustic::wavefieldUpdate()
{


}

void Acoustic::getSeismogram()
{


}

void Acoustic::exportSeismogram()
{


}
