# ifndef WAVE_FUNCTIONS_HPP
# define WAVE_FUNCTIONS_HPP

# include <cmath>
# include <iostream>

void progressMessage(int timeStep, int nt, float dt, int sId, float sx, float sy, float sz)
{
    if (timeStep % (nt/10) == 0)
    {    
        system("clear");
        std::cout<<"Running shot "<<sId+1<<" at position: (z = "<<sz<<", x = "<<sx<<", y = "<<sy<<")"<<std::endl;
        std::cout<<"Time step: "<<timeStep*dt<<std::endl;
    }
}

void pml_dampers(float * damp1D, float * damp2D, float * damp3D, float factor, int nb)
{
    float pi = 4.0f * atan(1.0f);

    /* 1D damp construction */
    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = factor * (1.0 - cos(pi * (nb - i) / (2.0 * nb)));
    }

    /* 2D damp construction */
    for(int i = 0; i < nb; i++) 
    {
        for (int j = 0; j < nb; j++)
        {   
            damp2D[j + i*nb] += damp1D[i]; // up to bottom
            damp2D[i + j*nb] += damp1D[i]; // left to right
        }
    }

    /* 3D damp construction */
    for (int i  = 0; i < nb; i++)
    {
        for(int j = 0; j < nb; j++)
        {
            for(int k = 0; k < nb; k++)
            {
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + j*nb]; // XY plane
                damp3D[i + j*nb + k*nb*nb] += damp2D[j + k*nb]; // ZX plane
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + k*nb]; // ZY plane
            }
        }
    }        
}

void setWaveField(float * U_pas, float * U_pre, float * U_fut, int nPoints)
{
    // pragma
    for (int index = 0; index < nPoints; index++)
    {
        U_pas[index] = 0.0f;
        U_pre[index] = 0.0f;
        U_fut[index] = 0.0f;
    }
}

void applyWavelet(float * U_pre, float * wavelet, int timeStep, int nsrc, int sId)
{
    # pragma acc kernels 
    {
        if (timeStep < nsrc) 
            U_pre[sId] += wavelet[timeStep];
    }
}

void pml_wavePropagation(float * V, float * U_pas, float * U_pre, float * U_fut, float * damp1D, float * damp2D, float * damp3D, int nb, int nxx, int nyy, int nzz, float dx, float dy, float dz, float dt)
{
    int nPoints = nxx*nyy*nzz;

    float alpha1, alpha2;
    float d2_Px2, d2_Py2, d2_Pz2;

    # pragma acc parallel loop present(V[0:nPoints],U_pas[0:nPoints],U_pre[0:nPoints],U_fut[0:nPoints],damp1D[0:nb],damp2D[0:nb*nb],damp3D[0:nb*nb*nb])
    for (int index = 0; index < nPoints; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        float damper = 0.0f;

        if((i >= 4) && (i < nzz-4) && (j >= 4) && (j < nxx-4) && (k >= 4) && (k < nyy-4)) 
        {
            d2_Px2 = (- 9.0f*(U_pre[i + (j-4)*nzz + k*nxx*nzz] + U_pre[i + (j+4)*nzz + k*nxx*nzz])
                  +   128.0f*(U_pre[i + (j-3)*nzz + k*nxx*nzz] + U_pre[i + (j+3)*nzz + k*nxx*nzz])
                  -  1008.0f*(U_pre[i + (j-2)*nzz + k*nxx*nzz] + U_pre[i + (j+2)*nzz + k*nxx*nzz])
                  +  8064.0f*(U_pre[i + (j-1)*nzz + k*nxx*nzz] + U_pre[i + (j+1)*nzz + k*nxx*nzz])
                  - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

            d2_Py2 = (- 9.0f*(U_pre[i + j*nzz + (k-4)*nxx*nzz] + U_pre[i + j*nzz + (k+4)*nxx*nzz])
                  +   128.0f*(U_pre[i + j*nzz + (k-3)*nxx*nzz] + U_pre[i + j*nzz + (k+3)*nxx*nzz])
                  -  1008.0f*(U_pre[i + j*nzz + (k-2)*nxx*nzz] + U_pre[i + j*nzz + (k+2)*nxx*nzz])
                  +  8064.0f*(U_pre[i + j*nzz + (k-1)*nxx*nzz] + U_pre[i + j*nzz + (k+1)*nxx*nzz])
                  - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

            d2_Pz2 = (- 9.0f*(U_pre[(i-4) + j*nzz + k*nxx*nzz] + U_pre[(i+4) + j*nzz + k*nxx*nzz])
                  +   128.0f*(U_pre[(i-3) + j*nzz + k*nxx*nzz] + U_pre[(i+3) + j*nzz + k*nxx*nzz])
                  -  1008.0f*(U_pre[(i-2) + j*nzz + k*nxx*nzz] + U_pre[(i+2) + j*nzz + k*nxx*nzz])
                  +  8064.0f*(U_pre[(i-1) + j*nzz + k*nxx*nzz] + U_pre[(i+1) + j*nzz + k*nxx*nzz])
                  - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));
        
            // 1D damping
            if((i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
            {
                damper = damp1D[i];
            }         
            else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
            {
                damper = damp1D[nb-(i-(nzz-nb))-1];
            }         
            else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
            {
                damper = damp1D[j];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb)) 
            {
                damper = damp1D[nb-(j-(nxx-nb))-1];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
            {
                damper = damp1D[k];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy)) 
            {
                damper = damp1D[nb-(k-(nyy-nb))-1];
            }

            // 2D damping 
            else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
            {
                damper = damp2D[j + k*nb];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
            {
                damper = damp2D[nb-(j-(nxx-nb))-1 + k*nb];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp2D[j + (nb-(k-(nyy-nb))-1)*nb];
            }
            else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp2D[nb-(j-(nxx-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
            }


            else if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
            {
                damper = damp2D[i + k*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
            {
                damper = damp2D[nb-(i-(nzz-nb))-1 + k*nb];
            }
            else if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp2D[i + (nb-(k-(nyy-nb))-1)*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp2D[nb-(i-(nzz-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
            }


            else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
            {
                damper = damp2D[i + j*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
            {
                damper = damp2D[nb-(i-(nzz-nb))-1 + j*nb];
            }
            else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
            {
                damper = damp2D[i + (nb-(j-(nxx-nb))-1)*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
            {
                damper = damp2D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb];
            }

            // 3D damping
            else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
            {
                damper = damp3D[i + j*nb + k*nb*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
            {
                damper = damp3D[nb-(i-(nzz-nb))-1 + j*nb + k*nb*nb];
            }
            else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
            {
                damper = damp3D[i + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
            }
            else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp3D[i + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
            {
                damper = damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp3D[nb-(i-(nzz-nb))-1 + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
            }
            else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp3D[i + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
            }
            else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
            {
                damper = damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
            }


            alpha1 = (damper*dt - 1.0)*U_pas[index];
            alpha2 = (2.0 - damper*damper*dt*dt) * U_pre[index];

            U_fut[index] = 1.0/(1.0 + damper*dt)*(alpha2 + alpha1 + V[index]*V[index]*dt*dt*(d2_Pz2 + d2_Px2 + d2_Py2));
        } 
    }
}    

void dampApplication(float * U_pre, float * U_fut, float * damp1D, float * damp2D, float * damp3D, int nxx, int nyy, int nzz, int nb)
{
    int nPoints = nxx*nyy*nzz;

    # pragma acc parallel loop present(U_pre[0:nPoints],U_fut[0:nPoints],damp1D[0:nb],damp2D[0:nb*nb],damp3D[0:nb*nb*nb])
    for (int index = 0; index < nPoints; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        // Damping 1D ------------------------------------------------------------------------

        if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
        {    
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[i];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[i];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp1D[i];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp1D[i];            
        }

        if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[j];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[j];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp1D[j];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp1D[j];            
        }

        if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[k];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[k];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp1D[k];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp1D[k];            
        }

        // Damping 2D ------------------------------------------------------------------------

        if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];            
        }

        if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];            
        }

        if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            
        }

        // Damping 3D ------------------------------------------------------------------------
        if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
        }
    }
}

void wavefieldUpdate(float * U_pas, float * U_pre, float * U_fut, int nPoints)
{
    # pragma acc parallel loop present(U_pas[0:nPoints], U_pre[0:nPoints], U_fut[0:nPoints])
    for (int index = 0; index < nPoints; index++)
    {
        U_pas[index] = U_pre[index];
        U_pre[index] = U_fut[index]; 
    }
}

void buildSeismogram(float * U_fut, int nxx, int nyy, int nzz, float * seismogram, int timeStep, int nt, int * rx, int * ry, int * rz, int nrec)
{
    # pragma acc parallel loop present(U_fut[0:nxx*nyy*nzz], seismogram[0:nt*nrec])
    for (int receiver = 0; receiver < nrec; receiver++)
    {
        int index = rz[receiver] + rx[receiver]*nzz + ry[receiver]*nxx*nzz;

        seismogram[timeStep + receiver*nt] = U_fut[index];
    } 
}

# endif