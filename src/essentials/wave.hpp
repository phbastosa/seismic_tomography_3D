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

void dampingGeneration(float * damp1D, float * damp2D, float * damp3D, float factor, int nb)
{
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

void wavePropagation(float * V, float * U_pas, float * U_pre, float * U_fut, int nxx, int nyy, int nzz, float dx, float dy, float dz, float dt)
{
    int nPoints = nxx*nyy*nzz;

    # pragma acc parallel loop present(V[0:nPoints],U_pas[0:nPoints],U_pre[0:nPoints],U_fut[0:nPoints])
    for (int index = 0; index < nPoints; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        if((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4)) 
        {
            float d2_Px2 = (- 9.0f*(U_pre[i + (j-4)*nzz + k*nxx*nzz] + U_pre[i + (j+4)*nzz + k*nxx*nzz])
                        +   128.0f*(U_pre[i + (j-3)*nzz + k*nxx*nzz] + U_pre[i + (j+3)*nzz + k*nxx*nzz])
                        -  1008.0f*(U_pre[i + (j-2)*nzz + k*nxx*nzz] + U_pre[i + (j+2)*nzz + k*nxx*nzz])
                        +  8064.0f*(U_pre[i + (j-1)*nzz + k*nxx*nzz] + U_pre[i + (j+1)*nzz + k*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

            float d2_Py2 = (- 9.0f*(U_pre[i + j*nzz + (k-4)*nxx*nzz] + U_pre[i + j*nzz + (k+4)*nxx*nzz])
                        +   128.0f*(U_pre[i + j*nzz + (k-3)*nxx*nzz] + U_pre[i + j*nzz + (k+3)*nxx*nzz])
                        -  1008.0f*(U_pre[i + j*nzz + (k-2)*nxx*nzz] + U_pre[i + j*nzz + (k+2)*nxx*nzz])
                        +  8064.0f*(U_pre[i + j*nzz + (k-1)*nxx*nzz] + U_pre[i + j*nzz + (k+1)*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

            float d2_Pz2 = (- 9.0f*(U_pre[(i-4) + j*nzz + k*nxx*nzz] + U_pre[(i+4) + j*nzz + k*nxx*nzz])
                        +   128.0f*(U_pre[(i-3) + j*nzz + k*nxx*nzz] + U_pre[(i+3) + j*nzz + k*nxx*nzz])
                        -  1008.0f*(U_pre[(i-2) + j*nzz + k*nxx*nzz] + U_pre[(i+2) + j*nzz + k*nxx*nzz])
                        +  8064.0f*(U_pre[(i-1) + j*nzz + k*nxx*nzz] + U_pre[(i+1) + j*nzz + k*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));

            U_fut[index] = powf(dt, 2.0f) * powf(V[index], 2.0f) * (d2_Px2 + d2_Py2 + d2_Pz2) + 2.0f*U_pre[index] - U_pas[index];
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