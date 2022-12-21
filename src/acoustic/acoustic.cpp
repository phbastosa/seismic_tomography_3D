# include <cmath>
# include <iostream>

# include "acoustic.hpp"

void Acoustic::sourceGenerator()
{
    nsrc = (int) (2.0f * tlag / dt + 1);

    source = new float[nsrc]();

    float pi = 4.0f * atanf(1.0f);
    float fc = fcut / (3.0f * sqrtf(pi));

    int s = (int) (nsrc / 2);

    for (int t = -s; t < s; t++)
    {
        float aux1 = 1.0f - 2.0f * pi * powf(t*dt, 2.0f) * powf(fc, 2.0f) * pow(pi, 2.0f);
        float aux2 = expf(-pi * powf(t*dt, 2.0f)*powf(fc, 2.0f)*powf(pi, 2.0f));

        source[t + s] = aux1 * aux2;
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

void Acoustic::setWaveField()
{
    for (int index = 0; index < nPointsB; index++)
    {
        U_pas[index] = 0.0f;
        U_pre[index] = 0.0f;
        U_fut[index] = 0.0f;
    }
}

void Acoustic::forwardModeling()
{
    U_pas = new float[nPointsB]();
    U_pre = new float[nPointsB]();
    U_fut = new float[nPointsB]();

    for (shotId = 0; shotId < shots.all; shotId++)
    {
        sIdx = (int) (shots.x[shotId]/dx + nb);    
        sIdy = (int) (shots.y[shotId]/dy + nb);    
        sIdz = (int) (shots.z[shotId]/dz + nb);    

        setWaveField();

        # pragma acc enter data copyin(this[0:1], V[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], source[0:nsrc])
        # pragma acc enter data copyin(this[0:1], U_pas[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], U_pre[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], U_fut[0:nPointsB])
        for (timeStep = 0; timeStep < nt; timeStep++)
        {
            progressMessage();
            
            wavePropagation();
            // dampApplication();
            wavefieldUpdate();
            // buildSeismogram();
        }
        # pragma acc exit data delete(V[0:nPointsB], this[0:1])
        # pragma acc exit data delete(source[0:nsrc], this[0:1])
        # pragma acc exit data delete(U_pas[0:nPointsB], this[0:1])
        # pragma acc exit data delete(U_pre[0:nPointsB], this[0:1])
        # pragma acc exit data copyout(U_fut[0:nPointsB], this[0:1])
    
        writeBinaryFloat("outputs/waveFiels3D.bin", U_fut, nPointsB);

        // exportSeismogram();
    }
}

void Acoustic::wavePropagation()
{
    # pragma acc parallel loop present(source[0:nsrc], V[0:nPointsB], U_pas[0:nPointsB], U_pre[0:nPointsB], U_fut[0:nPointsB])
    for (int index = 0; index < nPointsB; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction
        
        if ((index == 0) && (timeStep < nsrc)) 
            U_pre[sIdz + sIdx*nzz + sIdy*nxx*nzz] += source[timeStep] / (dx * dy * dz);

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

void Acoustic::dampApplication()
{


}


void Acoustic::wavefieldUpdate()
{
    # pragma acc parallel loop present(U_pas[0:nPointsB], U_pre[0:nPointsB], U_fut[0:nPointsB])
    for (int index = 0; index < nPointsB; index++)
    {
        U_pas[index] = U_pre[index];
        U_pre[index] = U_fut[index]; 
    }
}

void Acoustic::buildSeismogram()
{


}

void Acoustic::exportSeismogram()
{


}

void Acoustic::progressMessage()
{
    if (timeStep % (nt/5) == 0)
    {    
        system("clear");
        std::cout<<"Running shot "<<shotId+1<<" at position: (z = "<<shots.z[shotId]<<", x = "<<shots.x[shotId]<<", y = "<<shots.z[shotId]<<")"<<std::endl;
        std::cout<<"Time step: "<<timeStep*dt<<std::endl;
    }
}
