# include <cmath>
# include <iostream>

# include "elastic.hpp"

void Elastic::forwardModeling()
{
    allocateVolumes();

    for (int shot = 0; shot < shots.all; shot++)
    {
        shotId = shot;

        setWavefields();
 
        for (timeId = 0; timeId < nt; timeId++)
        {
            if (timeId % 100 == 0) std::cout << "propagation "<<timeId*dt<<std::endl;
            
            # pragma omp parallel 
            {
                elasticIsotropic_FD8E2T();

                cerjanAbsorbingCondition();

                getPressureSeismogram();
            }
        }
    }

    deleteVolumes();
}

void Elastic::allocateVolumes()
{
    M = new float[nPointsB]();
    L = new float[nPointsB]();
    Vx = new float[nPointsB]();
    Vy = new float[nPointsB]();
    Vz = new float[nPointsB]();
    Txx = new float[nPointsB]();
    Tyy = new float[nPointsB]();
    Tzz = new float[nPointsB]();
    Txy = new float[nPointsB]();
    Txz = new float[nPointsB]();
    Tyz = new float[nPointsB]();

    seismogram = new float[nt * nodes.all]();

    for (int index = 0; index < nPointsB; index++)
    {
        M[index] = rho[index]*powf(vs[index],2.0f);
        L[index] = rho[index]*powf(vp[index],2.0f) - 2.0f*M[index];
    }
}

void Elastic::deleteVolumes()
{
    delete[] Vx;
    delete[] Vy;
    delete[] Vz;
    delete[] Txx;
    delete[] Tyy;
    delete[] Tzz;
    delete[] Txy;
    delete[] Txz;
    delete[] Tyz;
    delete[] M;
    delete[] L;
}

void Elastic::setWavefields()
{
    for (int i = 0; i < nPointsB; i++)
    {
        Vx[i] = 0.0f;
        Vy[i] = 0.0f;
        Vz[i] = 0.0f;
        Txx[i] = 0.0f;
        Tyy[i] = 0.0f;
        Tzz[i] = 0.0f;
        Txy[i] = 0.0f;
        Txz[i] = 0.0f;
        Tyz[i] = 0.0f;
    }
}

void Elastic::cerjanAbsorbingCondition()
{
    # pragma omp for
    for(int index = 0; index < nPointsB; index++)    
    {
        Vx[index] *= cerjan[index];
        Vy[index] *= cerjan[index];
        Vz[index] *= cerjan[index];
        Txx[index] *= cerjan[index];
        Txy[index] *= cerjan[index];
        Txz[index] *= cerjan[index];
        Tyy[index] *= cerjan[index];
        Tyz[index] *= cerjan[index];
        Tzz[index] *= cerjan[index];
    }
}

void Elastic::getPressureSeismogram()
{
    # pragma omp for
    for (int node = 0; node < nodes.all; node++)
    {
        int x = (int)(nodes.x[node] / dx) + nb;
        int y = (int)(nodes.y[node] / dy) + nb;
        int z = (int)(nodes.z[node] / dz) + nb;

        int index = z + x*nzz + y*nxx*nzz; 

        seismogram[timeId + node*nt] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
    }
}

void Elastic::elasticIsotropic_FD8E2T() 
{
    # pragma omp for 
    for(int index = 0; index < nPointsB; index++) 
    {    
        int k = (int) (index / (nxx*nzz));             // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;       // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz);     // z direction

        if((index == 0) && (timeId < nsrc))
        {
            int x = (int)(shots.x[shotId] / dx) + nb;
            int y = (int)(shots.y[shotId] / dy) + nb;
            int z = (int)(shots.z[shotId] / dz) + nb;

            Txx[z + x*nzz + y*nxx*nzz] += source[timeId] / (dx * dy * dz);        
            Tyy[z + x*nzz + y*nxx*nzz] += source[timeId] / (dx * dy * dz);        
            Tzz[z + x*nzz + y*nxx*nzz] += source[timeId] / (dx * dy * dz);         
        }

        if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {    
            float d_Vx_dx = (75.0f*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                           1029.0f*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                           8575.0f*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                         128625.0f*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float d_Vy_dy = (75.0f*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                           1029.0f*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                           8575.0f*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                         128625.0f*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz]))/(107520.0f *dy);

            float d_Vz_dz = (75.0f*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                           1029.0f*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                           8575.0f*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                         128625.0f*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            Txx[index] += dt*((L[index] + 2*M[index])*d_Vx_dx + L[index]*(d_Vy_dy + d_Vz_dz));
            Tyy[index] += dt*((L[index] + 2*M[index])*d_Vy_dy + L[index]*(d_Vx_dx + d_Vz_dz));
            Tzz[index] += dt*((L[index] + 2*M[index])*d_Vz_dz + L[index]*(d_Vx_dx + d_Vy_dy));                    
        }

        if((i >= 0) && (i < nzz) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {    
            float d_Vx_dy = (75.0f*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                           1029.0f*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                           8575.0f*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                         128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz]))/(107520.0f *dy);

            float d_Vy_dx = (75.0f*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                           1029.0f*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                           8575.0f*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                         128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float M_xy = powf(0.25*(1/M[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1/M[i + j*nzz + (k+1)*nxx*nzz]     + 1/M[i + j*nzz + k*nxx*nzz]),-1.0f);

            Txy[index] += dt*M_xy*(d_Vx_dy + d_Vy_dx);
        }

        if((i > 3) && (i < nzz-3) &&  (j > 3) && (j < nxx-3) && (k >= 0) && (k < nyy)) 
        {
            float d_Vx_dz = (75.0f*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                           1029.0f*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                           8575.0f*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                         128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            float d_Vz_dx = (75.0f*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                           1029.0f*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                           8575.0f*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                         128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float M_xz = powf(0.25*(1/M[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1/M[(i+1) + j*nzz + k*nxx*nzz]     + 1/M[i + j*nzz + k*nxx*nzz]),-1.0f);

            Txz[index] += dt*M_xz*(d_Vx_dz + d_Vz_dx);
        }
    
        if((i > 3) && (i < nzz-3) &&  (j >= 0) && (j < nxx) && (k > 3) && (k < nyy-3)) 
        {
            float d_Vy_dz = (75.0f*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                           1029.0f*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                           8575.0f*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                         128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            float d_Vz_dy = (75.0f*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                           1029.0f*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                           8575.0f*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                         128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz]))/(107520.0f *dy);

            float M_yz = powf(0.25*(1/M[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1/M[i + j*nzz + (k+1)*nxx*nzz] + 
                                    1/M[(i+1) + j*nzz + k*nxx*nzz] +     1/M[i + j*nzz + k*nxx*nzz]),-1.0f);

            Tyz[index] += dt*M_yz*(d_Vy_dz + d_Vz_dy);
        }
    }

    # pragma omp for 
    for(int index = 0; index < nPointsB; index++) 
    {    
        int k = (int) (index / (nxx*nzz));             // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;       // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz);     // z direction

        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float d_Txx_dx = (75.0f*(Txx[i + (j-4)*nzz + k*nxx*nzz] - Txx[i + (j+3)*nzz + k*nxx*nzz]) +
                            1029.0f*(Txx[i + (j+2)*nzz + k*nxx*nzz] - Txx[i + (j-3)*nzz + k*nxx*nzz]) +
                            8575.0f*(Txx[i + (j-2)*nzz + k*nxx*nzz] - Txx[i + (j+1)*nzz + k*nxx*nzz]) +
                          128625.0f*(Txx[i + j*nzz + k*nxx*nzz]     - Txx[i + (j-1)*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float d_Txy_dy = (75.0f*(Txy[i + j*nzz + (k-3)*nxx*nzz] - Txy[i + j*nzz + (k+4)*nxx*nzz]) +
                            1029.0f*(Txy[i + j*nzz + (k+3)*nxx*nzz] - Txy[i + j*nzz + (k-2)*nxx*nzz]) +
                            8575.0f*(Txy[i + j*nzz + (k-1)*nxx*nzz] - Txy[i + j*nzz + (k+2)*nxx*nzz]) +
                          128625.0f*(Txy[i + j*nzz + (k+1)*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz]))/(107520.0f *dy);

            float d_Txz_dz = (75.0f*(Txz[(i-3) + j*nzz + k*nxx*nzz] - Txz[(i+4) + j*nzz + k*nxx*nzz]) +
                            1029.0f*(Txz[(i+3) + j*nzz + k*nxx*nzz] - Txz[(i-2) + j*nzz + k*nxx*nzz]) +
                            8575.0f*(Txz[(i-1) + j*nzz + k*nxx*nzz] - Txz[(i+2) + j*nzz + k*nxx*nzz]) +
                          128625.0f*(Txz[(i+1) + j*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            float rhox = 0.5f*(rho[i + (j+1)*nzz + k*nxx*nzz] + rho[i + j*nzz + k*nxx*nzz]);

            Vx[index] += dt/rhox*(d_Txx_dx + d_Txy_dy + d_Txz_dz); 
        }
    
        if((i >= 3) && (i <nzz-3) && (j >= 3) && (j <nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float d_Txy_dx = (75.0f*(Txy[i + (j-3)*nzz + k*nxx*nzz] - Txy[i + (j+4)*nzz + k*nxx*nzz]) +
                            1029.0f*(Txy[i + (j+3)*nzz + k*nxx*nzz] - Txy[i + (j-2)*nzz + k*nxx*nzz]) +
                            8575.0f*(Txy[i + (j-1)*nzz + k*nxx*nzz] - Txy[i + (j+2)*nzz + k*nxx*nzz]) +
                          128625.0f*(Txy[i + (j+1)*nzz + k*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float  d_Tyy_dy = (75.0f*(Tyy[i + j*nzz + (k-4)*nxx*nzz] - Tyy[i + j*nzz + (k+3)*nxx*nzz]) +
                             1029.0f*(Tyy[i + j*nzz + (k+2)*nxx*nzz] - Tyy[i + j*nzz + (k-3)*nxx*nzz]) +
                             8575.0f*(Tyy[i + j*nzz + (k-2)*nxx*nzz] - Tyy[i + j*nzz + (k+1)*nxx*nzz]) +
                           128625.0f*(Tyy[i + j*nzz + k*nxx*nzz]     - Tyy[i + j*nzz + (k-1)*nxx*nzz]))/(107520.0f *dy);

            float d_Tyz_dz = (75.0f*(Tyz[(i-3) + j*nzz + k*nxx*nzz] - Tyz[(i+4) + j*nzz + k*nxx*nzz]) +
                            1029.0f*(Tyz[(i+3) + j*nzz + k*nxx*nzz] - Tyz[(i-2) + j*nzz + k*nxx*nzz]) +
                            8575.0f*(Tyz[(i-1) + j*nzz + k*nxx*nzz] - Tyz[(i+2) + j*nzz + k*nxx*nzz]) +
                          128625.0f*(Tyz[(i+1) + j*nzz + k*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            float rhoy = 0.5f*(rho[i + j*nzz + (k+1)*nxx*nzz] + rho[i + j*nzz + k*nxx*nzz]);

            Vy[index] += dt/rhoy*(d_Txy_dx + d_Tyy_dy + d_Tyz_dz); 
        }    

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {
            float d_Txz_dx = (75.0f*(Txz[i + (j-3)*nzz + k*nxx*nzz] - Txz[i + (j+4)*nzz + k*nxx*nzz]) +
                            1029.0f*(Txz[i + (j+3)*nzz + k*nxx*nzz] - Txz[i + (j-2)*nzz + k*nxx*nzz]) +
                            8575.0f*(Txz[i + (j-1)*nzz + k*nxx*nzz] - Txz[i + (j+2)*nzz + k*nxx*nzz]) +
                          128625.0f*(Txz[i + (j+1)*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz]))/(107520.0f *dx);

            float d_Tyz_dy = (75.0f*(Tyz[i + j*nzz + (k-3)*nxx*nzz] - Tyz[i + j*nzz + (k+4)*nxx*nzz]) +
                            1029.0f*(Tyz[i + j*nzz + (k+3)*nxx*nzz] - Tyz[i + j*nzz + (k-2)*nxx*nzz]) +
                            8575.0f*(Tyz[i + j*nzz + (k-1)*nxx*nzz] - Tyz[i + j*nzz + (k+2)*nxx*nzz]) +
                          128625.0f*(Tyz[i + j*nzz + (k+1)*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz]))/(107520.0f *dy);

            float d_Tzz_dz = (75.0f*(Tzz[(i-4) + j*nzz + k*nxx*nzz] - Tzz[(i+3) + j*nzz + k*nxx*nzz]) +
                            1029.0f*(Tzz[(i+2) + j*nzz + k*nxx*nzz] - Tzz[(i-3) + j*nzz + k*nxx*nzz]) +
                            8575.0f*(Tzz[(i-2) + j*nzz + k*nxx*nzz] - Tzz[(i+1) + j*nzz + k*nxx*nzz]) +
                          128625.0f*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i-1) + j*nzz + k*nxx*nzz]))/(107520.0f *dz);

            float rhoz = 0.5f*(rho[(i+1) + j*nzz + k*nxx*nzz] + rho[i + j*nzz + k*nxx*nzz]);

            Vz[index] += dt/rhoz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
        }
    }    
}

