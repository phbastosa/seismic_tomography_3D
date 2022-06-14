# include <cmath>

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

# include "elastic.hpp"

void Elastic::forwardModeling()
{
    allocateVolumes();

    for (int shot = 0; shot < g3D.shots.n; shot++)
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
    M = new float[m3D.nPointsB]();
    L = new float[m3D.nPointsB]();
    Vx = new float[m3D.nPointsB]();
    Vy = new float[m3D.nPointsB]();
    Vz = new float[m3D.nPointsB]();
    Txx = new float[m3D.nPointsB]();
    Tyy = new float[m3D.nPointsB]();
    Tzz = new float[m3D.nPointsB]();
    Txy = new float[m3D.nPointsB]();
    Txz = new float[m3D.nPointsB]();
    Tyz = new float[m3D.nPointsB]();

    seismogram = new float[nt * g3D.nodes.n]();

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        M[index] = m3D.rho[index]*powf(m3D.vs[index],2.0f);
        L[index] = m3D.rho[index]*powf(m3D.vp[index],2.0f) - 2.0f*M[index];
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
    for (int i = 0; i < m3D.nPointsB; i++)
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
    for(int index = 0; index < m3D.nPointsB; index++)    
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
    for (int node = 0; node < g3D.nodes.n; node++)
    {
        int x = (int)(g3D.nodes.x[node] / m3D.dx) + m3D.nb;
        int y = (int)(g3D.nodes.y[node] / m3D.dy) + m3D.nb;
        int z = (int)(g3D.nodes.z[node] / m3D.dz) + m3D.nb;

        int index = z + x*m3D.nzz + y*m3D.nxx*m3D.nzz; 

        seismogram[timeId + node*nt] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
    }
}

void Elastic::elasticIsotropic_FD8E2T() 
{
    # pragma omp for 
    for(int index = 0; index < m3D.nPointsB; index++) 
    {    
        int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
        int j = (int) (index - k*m3D.nxx*m3D.nzz) / m3D.nzz;   // x direction
        int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

        if((index == 0) && (timeId < nsrc))
        {
            int x = (int)(g3D.shots.x[shotId] / m3D.dx) + m3D.nb;
            int y = (int)(g3D.shots.y[shotId] / m3D.dy) + m3D.nb;
            int z = (int)(g3D.shots.z[shotId] / m3D.dz) + m3D.nb;

            Txx[z + x*m3D.nzz + y*m3D.nxx*m3D.nzz] += source[timeId] / (m3D.dx * m3D.dy * m3D.dz);        
            Tyy[z + x*m3D.nzz + y*m3D.nxx*m3D.nzz] += source[timeId] / (m3D.dx * m3D.dy * m3D.dz);        
            Tzz[z + x*m3D.nzz + y*m3D.nxx*m3D.nzz] += source[timeId] / (m3D.dx * m3D.dy * m3D.dz);         
        }

        if((i >= 3) && (i < m3D.nzz-4) && (j >= 3) && (j < m3D.nxx-4) && (k >= 3) && (k < m3D.nyy-4)) 
        {    
            float d_Vx_dx = (75.0f*(Vx[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[i + (j+4)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vx[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vx[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vx[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float d_Vy_dy = (75.0f*(Vy[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz] - Vy[i + j*m3D.nzz + (k+4)*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vy[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz] - Vy[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vy[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz] - Vy[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vy[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] - Vy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float d_Vz_dz = (75.0f*(Vz[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[(i+4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vz[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vz[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vz[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            Txx[index] += dt*((L[index] + 2*M[index])*d_Vx_dx + L[index]*(d_Vy_dy + d_Vz_dz));
            Tyy[index] += dt*((L[index] + 2*M[index])*d_Vy_dy + L[index]*(d_Vx_dx + d_Vz_dz));
            Tzz[index] += dt*((L[index] + 2*M[index])*d_Vz_dz + L[index]*(d_Vx_dx + d_Vy_dy));                    
        }

        if((i >= 0) && (i <m3D.nzz) && (j > 3) && (j <m3D.nxx-3) && (k > 3) && (k < m3D.nyy-3)) 
        {    
            float d_Vx_dy = (75.0f*(Vx[i + j*m3D.nzz + (k-4)*m3D.nxx*m3D.nzz] - Vx[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vx[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz] - Vx[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vx[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz] - Vx[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vx[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vx[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float d_Vy_dx = (75.0f*(Vy[i + (j-4)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vy[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vy[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vy[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float M_xy = powf(0.25*(1/M[i + (j+1)*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] + 1/M[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] + 
                                    1/M[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz]     + 1/M[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]),-1.0f);

            Txy[index] += dt*M_xy*(d_Vx_dy + d_Vy_dx);
        }

        if((i > 3) && (i <m3D.nzz-3) &&  (j > 3) && (j <m3D.nxx-3) && (k >= 0) && (k < m3D.nyy)) 
        {
            float d_Vx_dz = (75.0f*(Vx[(i-4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vx[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vx[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vx[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vx[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vx[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            float d_Vz_dx = (75.0f*(Vz[i + (j-4)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vz[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vz[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vz[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vz[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float M_xz = powf(0.25*(1/M[(i+1) + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] + 1/M[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] + 
                                    1/M[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     + 1/M[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]),-1.0f);

            Txz[index] += dt*M_xz*(d_Vx_dz + d_Vz_dx);
        }
    
        if((i > 3) && (i <m3D.nzz-3) &&  (j >= 0) && (j <m3D.nxx) && (k > 3) && (k < m3D.nyy-3)) 
        {
            float d_Vy_dz = (75.0f*(Vy[(i-4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vy[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vy[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Vy[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vy[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            float d_Vz_dy = (75.0f*(Vz[i + j*m3D.nzz + (k-4)*m3D.nxx*m3D.nzz] - Vz[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz]) +
                           1029.0f*(Vz[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz] - Vz[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz]) +
                           8575.0f*(Vz[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz] - Vz[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz]) +
                         128625.0f*(Vz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Vz[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float M_yz = powf(0.25*(1/M[(i+1) + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] + 1/M[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] + 
                                    1/M[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] +     1/M[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]),-1.0f);

            Tyz[index] += dt*M_yz*(d_Vy_dz + d_Vz_dy);
        }
    }

    # pragma omp for 
    for(int index = 0; index < m3D.nPointsB; index++) 
    {    
        int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
        int j = (int) (index - k*m3D.nxx*m3D.nzz) /m3D.nzz;    // x direction
        int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

        if((i >= 3) && (i <m3D.nzz-4) && (j > 3) && (j <m3D.nxx-3) && (k >= 3) && (k < m3D.nyy-4)) 
        {
            float d_Txx_dx = (75.0f*(Txx[i + (j-4)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txx[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Txx[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txx[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Txx[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txx[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Txx[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Txx[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float d_Txy_dy = (75.0f*(Txy[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz] - Txy[i + j*m3D.nzz + (k+4)*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Txy[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz] - Txy[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Txy[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz] - Txy[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Txy[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] - Txy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float d_Txz_dz = (75.0f*(Txz[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[(i+4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Txz[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Txz[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Txz[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            float rhox = 0.5f*(m3D.rho[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] + m3D.rho[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]);

            Vx[index] += dt/rhox*(d_Txx_dx + d_Txy_dy + d_Txz_dz); 
        }
    
        if((i >= 3) && (i <m3D.nzz-3) && (j >= 3) && (j <m3D.nxx-4) && (k > 3) && (k < m3D.nyy-3)) 
        {
            float d_Txy_dx = (75.0f*(Txy[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txy[i + (j+4)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Txy[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txy[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Txy[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txy[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Txy[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float  d_Tyy_dy = (75.0f*(Tyy[i + j*m3D.nzz + (k-4)*m3D.nxx*m3D.nzz] - Tyy[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz]) +
                             1029.0f*(Tyy[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz] - Tyy[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz]) +
                             8575.0f*(Tyy[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz] - Tyy[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz]) +
                           128625.0f*(Tyy[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Tyy[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float d_Tyz_dz = (75.0f*(Tyz[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tyz[(i+4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Tyz[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tyz[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Tyz[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tyz[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Tyz[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tyz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            float rhoy = 0.5f*(m3D.rho[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] + m3D.rho[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]);

            Vy[index] += dt/rhoy*(d_Txy_dx + d_Tyy_dy + d_Tyz_dz); 
        }    

        if((i > 3) && (i <m3D.nzz-3) && (j >= 3) && (j <m3D.nxx-4) && (k >= 3) && (k < m3D.nyy-4)) 
        {
            float d_Txz_dx = (75.0f*(Txz[i + (j-3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[i + (j+4)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Txz[i + (j+3)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[i + (j-2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Txz[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[i + (j+2)*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Txz[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - Txz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dx);

            float d_Tyz_dy = (75.0f*(Tyz[i + j*m3D.nzz + (k-3)*m3D.nxx*m3D.nzz] - Tyz[i + j*m3D.nzz + (k+4)*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Tyz[i + j*m3D.nzz + (k+3)*m3D.nxx*m3D.nzz] - Tyz[i + j*m3D.nzz + (k-2)*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Tyz[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz] - Tyz[i + j*m3D.nzz + (k+2)*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Tyz[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] - Tyz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dy);

            float d_Tzz_dz = (75.0f*(Tzz[(i-4) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tzz[(i+3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            1029.0f*(Tzz[(i+2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tzz[(i-3) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                            8575.0f*(Tzz[(i-2) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - Tzz[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) +
                          128625.0f*(Tzz[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]     - Tzz[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]))/(107520.0f *m3D.dz);

            float rhoz = 0.5f*(m3D.rho[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] + m3D.rho[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz]);

            Vz[index] += dt/rhoz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
        }
    }    
}

