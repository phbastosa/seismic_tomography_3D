# include <cmath>
# include <algorithm>

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

# include "eikonal.hpp"

float Eikonal::min(float v1, float v2)
{
    if (v1 < v2)
        return v1;
    else
        return v2;
}

void Eikonal::allocateVolumes()
{
    T = new float[m3D.nPointsB]();    
    S = new float[m3D.nPointsB]();    
    K = new float[m3D.nPointsB]();    
    nT = new float[m3D.nPointsB]();    
    nK = new float[m3D.nPointsB]();    
}

void Eikonal::deleteVolumes()
{
    delete[] T;
    delete[] S;
    delete[] K;
    delete[] nT;
    delete[] nK;
}

void Eikonal::podvin()
{
    g3D.shots.idx = (int)(g3D.shots.x[shotId] / m3D.dx) + m3D.nb;
    g3D.shots.idy = (int)(g3D.shots.y[shotId] / m3D.dy) + m3D.nb;
    g3D.shots.idz = (int)(g3D.shots.z[shotId] / m3D.dz) + m3D.nb;

    int sId = g3D.shots.idz + g3D.shots.idx*m3D.nzz + g3D.shots.idy*m3D.nxx*m3D.nzz; 

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        S[index] = 1.0f / m3D.vp[index];

        if (index == sId)
        {
            float sx = floorf(g3D.shots.x[shotId] / m3D.dx) * m3D.dx;
            float sy = floorf(g3D.shots.y[shotId] / m3D.dy) * m3D.dy;
            float sz = floorf(g3D.shots.z[shotId] / m3D.dz) * m3D.dz;

            float dist = sqrtf(powf(sx - g3D.shots.x[shotId],2.0f) + powf(sy - g3D.shots.y[shotId],2.0f) + powf(sz - g3D.shots.z[shotId],2.0f));

            T[sId] = dist * S[sId];
            nT[sId] = T[sId]; 
        }
        else
        {
            T[index] = 1e6f;
            nT[index] = 1e6f;
        }
        
        K[index] = 0.0f;
        nK[index] = 0.0f;
    }

    K[sId - 1] = 1.0f;
    K[sId + 1] = 1.0f;
    K[sId - m3D.nzz] = 1.0f;
    K[sId + m3D.nzz] = 1.0f;
    K[sId - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 - m3D.nzz] = 1.0f;
    K[sId - 1 - m3D.nzz] = 1.0f;
    K[sId + 1 + m3D.nzz] = 1.0f;
    K[sId - 1 + m3D.nzz] = 1.0f;
    K[sId + 1 + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + 1 - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId - 1 + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;

    int aux = 0;
    int nItEikonal = 0;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    float sqrt2 = sqrtf(2.0f);
    float sqrt3 = sqrtf(3.0f);

    # pragma acc enter data copyin(this[0:1], S[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], K[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], nT[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], nK[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], T[0:m3D.nPointsB])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {  
            # pragma acc parallel loop present(S[0:m3D.nPointsB],T[0:m3D.nPointsB],K[0:m3D.nPointsB],nT[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
                    int j = (int) (index - k*m3D.nxx*m3D.nzz) / m3D.nzz;   // x direction
                    int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

                    if ((i > 0) && (i < m3D.nzz-1) && (j > 0) && (j < m3D.nxx-1) && (k > 0) && (k < m3D.nyy-1))
                    {
                        float h = m3D.dx;
                        float lowest = T[index];
                        float Tijk, T1, T2, Sref, M, N, P, Q, hs2; 

                        /* 1D operator head wave: i,j-1,k -> i,j,k (x direction) */
                        Tijk = T[index - m3D.nzz] + h*S[index - m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j+1,k -> i,j,k (x direction) */
                        Tijk = T[index + m3D.nzz] + h*S[index];
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k-1 -> i,j,k (y direction) */
                        Tijk = T[index - m3D.nxx*m3D.nzz] + h*S[index - m3D.nxx*m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k+1 -> i,j,k (y direction) */
                        Tijk = T[index + m3D.nxx*m3D.nzz] + h*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i-1,j,k -> i,j,k (z direction) */
                        Tijk = T[index - 1] + h*S[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i+1,j,k -> i,j,k (z direction) */
                        Tijk = T[index + 1] + h*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;
                    
                        /* 1D operator diffraction XZ plane */
                        
                        // i-1,j-1,k -> i,j,k
                        Tijk = T[index - 1 - m3D.nzz] + h*sqrt2*min(S[index - 1 - m3D.nzz], S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k -> i,j,k
                        Tijk = T[index - 1 + m3D.nzz] + h*sqrt2*min(S[index - 1], S[index - 1 - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j-1,k -> i,j,k
                        Tijk = T[index + 1 - m3D.nzz] + h*sqrt2*min(S[index - m3D.nzz], S[index - m3D.nzz - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j+1,k -> i,j,k
                        Tijk = T[index + 1 + m3D.nzz] + h*sqrt2*min(S[index], S[index - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction YZ plane */

                        // i-1,j,k-1 -> i,j,k
                        Tijk = T[index - 1 - m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index - 1], S[index - 1 - m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j,k+1 -> i,j,k
                        Tijk = T[index - 1 + m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index - 1 + m3D.nxx*m3D.nzz], S[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k-1 -> i,j,k
                        Tijk = T[index + 1 - m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index], S[index - m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k+1 -> i,j,k
                        Tijk = T[index + 1 + m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index + m3D.nxx*m3D.nzz], S[index - m3D.nzz + m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction XY plane */
                        
                        // i,j-1,k-1 -> i,j,k
                        Tijk = T[index - m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index - m3D.nzz - m3D.nxx*m3D.nzz], S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j-1,k+1 -> i,j,k
                        Tijk = T[index - m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index - m3D.nzz], S[index - 1 - m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k-1 -> i,j,k
                        Tijk = T[index + m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index - m3D.nxx*m3D.nzz], S[index - 1 - m3D.nxx*m3D.nzz]); 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k+1 -> i,j,k
                        Tijk = T[index + m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt2*min(S[index], S[index - 1]); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator corner diffractions */

                        // i-1,j-1,k-1 -> i,j,k
                        Tijk = T[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt3*S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j-1,k+1 -> i,j,k
                        Tijk = T[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt3*S[index - 1 - m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k-1 -> i,j,k
                        Tijk = T[index + 1 - m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt3*S[index - m3D.nzz - m3D.nxx*m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k+1 -> i,j,k
                        Tijk = T[index + 1 - m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt3*S[index - m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k-1 -> i,j,k
                        Tijk = T[index - 1 + m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt3*S[index - 1 - m3D.nxx*m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k+1 -> i,j,k
                        Tijk = T[index - 1 + m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt3*S[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k-1 -> i,j,k
                        Tijk = T[index + 1 + m3D.nzz - m3D.nxx*m3D.nzz] + h*sqrt3*S[index - m3D.nxx*m3D.nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k+1 -> i,j,k
                        Tijk = T[index + 1 + m3D.nzz + m3D.nxx*m3D.nzz] + h*sqrt3*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 2D operator XZ plane: First Quadrant*/

                        Sref = min(S[index - 1 - m3D.nzz], S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]);

                        // i,j-1,k - i-1,j-1,k -> i,j,k
                        T1 = T[index - m3D.nzz];
                        T2 = T[index - 1 - m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j-1,k -> i,j,k
                        T1 = T[index - 1];
                        T2 = T[index - 1 - m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Second Quadrant*/                        

                        Sref = min(S[index - m3D.nzz],S[index - m3D.nzz - m3D.nxx*m3D.nzz]);

                        // i,j-1,k - i+1,j-1,k -> i,j,k
                        T1 = T[index - m3D.nzz];
                        T2 = T[index + 1 - m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i+1,j,k - i+1,j-1,k -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 - m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Third Quadrant*/                        

                        Sref = min(S[index], S[index - m3D.nxx*m3D.nzz]);

                        // i+1,j,k - i+1,j+1,k -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 + m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i+1,j+1,k -> i,j,k
                        T1 = T[index + m3D.nzz];
                        T2 = T[index + 1 + m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Fourth Quadrant*/                        

                        Sref = min(S[index - 1], S[index - 1 - m3D.nxx*m3D.nzz]);

                        // i,j+1,k - i-1,j+1,k -> i,j,k
                        T1 = T[index + m3D.nzz];
                        T2 = T[index - 1 + m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j+1,k -> i,j,k
                        T1 = T[index - 1];
                        T2 = T[index - 1 + m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: First Quadrant */                        

                        Sref = min(S[index - 1 - m3D.nxx*m3D.nzz], S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]);

                        // i,j,k-1 - i-1,j,k-1 -> i,j,k
                        T1 = T[index - m3D.nxx*m3D.nzz];
                        T2 = T[index - 1 - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j,k-1 -> i,j,k
                        T1 = T[index - 1];
                        T2 = T[index - 1 - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Second Quadrant */                        

                        Sref = min(S[index - m3D.nxx*m3D.nzz], S[index - m3D.nzz - m3D.nxx*m3D.nzz]);

                        // i,j,k-1 - i+1,j,k-1 -> i,j,k
                        T1 = T[index - m3D.nxx*m3D.nzz];
                        T2 = T[index + 1 - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i+1,j,k - i+1,j,k-1 -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Third Quadrant*/                        

                        Sref = min(S[index], S[index - m3D.nzz]);

                        // i+1,j,k - i+1,j,k+1 -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i+1,j,k+1 -> i,j,k
                        T1 = T[index + m3D.nxx*m3D.nzz];
                        T2 = T[index + 1 + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Fourth Quadrant*/                        

                        Sref = min(S[index - 1], S[index - 1 - m3D.nzz]);

                        // i,j,k+1 - i-1,j,k+1 -> i,j,k
                        T1 = T[index + m3D.nxx*m3D.nzz];
                        T2 = T[index - 1 + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j,k+1 -> i,j,k
                        T1 = T[index - 1];
                        T2 = T[index - 1 + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: First Quadrant*/                        

                        Sref = min(S[index - m3D.nzz - m3D.nxx*m3D.nzz],S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz]);

                        // i,j-1,k - i,j-1,k-1 -> i,j,k
                        T1 = T[index - m3D.nzz];
                        T2 = T[index - m3D.nzz - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j-1,k-1 -> i,j,k
                        T1 = T[index - m3D.nxx*m3D.nzz];
                        T2 = T[index - m3D.nzz - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Second Quadrant*/                        

                        Sref = min(S[index - m3D.nzz],S[index - 1 - m3D.nzz]);

                        // i,j-1,k - i,j-1,k+1 -> i,j,k
                        T1 = T[index - m3D.nzz];
                        T2 = T[index - m3D.nzz + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i,j-1,k+1 -> i,j,k
                        T1 = T[index + m3D.nxx*m3D.nzz];
                        T2 = T[index - m3D.nzz + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Third Quadrant*/                        

                        Sref = min(S[index],S[index - 1]);

                        // i,j,k+1 - i,j+1,k+1 -> i,j,k
                        T1 = T[index + m3D.nxx*m3D.nzz];
                        T2 = T[index + m3D.nzz + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i,j+1,k+1 -> i,j,k
                        T1 = T[index + m3D.nzz];
                        T2 = T[index + m3D.nzz + m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Fourth Quadrant*/                        

                        Sref = min(S[index - m3D.nxx*m3D.nzz], S[index - 1 - m3D.nxx*m3D.nzz]);

                        // i,j+1,k - i,j+1,k-1 -> i,j,k
                        T1 = T[index + m3D.nzz];
                        T2 = T[index + m3D.nzz - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j+1,k-1 -> i,j,k
                        T1 = T[index - m3D.nxx*m3D.nzz];
                        T2 = T[index + m3D.nzz - m3D.nxx*m3D.nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 3D operator - First octant: XY plane */

                        Sref = S[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k-1 */   M = T[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - m3D.nzz];             
    /* i-1, j ,k-1 */   P = T[index - 1 - m3D.nxx*m3D.nzz];       
    /* i-1, j , k  */   Q = T[index - 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - First octant: YZ plane */

    /* i-1,j-1,k-1 */   M = T[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - m3D.nzz];             
    /*  i ,j-1,k-1 */   P = T[index - m3D.nzz - m3D.nxx*m3D.nzz];       
    /*  i ,j-1, k  */   Q = T[index - m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - First octant: XZ plane */

    /* i-1,j-1,k-1 */   M = T[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /*  i ,j-1,k-1 */   N = T[index - m3D.nzz - m3D.nxx*m3D.nzz];             
    /* i-1, j ,k-1 */   P = T[index - 1 - m3D.nxx*m3D.nzz];       
    /*  i , j ,k-1 */   Q = T[index - m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Second octant: XY plane */

                        Sref = S[index - 1 - m3D.nxx*m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j+1,k-1 */   M = T[index - 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i-1, j ,k-1 */   N = T[index - 1 - m3D.nxx*m3D.nzz];             
    /* i-1,j+1, k  */   P = T[index - 1 + m3D.nzz];       
    /* i-1, j , k  */   Q = T[index - 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Second octant: YZ plane */

    /* i-1,j+1,k-1 */   M = T[index - 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i-1,j+1, k  */   N = T[index - 1 + m3D.nzz];             
    /*  i ,j+1,k-1 */   P = T[index + m3D.nzz - m3D.nxx*m3D.nzz];       
    /*  i ,j+1, k  */   Q = T[index + m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Second octant: XZ plane */

    /* i-1,j+1,k-1 */   M = T[index - 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i-1, j ,k-1 */   N = T[index - 1 - m3D.nxx*m3D.nzz];             
    /*  i ,j+1,k-1 */   P = T[index + m3D.nzz - m3D.nxx*m3D.nzz];       
    /*  i , j ,k-1 */   Q = T[index - m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Third octant: XY plane */

                        Sref = S[index - 1];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j+1,k+1 */   M = T[index - 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i-1,j+1, k  */   N = T[index - 1 + m3D.nzz];             
    /* i-1, j ,k+1 */   P = T[index - 1 + m3D.nxx*m3D.nzz];       
    /* i-1, j , k  */   Q = T[index - 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Third octant: YZ plane */

    /* i-1,j+1,k+1 */   M = T[index - 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /*  i ,j+1,k+1 */   N = T[index + m3D.nzz + m3D.nxx*m3D.nzz];             
    /* i-1,j+1, k  */   P = T[index - 1 + m3D.nzz];       
    /*  i ,j+1, k  */   Q = T[index + m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Third octant: XZ plane */

    /* i-1,j+1,k+1 */   M = T[index - 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i-1, j ,k+1 */   N = T[index - 1 + m3D.nxx*m3D.nzz];             
    /*  i ,j+1,k+1 */   P = T[index + m3D.nzz + m3D.nxx*m3D.nzz];       
    /*  i , j ,k+1 */   Q = T[index + m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fourth octant: XY plane */

                        Sref = S[index - 1 - m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k+1 */   M = T[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i-1, j ,k+1 */   N = T[index - 1 + m3D.nxx*m3D.nzz];             
    /* i-1,j-1, k  */   P = T[index - 1 - m3D.nzz];       
    /* i-1, j , k  */   Q = T[index - 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fourth octant: YZ plane */

    /* i-1,j-1,k+1 */   M = T[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - m3D.nzz];             
    /*  i ,j-1,k+1 */   P = T[index - m3D.nzz + m3D.nxx*m3D.nzz];       
    /*  i ,j-1, k  */   Q = T[index - m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fourth octant: XZ plane */

    /* i-1,j-1,k+1 */   M = T[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /*  i ,j-1,k+1 */   N = T[index - m3D.nzz + m3D.nxx*m3D.nzz];             
    /* i-1, j ,k+1 */   P = T[index - 1 + m3D.nxx*m3D.nzz];       
    /*  i , j ,k+1 */   Q = T[index + m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fifth octant: XY plane */

                        Sref = S[index - m3D.nzz - m3D.nxx*m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k-1 */   M = T[index + 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i+1, j ,k-1 */   N = T[index + 1 - m3D.nxx*m3D.nzz];             
    /* i+1,j-1, k  */   P = T[index + 1 - m3D.nzz];       
    /* i+1, j , k  */   Q = T[index + 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fifth octant: YZ plane */

    /* i+1,j-1,k-1 */   M = T[index + 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i+1,j-1, k  */   N = T[index + 1 - m3D.nzz];             
    /*  i ,j-1,k-1 */   P = T[index - m3D.nzz - m3D.nxx*m3D.nzz];       
    /*  i ,j-1, k  */   Q = T[index - m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fifth octant: XZ plane */

    /* i+1,j-1,k-1 */   M = T[index + 1 - m3D.nzz - m3D.nxx*m3D.nzz];   
    /*  i ,j-1,k-1 */   N = T[index - m3D.nzz - m3D.nxx*m3D.nzz];             
    /* i+1, j ,k-1 */   P = T[index + 1 - m3D.nxx*m3D.nzz];       
    /*  i , j ,k-1 */   Q = T[index - m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Sixth octant: XY plane */

                        Sref = S[index - m3D.nxx*m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j+1,k-1 */   M = T[index + 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i+1,j+1, k  */   N = T[index + 1 + m3D.nzz];             
    /* i+1, j ,k-1 */   P = T[index + 1 - m3D.nxx*m3D.nzz];       
    /* i+1, j , k  */   Q = T[index + 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Sixth octant: YZ plane */

    /* i+1,j+1,k-1 */   M = T[index + 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /*  i ,j+1,k-1 */   N = T[index + m3D.nzz - m3D.nxx*m3D.nzz];             
    /* i+1,j+1, k  */   P = T[index + 1 + m3D.nzz];       
    /*  i ,j+1, k  */   Q = T[index + m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Sixth octant: XZ plane */

    /* i+1,j+1,k-1 */   M = T[index + 1 + m3D.nzz - m3D.nxx*m3D.nzz];   
    /* i+1, j ,k-1 */   N = T[index + 1 - m3D.nxx*m3D.nzz];             
    /*  i ,j+1,k-1 */   P = T[index + m3D.nzz - m3D.nxx*m3D.nzz];       
    /*  i , j ,k-1 */   Q = T[index - m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Seventh octant: XY plane */
                        
                        Sref = S[index - m3D.nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k+1 */   M = T[index + 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i+1,j-1, k  */   N = T[index + 1 - m3D.nzz];             
    /* i+1, j ,k+1 */   P = T[index + 1 + m3D.nxx*m3D.nzz];       
    /* i+1, j , k  */   Q = T[index + 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Seventh octant: YZ plane */

    /* i+1,j-1,k+1 */   M = T[index + 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /*  i ,j-1,k+1 */   N = T[index - m3D.nzz + m3D.nxx*m3D.nzz];             
    /* i+1,j-1, k  */   P = T[index + 1 - m3D.nzz];       
    /*  i ,j-1, k  */   Q = T[index - m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Seventh octant: XZ plane */

    /* i+1,j-1,k+1 */   M = T[index + 1 - m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i+1, j ,k+1 */   N = T[index + 1 + m3D.nxx*m3D.nzz];             
    /*  i ,j-1,k+1 */   P = T[index - m3D.nzz + m3D.nxx*m3D.nzz];       
    /*  i , j ,k+1 */   Q = T[index + m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Eighth octant: XY plane */

                        Sref = S[index];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j+1,k+1 */   M = T[index + 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i+1, j ,k+1 */   N = T[index + 1 + m3D.nxx*m3D.nzz];             
    /* i+1,j+1, k  */   P = T[index + 1 + m3D.nzz];       
    /* i+1, j , k  */   Q = T[index + 1];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Eighth octant: YZ plane */

    /* i+1,j+1,k+1 */   M = T[index + 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /* i+1,j+1, k  */   N = T[index + 1 + m3D.nzz];             
    /*  i ,j+1,k+1 */   P = T[index + m3D.nzz + m3D.nxx*m3D.nzz];       
    /*  i ,j+1, k  */   Q = T[index + m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Eighth octant: XZ plane */

    /* i+1,j+1,k+1 */   M = T[index + 1 + m3D.nzz + m3D.nxx*m3D.nzz];   
    /*  i ,j+1,k+1 */   N = T[index + m3D.nzz + m3D.nxx*m3D.nzz];             
    /* i+1, j ,k+1 */   P = T[index + 1 + m3D.nxx*m3D.nzz];       
    /*  i , j ,k+1 */   Q = T[index + m3D.nxx*m3D.nzz];                 

                        // MNP -> R 
                        if ((M <= N) && (M <= P) && 
                           ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                           ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                           ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
                        {
                            Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* Time atualization */
                        if (lowest == T[index]) K[index] = 0.0f;

                        nT[index] = lowest;
                    }
                }
            }

            # pragma acc parallel loop present(nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++) nK[index] = 0.0f;

            # pragma acc parallel loop present(K[0:m3D.nPointsB], nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
                    int j = (int) (index - k*m3D.nxx*m3D.nzz) / m3D.nzz;   // x direction
                    int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

                    if ((i > 0) && (i < m3D.nzz-1) && (j > 0) && (j < m3D.nxx-1) && (k > 0) && (k < m3D.nyy-1))
                    {
                        nK[index - 1] = 1.0f;
                        nK[index + 1] = 1.0f;
                        nK[index - m3D.nzz] = 1.0f;
                        nK[index + m3D.nzz] = 1.0f;
                        nK[index - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 - m3D.nzz] = 1.0f;
                        nK[index - 1 - m3D.nzz] = 1.0f;
                        nK[index + 1 + m3D.nzz] = 1.0f;
                        nK[index - 1 + m3D.nzz] = 1.0f;
                        nK[index + 1 + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + 1 - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 - m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 - m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 + m3D.nzz - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index - 1 + m3D.nzz + m3D.nxx*m3D.nzz] = 1.0f;
                    }
                }
            }

            # pragma acc parallel loop present(T[0:m3D.nPointsB],nT[0:m3D.nPointsB],K[0:m3D.nPointsB],nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(S[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(K[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(nT[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(nK[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data copyout(T[0:m3D.nPointsB], this[0:1])

    writeTravelTimes();
    writeFirstArrivals();
}

void Eikonal::jeongFIM()
{
    g3D.shots.idx = (int)(g3D.shots.x[shotId] / m3D.dx) + m3D.nb;
    g3D.shots.idy = (int)(g3D.shots.y[shotId] / m3D.dy) + m3D.nb;
    g3D.shots.idz = (int)(g3D.shots.z[shotId] / m3D.dz) + m3D.nb;

    int sId = g3D.shots.idz + g3D.shots.idx*m3D.nzz + g3D.shots.idy*m3D.nxx*m3D.nzz; 

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        S[index] = 1.0f / m3D.vp[index];

        if (index == sId)
        {
            float sx = floorf(g3D.shots.x[shotId] / m3D.dx) * m3D.dx;
            float sy = floorf(g3D.shots.y[shotId] / m3D.dy) * m3D.dy;
            float sz = floorf(g3D.shots.z[shotId] / m3D.dz) * m3D.dz;

            float dist = sqrtf(powf(sx - g3D.shots.x[shotId],2.0f) + powf(sy - g3D.shots.y[shotId],2.0f) + powf(sz - g3D.shots.z[shotId],2.0f));

            T[sId] = dist * S[sId];
            nT[sId] = T[sId]; 
        }
        else
        {
            T[index] = 1e6f;
            nT[index] = 1e6f;
        }
        
        K[index] = 0.0f;
        nK[index] = 0.0f;
    }

    K[sId - 1] = 1.0f;
    K[sId + 1] = 1.0f;
    K[sId - m3D.nzz] = 1.0f;
    K[sId + m3D.nzz] = 1.0f;
    K[sId - m3D.nxx*m3D.nzz] = 1.0f;
    K[sId + m3D.nxx*m3D.nzz] = 1.0f;      
 
    int aux = 0;
    int nItEikonal = 0;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(m3D.nxx - g3D.shots.idx,2.0f) + powf(m3D.nyy - g3D.shots.idy,2.0f) + powf(m3D.nzz - g3D.shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    # pragma acc enter data copyin(this[0:1], S[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], K[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], nT[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], nK[0:m3D.nPointsB])
    # pragma acc enter data copyin(this[0:1], T[0:m3D.nPointsB])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {
            # pragma acc parallel loop present(S[0:m3D.nPointsB],T[0:m3D.nPointsB],K[0:m3D.nPointsB],nT[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
                    int j = (int) (index - k*m3D.nxx*m3D.nzz) / m3D.nzz;   // x direction
                    int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

                    if ((i > 0) && (i < m3D.nzz-1) && (j > 0) && (j < m3D.nxx-1) && (k > 0) && (k < m3D.nyy-1))
                    {
                        float h = m3D.dx;
                        float a, b, c, tmp, Tijk;
                        float tlag = - 0.4f * h * S[sId];

                        a = min(T[index - m3D.nzz],T[index + m3D.nzz]);                 // Tx min        
                        b = min(T[index - m3D.nxx*m3D.nzz],T[index + m3D.nxx*m3D.nzz]); // Ty min        
                        c = min(T[index - 1],T[index + 1]);                             // Tz min        

                        // a,b,c <------- sort(Tx,Ty,Tz), where a > b > c
                        if (a < b) {tmp = a; a = b; b = tmp;}
                        if (b < c) {tmp = b; b = c; c = tmp;}
                        if (a < b) {tmp = a; a = b; b = tmp;}

                        Tijk = 1e6;

                        if (c < 1e6)
                        {
                            Tijk = c + h*S[index];

                            if (Tijk > b)
                            {
                                tmp = 0.5f * (b + c + sqrtf(2.0f*h*h*S[index]*S[index] - (b - c)*(b - c)));           

                                if (tmp > b) Tijk = tmp;

                                if (Tijk > a)
                                {
                                    tmp = (a + b + c)/3.0f + sqrtf(2.0f*(a*(b - a) + b*(c - b) + c*(a - c)) + 3.0f*h*h*S[index]*S[index])/3.0f;

                                    if (tmp > a) Tijk = tmp;
                                }
                            }
                        }

                        /* Time atualization */
                        float lowest = min(Tijk,T[index]);    

                        if (lowest == T[index]) K[index] = 0.0f;

                        nT[index] = lowest;
                    }
                }
            }

            # pragma acc parallel loop present(nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++) nK[index] = 0.0f;

            # pragma acc parallel loop present(K[0:m3D.nPointsB], nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (m3D.nxx*m3D.nzz));             // y direction
                    int j = (int) (index - k*m3D.nxx*m3D.nzz) / m3D.nzz;   // x direction
                    int i = (int) (index - j*m3D.nzz - k*m3D.nxx*m3D.nzz); // z direction

                    if ((i > 0) && (i < m3D.nzz-1) && (j > 0) && (j < m3D.nxx-1) && (k > 0) && (k < m3D.nyy-1))
                    {
                        nK[index - 1] = 1.0f;
                        nK[index + 1] = 1.0f;
                        nK[index - m3D.nzz] = 1.0f;
                        nK[index + m3D.nzz] = 1.0f;
                        nK[index - m3D.nxx*m3D.nzz] = 1.0f;
                        nK[index + m3D.nxx*m3D.nzz] = 1.0f;      
                    }
                }
            }

            # pragma acc parallel loop present(T[0:m3D.nPointsB],nT[0:m3D.nPointsB],K[0:m3D.nPointsB],nK[0:m3D.nPointsB])
            for (int index = 0; index < m3D.nPointsB; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(S[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(K[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(nT[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data delete(nK[0:m3D.nPointsB], this[0:1])
    # pragma acc exit data copyout(T[0:m3D.nPointsB], this[0:1])

    writeTravelTimes();
    writeFirstArrivals();
} 

void Eikonal::innerSweep()
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = fsm.i - fsm.sgnvz; 
    int j1 = fsm.j - fsm.sgnvx; 
    int k1 = fsm.k - fsm.sgnvy;

    // Get local times of surrounding points
    float tv = T[(fsm.i - fsm.sgntz) + fsm.j*m3D.nzz + fsm.k*m3D.nxx*m3D.nzz];
    float te = T[fsm.i + (fsm.j - fsm.sgntx)*m3D.nzz + fsm.k*m3D.nxx*m3D.nzz];
    float tn = T[fsm.i + fsm.j*m3D.nzz + (fsm.k - fsm.sgnty)*m3D.nxx*m3D.nzz];
    float tev = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*m3D.nzz + fsm.k*m3D.nxx*m3D.nzz];
    float ten = T[fsm.i + (fsm.j - fsm.sgntx)*m3D.nzz + (fsm.k - fsm.sgnty)*m3D.nxx*m3D.nzz];
    float tnv = T[(fsm.i - fsm.sgntz) + fsm.j*m3D.nzz + (fsm.k - fsm.sgnty)*m3D.nxx*m3D.nzz];
    float tnve = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*m3D.nzz + (fsm.k - fsm.sgnty)*m3D.nxx*m3D.nzz];     

    float Tijk = T[fsm.i + fsm.j*m3D.nzz + fsm.k*m3D.nxx*m3D.nzz];

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + m3D.dz * utils.min4(S[i1 + utils.imax(fsm.j-1,1)*m3D.nzz       + utils.imax(fsm.k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[i1 + utils.imax(fsm.j-1,1)*m3D.nzz       + utils.imin(fsm.k,m3D.nyy-1)*m3D.nxx*m3D.nzz],
                                    S[i1 + utils.imin(fsm.j,m3D.nxx-1)*m3D.nzz + utils.imax(fsm.k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[i1 + utils.imin(fsm.j,m3D.nxx-1)*m3D.nzz + utils.imin(fsm.k,m3D.nyy-1)*m3D.nxx*m3D.nzz]);

    // X direction
    t1D2 = te + m3D.dx * utils.min4(S[utils.imax(fsm.i-1,1)       + j1*m3D.nzz + utils.imax(fsm.k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(fsm.i,m3D.nzz-1) + j1*m3D.nzz + utils.imax(fsm.k-1,1)*m3D.nxx*m3D.nzz],
                                    S[utils.imax(fsm.i-1,1)       + j1*m3D.nzz + utils.imin(fsm.k,m3D.nyy-1)*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(fsm.i,m3D.nzz-1) + j1*m3D.nzz + utils.imin(fsm.k,m3D.nyy-1)*m3D.nxx*m3D.nzz]);

    // Y direction
    t1D3 = tn + m3D.dy * utils.min4(S[utils.imax(fsm.i-1,1)       + utils.imax(fsm.j-1,1)*m3D.nzz       + k1*m3D.nxx*m3D.nzz], 
                                    S[utils.imax(fsm.i-1,1)       + utils.imin(fsm.j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz],
                                    S[utils.imin(fsm.i,m3D.nzz-1) + utils.imax(fsm.j-1,1)*m3D.nzz       + k1*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(fsm.i,m3D.nzz-1) + utils.imin(fsm.j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    t1D = utils.min3(t1D1,t1D2,t1D3);

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = utils.min(S[i1 + j1*m3D.nzz + utils.imax(fsm.k-1,1)*m3D.nxx*m3D.nzz], S[i1 + j1*m3D.nzz + utils.imin(fsm.k, m3D.nyy-1)*m3D.nxx*m3D.nzz]);
    
    if ((tv < te + m3D.dx*Sref) && (te < tv + m3D.dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*fsm.dz2i + ta*fsm.dx2i) + sqrt(4.0f*Sref*Sref*(fsm.dz2i + fsm.dx2i) - fsm.dz2i*fsm.dx2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = utils.min(S[i1 + utils.imax(fsm.j-1,1)*m3D.nzz + k1*m3D.nxx*m3D.nzz], S[i1 + utils.imin(fsm.j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    if((tv < tn + m3D.dy*Sref) && (tn < tv + m3D.dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*fsm.dz2i + tb*fsm.dy2i) + sqrt(4.0f*Sref*Sref*(fsm.dz2i + fsm.dy2i) - fsm.dz2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = utils.min(S[utils.imax(fsm.i-1,1) + j1*m3D.nzz + k1*m3D.nxx*m3D.nzz],S[utils.imin(fsm.i,m3D.nzz-1) + j1*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    if((te < tn + m3D.dy*Sref) && (tn < te + m3D.dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*fsm.dx2i + tb*fsm.dy2i) + sqrt(4.0f*Sref*Sref*(fsm.dx2i + fsm.dy2i) - fsm.dx2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dx2i + fsm.dy2i);
    }

    t2D = utils.min3(t2D1,t2D2,t2D3);

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[i1 + j1*m3D.nxx + k1*m3D.nxx*m3D.nzz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (utils.min(t1D,t2D) > utils.max3(tv,te,tn))
    {
        t2 = 9.0f*Sref*Sref*fsm.dsum;
        
        t3 = fsm.dz2dx2*(ta - tb)*(ta - tb) + fsm.dz2dy2*(tb - tc)*(tb - tc) + fsm.dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*fsm.dz2i + ta*fsm.dx2i + tc*fsm.dy2i;        
            
            t3D = (t1 + sqrtf(t2 - t3)) / fsm.dsum;
        }
    }
   
    T[fsm.i + fsm.j*m3D.nzz + fsm.k*m3D.nxx*m3D.nzz] = utils.min4(Tijk,t1D,t2D,t3D);
}

void Eikonal::initSweep()
{
    // First sweeping: Top->Bottom; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = 1; 
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = utils.imax(1, g3D.shots.idy); fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = utils.imax(1, g3D.shots.idx); fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = utils.imax(1, g3D.shots.idz); fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = utils.imax(1, g3D.shots.idy); fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = utils.imax(1, g3D.shots.idx); fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = g3D.shots.idz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = g3D.shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = utils.imax(1, g3D.shots.idx); fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = utils.imax(1, g3D.shots.idz); fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = g3D.shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = utils.imax(1, g3D.shots.idx); fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = g3D.shots.idz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = utils.imax(1, g3D.shots.idy); fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = g3D.shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = utils.imax(1, g3D.shots.idz); fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = utils.imax(1, g3D.shots.idy); fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = g3D.shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = g3D.shots.idz + 1; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = g3D.shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = g3D.shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = utils.imax(1, g3D.shots.idz); fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = g3D.shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = g3D.shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = g3D.shots.idz + 1; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }
}

void Eikonal::fullSweep()
{
    // First sweeping: Top->Bottom; West->East; South->North 
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = 1; 
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = m3D.nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = m3D.nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = m3D.nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < m3D.nxx; fsm.j++)
        {
            for (fsm.i = m3D.nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = m3D.nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < m3D.nyy; fsm.k++)
    {
        for (fsm.j = m3D.nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = m3D.nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = m3D.nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = m3D.nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < m3D.nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = m3D.nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = m3D.nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = m3D.nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }
}

void Eikonal::nobleFSM()
{
    g3D.shots.idx = (int)(g3D.shots.x[shotId] / m3D.dx);
    g3D.shots.idy = (int)(g3D.shots.y[shotId] / m3D.dy);
    g3D.shots.idz = (int)(g3D.shots.z[shotId] / m3D.dz);

    int sId = (g3D.shots.idz + m3D.nb) + (g3D.shots.idx + m3D.nb)*m3D.nzz + (g3D.shots.idy + m3D.nb)*m3D.nxx*m3D.nzz;     

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        T[index] = 1e6f;
        S[index] = 1.0f / m3D.vp[index];
    }

    T[sId] = S[sId] * sqrtf(powf(g3D.shots.idx*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf(g3D.shots.idy*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf(g3D.shots.idz*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + 1] = S[sId] * sqrtf(powf(g3D.shots.idx*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf(g3D.shots.idy*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf((g3D.shots.idz+1)*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + m3D.nzz] = S[sId] * sqrtf(powf((g3D.shots.idx+1)*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf(g3D.shots.idy*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf(g3D.shots.idz*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + m3D.nxx*m3D.nzz] = S[sId] * sqrtf(powf(g3D.shots.idx*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf((g3D.shots.idy+1)*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf(g3D.shots.idz*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + 1 + m3D.nzz] = S[sId] * sqrtf(powf((g3D.shots.idx+1)*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf(g3D.shots.idy*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf((g3D.shots.idz+1)*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + 1 + m3D.nxx*m3D.nzz] = S[sId] * sqrtf(powf(g3D.shots.idx*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf((g3D.shots.idy+1)*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf((g3D.shots.idz+1)*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + m3D.nzz + m3D.nxx*m3D.nzz] = S[sId] * sqrtf(powf((g3D.shots.idx+1)*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf((g3D.shots.idy+1)*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf(g3D.shots.idz*m3D.dz - g3D.shots.z[shotId], 2.0f));
    T[sId + 1 + m3D.nzz + m3D.nxx*m3D.nzz] = S[sId] * sqrtf(powf((g3D.shots.idx+1)*m3D.dx - g3D.shots.x[shotId], 2.0f) + powf((g3D.shots.idy+1)*m3D.dy - g3D.shots.y[shotId], 2.0f) + powf((g3D.shots.idz+1)*m3D.dz - g3D.shots.z[shotId], 2.0f));

    fsm.dzi = 1.0f / m3D.dz;
    fsm.dxi = 1.0f / m3D.dx;
    fsm.dyi = 1.0f / m3D.dy;
    fsm.dz2i = 1.0f / (m3D.dz*m3D.dz);
    fsm.dx2i = 1.0f / (m3D.dx*m3D.dx);
    fsm.dy2i = 1.0f / (m3D.dy*m3D.dy);
    fsm.dz2dx2 = fsm.dz2i * fsm.dx2i;
    fsm.dz2dy2 = fsm.dz2i * fsm.dy2i;
    fsm.dx2dy2 = fsm.dx2i * fsm.dy2i;
    fsm.dsum = fsm.dz2i + fsm.dx2i + fsm.dy2i;

    g3D.shots.idx += m3D.nb;
    g3D.shots.idy += m3D.nb;
    g3D.shots.idz += m3D.nb;

    initSweep();
    
    for (int s = 0; s < 5; s++) 
        fullSweep();

    writeTravelTimes();
    writeFirstArrivals();
}

void Eikonal::writeTravelTimes()
{
    if (exportTimesVolume)
    {    
        float * TT = new float[m3D.nPointsB];

        for (int indb = 0; indb < m3D.nPointsB; indb++)
        {
            int yb = (int) (indb / (m3D.nxx*m3D.nzz));              // y direction
            int xb = (int) (indb - yb*m3D.nxx*m3D.nzz) / m3D.nzz;    // x direction
            int zb = (int) (indb - xb*m3D.nzz - yb*m3D.nxx*m3D.nzz); // z direction

            if ((zb >= m3D.nb) && (zb < m3D.nzz - m3D.nb) && (yb >= m3D.nb) && (yb < m3D.nyy - m3D.nb) && (xb >= m3D.nb) && (xb < m3D.nxx - m3D.nb))
            {
                TT[(zb - m3D.nb) + (xb - m3D.nb)*m3D.nz + (yb - m3D.nb)*m3D.nx*m3D.nz] = T[indb];
            }
        }
        
        io.writeBinaryFloat(eikonalPath + "eikonal_nz" + std::to_string(m3D.nz) + "_nx" + std::to_string(m3D.nx) + "_ny" + std::to_string(m3D.ny) + "_shot_" + std::to_string(shotId+1) + ".bin", TT, m3D.nPoints);

        delete[] TT;
    }
}

void::Eikonal::writeFirstArrivals()
{
    if (exportFirstArrivals) 
    {
        Utils::point3D p;    
        
        float * firstArrivals = new float[g3D.nodes.n]();
        
        for (int r = 0; r < g3D.nodes.n; r++)
        {
            p.x = g3D.nodes.x[r];
            p.y = g3D.nodes.y[r];
            p.z = g3D.nodes.z[r];

            firstArrivals[r] = utils.triLinearInterpolation(p, m3D, T);
            // firstArrivals[r] = utils.triCubicInterpolation(p, m3D, T);            
        }

        io.writeBinaryFloat(arrivalsPath + "times_nr" + std::to_string(g3D.nodes.n) + "_shot_" + std::to_string(shotId+1) + ".bin", firstArrivals, g3D.nodes.n);

        delete[] firstArrivals;
    }
}

