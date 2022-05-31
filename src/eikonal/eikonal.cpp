# include <cmath>
#include <algorithm>

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
    int sIdx = (int)(g3D.shots->x[shotId] / m3D.dx) + m3D.nb;
    int sIdy = (int)(g3D.shots->y[shotId] / m3D.dy) + m3D.nb;
    int sIdz = (int)(g3D.shots->z[shotId] / m3D.dz) + m3D.nb;

    int sId = sIdz + sIdx*m3D.nzz + sIdy*m3D.nxx*m3D.nzz; 

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        S[index] = 1.0f / m3D.vp[index];

        if (index == sId)
        {
            float sx = floorf(g3D.shots->x[shotId] / m3D.dx) * m3D.dx;
            float sy = floorf(g3D.shots->y[shotId] / m3D.dy) * m3D.dy;
            float sz = floorf(g3D.shots->z[shotId] / m3D.dz) * m3D.dz;

            float dist = sqrtf(powf(sx - g3D.shots->x[shotId],2.0f) + powf(sy - g3D.shots->y[shotId],2.0f) + powf(sz - g3D.shots->z[shotId],2.0f));

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

    aux = sqrtf(sIdx*sIdx + sIdy*sIdy + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
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
    int sIdx = (int)(g3D.shots->x[shotId] / m3D.dx) + m3D.nb;
    int sIdy = (int)(g3D.shots->y[shotId] / m3D.dy) + m3D.nb;
    int sIdz = (int)(g3D.shots->z[shotId] / m3D.dz) + m3D.nb;

    int sId = sIdz + sIdx*m3D.nzz + sIdy*m3D.nxx*m3D.nzz; 

    for (int index = 0; index < m3D.nPointsB; index++)
    {
        S[index] = 1.0f / m3D.vp[index];

        if (index == sId)
        {
            float sx = floorf(g3D.shots->x[shotId] / m3D.dx) * m3D.dx;
            float sy = floorf(g3D.shots->y[shotId] / m3D.dy) * m3D.dy;
            float sz = floorf(g3D.shots->z[shotId] / m3D.dz) * m3D.dz;

            float dist = sqrtf(powf(sx - g3D.shots->x[shotId],2.0f) + powf(sy - g3D.shots->y[shotId],2.0f) + powf(sz - g3D.shots->z[shotId],2.0f));

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

    aux = sqrtf(sIdx*sIdx + sIdy*sIdy + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrtf((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
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

void Eikonal::innerSweep(int i, int j, int k, int sx, int sy, int sz, int sgntz, int sgntx, int sgnty, int sgnvz, int sgnvx, int sgnvy, float dzi, float dxi, float dyi, float dz2i, float dx2i, float dy2i, float dz2dx2, float dz2dy2, float dx2dy2, float dsum)
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = i - sgnvz; 
    int j1 = j - sgnvx; 
    int k1 = k - sgnvy;

    // Get local times of surrounding points
    float tv = T[(i-sgntz) + j*m3D.nzz + k*m3D.nxx*m3D.nzz];
    float te = T[i + (j-sgntx)*m3D.nzz + k*m3D.nxx*m3D.nzz];
    float tn = T[i + j*m3D.nzz + (k-sgnty)*m3D.nxx*m3D.nzz];
    float tev = T[(i-sgntz) + (j-sgntx)*m3D.nzz + k*m3D.nxx*m3D.nzz];
    float ten = T[i + (j-sgntx)*m3D.nzz + (k-sgnty)*m3D.nxx*m3D.nzz];
    float tnv = T[(i-sgntz) + j*m3D.nzz + (k-sgnty)*m3D.nxx*m3D.nzz];
    float tnve = T[(i-sgntz) + (j-sgntx)*m3D.nzz + (k-sgnty)*m3D.nxx*m3D.nzz];     

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e6; t1D2 = 1e6; t1D3 = 1e6;     

    // Z direction
    t1D1 = tv + m3D.dz * utils.min4(S[i1 + utils.imax(j-1,1)*m3D.nzz       + utils.imax(k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[i1 + utils.imax(j-1,1)*m3D.nzz       + utils.imin(k,m3D.nyy-1)*m3D.nxx*m3D.nzz],
                                    S[i1 + utils.imin(j,m3D.nxx-1)*m3D.nzz + utils.imax(k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[i1 + utils.imin(j,m3D.nxx-1)*m3D.nzz + utils.imin(k,m3D.nyy-1)*m3D.nxx*m3D.nzz]);

    // X direction
    t1D2 = te + m3D.dx * utils.min4(S[utils.imax(i-1,1)       + j1*m3D.nzz + utils.imax(k-1,1)*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(i,m3D.nzz-1) + j1*m3D.nzz + utils.imax(k-1,1)*m3D.nxx*m3D.nzz],
                                    S[utils.imax(i-1,1)       + j1*m3D.nzz + utils.imin(k,m3D.nyy-1)*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(i,m3D.nzz-1) + j1*m3D.nzz + utils.imin(k,m3D.nyy-1)*m3D.nxx*m3D.nzz]);

    // Y direction
    t1D3 = tn + m3D.dy * utils.min4(S[utils.imax(i-1,1)       + utils.imax(j-1,1)*m3D.nzz       + k1*m3D.nxx*m3D.nzz], 
                                    S[utils.imax(i-1,1)       + utils.imin(j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz],
                                    S[utils.imin(i,m3D.nzz-1) + utils.imax(j-1,1)*m3D.nzz       + k1*m3D.nxx*m3D.nzz], 
                                    S[utils.imin(i,m3D.nzz-1) + utils.imin(j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    t1D = utils.min3(t1D1, t1D2, t1D3);

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + j1*m3D.nzz + utils.imax(k-1,1)*m3D.nxx*m3D.nzz],S[i1 + j1*m3D.nzz + utils.imin(k,m3D.nyy-1)*m3D.nxx*m3D.nzz]);
    
    if ((tv < te + m3D.dx*Sref) && (te < tv + m3D.dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*dz2i + ta*dx2i) + sqrt(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + utils.imax(j-1,1)*m3D.nzz + k1*m3D.nxx*m3D.nzz],S[i1 + utils.imin(j,m3D.nxx-1)*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    if((tv < tn + m3D.dy*Sref) && (tn < tv + m3D.dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*dz2i + tb*dy2i) + sqrt(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[utils.imax(i-1,1) + j1*m3D.nzz + k1*m3D.nxx*m3D.nzz],S[utils.imin(i,m3D.nzz-1) + j1*m3D.nzz + k1*m3D.nxx*m3D.nzz]);

    if((te < tn + m3D.dy*Sref) && (tn < te + m3D.dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*dx2i + tb*dy2i) + sqrt(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
    }

    t2D = utils.min3(t2D1,t2D2,t2D3);

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    // Sref = S[j1 + i1*m3D.nzz + k1*m3D.nxx*m3D.nzz];

    // ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    // tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    // tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    // if (min(t1D,t2D) > max3(tv,te,tn))
    // {
    //     t2 = 9.0f*Sref*Sref*dsum;
        
    //     t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
        
    //     if (t2 >= t3)
    //     {
    //         t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
            
    //         t3D = (t1 + sqrt(t2 - t3)) / dsum;
    //     }
    // }
   
    T[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz] = utils.min4(T[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz],t1D,t2D,t3D);
}

void Eikonal::initSweep()
{
    int sx = (int)(g3D.shots->x[shotId] / m3D.dx) + m3D.nb;
    int sy = (int)(g3D.shots->y[shotId] / m3D.dy) + m3D.nb;
    int sz = (int)(g3D.shots->z[shotId] / m3D.dz) + m3D.nb;

    int sgntz; int sgntx; int sgnty;
    int sgnvz; int sgnvx; int sgnvy;

    float dzi = 1.0f / m3D.dz;
    float dxi = 1.0f / m3D.dx;
    float dyi = 1.0f / m3D.dy;
    float dz2i = 1.0f / (m3D.dz*m3D.dz);
    float dx2i = 1.0f / (m3D.dx*m3D.dx);
    float dy2i = 1.0f / (m3D.dy*m3D.dy);
    float dz2dx2 = dz2i * dx2i;
    float dz2dy2 = dz2i * dy2i;
    float dx2dy2 = dx2i * dy2i;
    float dsum = dz2i + dx2i + dy2i;

    // First sweeping: Top->Bottom; West->East; South->North
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (int k = utils.imax(1,sy); k < m3D.nyy; k++)
    {
        for (int j = utils.imax(1,sx); j < m3D.nxx; j++)
        {
            for (int i = utils.imax(1,sz); i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (int k = utils.imax(1,sy); k < m3D.nyy; k++)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = utils.imax(1,sz); i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = utils.imax(1,sx); j < m3D.nxx; j++)
        {
            for (int i = utils.imax(1,sz); i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = utils.imax(1,sz); i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (int k = utils.imax(1,sy); k < m3D.nyy; k++)
    {
        for (int j = utils.imax(1,sx); j < m3D.nxx; j++)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (int k = utils.imax(1,sy); k < m3D.nyy; k++)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = utils.imax(1,sx); j < m3D.nxx; j++)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
}

void Eikonal::fullSweep()
{
    int sx = (int)(g3D.shots->x[shotId] / m3D.dx) + m3D.nb;
    int sy = (int)(g3D.shots->y[shotId] / m3D.dy) + m3D.nb;
    int sz = (int)(g3D.shots->z[shotId] / m3D.dz) + m3D.nb;

    int sgntz; int sgntx; int sgnty;
    int sgnvz; int sgnvx; int sgnvy;

    float dzi = 1.0f / m3D.dz;
    float dxi = 1.0f / m3D.dx;
    float dyi = 1.0f / m3D.dy;
    float dz2i = 1.0f / (m3D.dz*m3D.dz);
    float dx2i = 1.0f / (m3D.dx*m3D.dx);
    float dy2i = 1.0f / (m3D.dy*m3D.dy);
    float dz2dx2 = dz2i * dx2i;
    float dz2dy2 = dz2i * dy2i;
    float dx2dy2 = dx2i * dy2i;
    float dsum = dz2i + dx2i + dy2i;
        
    // First sweeping: Top->Bottom; West->East; South->North 
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (int k = 1; k < m3D.nyy; k++)
    {
        for (int j = 1; j < m3D.nxx; j++)
        {
            for (int i = 1; i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (int k = 1; k < m3D.nyy; k++)
    {
        for (int j = m3D.nxx - 2; j >= 0; j--)
        {
            for (int i = 1; i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (int k = m3D.nyy - 2; k >= 0; k--)
    {
        for (int j = 1; j < m3D.nxx; j++)
        {
            for (int i = 1; i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (int k = m3D.nyy - 2; k >= 0; k--)
    {
        for (int j = m3D.nxx - 2; j >= 0; j--)
        {
            for (int i = 1; i < m3D.nzz; i++)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (int k = 1; k < m3D.nyy; k++)
    {
        for (int j = 1; j < m3D.nxx; j++)
        {
            for (int i = m3D.nzz - 2; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (int k = 1; k < m3D.nyy; k++)
    {
        for (int j = m3D.nxx - 2; j >= 0; j--)
        {
            for (int i = m3D.nzz - 2; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (int k = m3D.nyy - 2; k >= 0; k--)
    {
        for (int j = 1; j < m3D.nxx; j++)
        {
            for (int i = m3D.nzz - 2; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (int k = m3D.nyy - 2; k >= 0; k--)
    {
        for (int j = m3D.nxx - 2; j >= 0; j--)
        {
            for (int i = m3D.nzz - 2; i >= 0; i--)
            {
                innerSweep(i,j,k,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
}

void Eikonal::nobleFSM()
{
    int sIdx = (int)(g3D.shots->x[shotId] / m3D.dx) + m3D.nb;
    int sIdy = (int)(g3D.shots->y[shotId] / m3D.dy) + m3D.nb;
    int sIdz = (int)(g3D.shots->z[shotId] / m3D.dz) + m3D.nb;

    int sId = sIdz + sIdx*m3D.nzz + sIdy*m3D.nxx*m3D.nzz; 

    for (int index = 0; index < m3D.nPointsB; index++) 
    {    
        S[index] = 1.0f / m3D.vp[index];
        T[index] = 1e-6f;
    }

    T[sId] = S[sId] * sqrt(powf(sIdx*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf(sIdy*m3D.dy - g3D.shots->y[shotId],2.0f) + powf(sIdz*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + 1] = S[sId] * sqrt(powf(sIdx*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf(sIdy*m3D.dy - g3D.shots->y[shotId],2.0f) + powf((sIdz + 1)*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + m3D.nzz] = S[sId] * sqrt(powf((sIdx+1)*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf(sIdy*m3D.dy - g3D.shots->y[shotId],2.0f) + powf(sIdz*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + m3D.nxx*m3D.nzz] = S[sId] * sqrt(powf(sIdx*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf((sIdy+1)*m3D.dy - g3D.shots->y[shotId],2.0f) + powf(sIdz*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + 1 + m3D.nzz] = S[sId] * sqrt(powf((sIdx+1)*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf(sIdy*m3D.dy - g3D.shots->y[shotId],2.0f) + powf((sIdz+1)*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + 1 + m3D.nxx*m3D.nzz] = S[sId] * sqrt(powf(sIdx*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf((sIdy+1)*m3D.dy - g3D.shots->y[shotId],2.0f) + powf((sIdz+1)*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + m3D.nzz + m3D.nxx*m3D.nzz] = S[sId] * sqrt(powf((sIdx+1)*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf((sIdy+1)*m3D.dy - g3D.shots->y[shotId],2.0f) + powf(sIdz*m3D.dz - g3D.shots->z[shotId], 2.0f));
    T[sId + 1 + m3D.nzz + m3D.nxx*m3D.nzz] = S[sId] * sqrt(powf((sIdx+1)*m3D.dx - g3D.shots->x[shotId], 2.0f) + powf((sIdy+1)*m3D.dy - g3D.shots->y[shotId],2.0f) + powf((sIdz+1)*m3D.dz - g3D.shots->z[shotId], 2.0f));

    initSweep();
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
        
        float * firstArrivals = new float[g3D.nr]();
        
        for (int r = 0; r < g3D.nr; r++)
        {
            p.x = g3D.nodes->x[r];
            p.y = g3D.nodes->y[r];
            p.z = g3D.nodes->z[r];

            firstArrivals[r] = utils.triLinearInterpolation(p, m3D, T);
            // firstArrivals[r] = utils.triCubicInterpolation(p, m3D, T);            
        }

        io.writeBinaryFloat(arrivalsPath + "times_nr" + std::to_string(g3D.nr) + "_shot_" + std::to_string(shotId+1) + ".bin", firstArrivals, g3D.nr);

        delete[] firstArrivals;
    }
}

