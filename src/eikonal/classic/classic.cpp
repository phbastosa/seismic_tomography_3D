# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "classic.hpp"

float Classic::min(float v1, float v2) { return !(v2 < v1) ? v1 : v2; }

void Classic::prepare_volumes()
{
    S = eiko_m.expand_fdm(slowness);

    T = new float[eiko_m.total_samples_b]();
    K = new float[eiko_m.total_samples_b]();    
    nT = new float[eiko_m.total_samples_b]();    
    nK = new float[eiko_m.total_samples_b]();  
}

void Classic::solve()
{
    int nxx = eiko_m.x_samples_b;
    int nyy = eiko_m.y_samples_b;
    int nzz = eiko_m.z_samples_b;

    int nPoints = eiko_m.total_samples_b;

    int nb = 1;

    float dx = eiko_m.x_spacing;
    float dy = eiko_m.y_spacing;
    float dz = eiko_m.z_spacing;

    float sx = geometry[shots_type]->shots.x[shot_id];  
    float sy = geometry[shots_type]->shots.y[shot_id];  
    float sz = geometry[shots_type]->shots.z[shot_id];  

    int sidx = (int)(sx / dx) + nb;
    int sidy = (int)(sy / dy) + nb;
    int sidz = (int)(sz / dz) + nb;

    int sId = sidz + sidx*nzz + sidy*nxx*nzz; 

    for (int index = 0; index < nPoints; index++)
    {
        if (index == sId)
        {
            float grid_sx = floorf(sx / dx) * dx; 
            float grid_sy = floorf(sy / dy) * dy; 
            float grid_sz = floorf(sz / dz) * dz; 

            float dist = sqrtf(powf(sx - grid_sx, 2.0f) + powf(sy - grid_sy, 2.0f) + powf(sz - grid_sz, 2.0f));

            T[sId] = dist * S[sId];
            nT[sId] = T[sId]; 
        }
        else
        {
            nT[index] = 1e6f;
            T[index] = 1e6f;
        }
        
        K[index] = 0.0f;
        nK[index] = 0.0f;
    }

    K[sId - 1] = 1.0f;
    K[sId + 1] = 1.0f;
    K[sId - nzz] = 1.0f;
    K[sId + nzz] = 1.0f;
    K[sId - nxx*nzz] = 1.0f;
    K[sId + nxx*nzz] = 1.0f;
    K[sId + 1 - nzz] = 1.0f;
    K[sId - 1 - nzz] = 1.0f;
    K[sId + 1 + nzz] = 1.0f;
    K[sId - 1 + nzz] = 1.0f;
    K[sId + 1 + nxx*nzz] = 1.0f;
    K[sId + 1 - nxx*nzz] = 1.0f;
    K[sId - 1 + nxx*nzz] = 1.0f;
    K[sId - 1 - nxx*nzz] = 1.0f;
    K[sId - nzz - nxx*nzz] = 1.0f;
    K[sId - nzz + nxx*nzz] = 1.0f;
    K[sId + nzz - nxx*nzz] = 1.0f;
    K[sId + nzz + nxx*nzz] = 1.0f;
    K[sId + 1 + nzz + nxx*nzz] = 1.0f;
    K[sId + 1 + nzz - nxx*nzz] = 1.0f;
    K[sId + 1 - nzz + nxx*nzz] = 1.0f;
    K[sId + 1 - nzz - nxx*nzz] = 1.0f;
    K[sId - 1 - nzz - nxx*nzz] = 1.0f;
    K[sId - 1 - nzz + nxx*nzz] = 1.0f;
    K[sId - 1 + nzz - nxx*nzz] = 1.0f;
    K[sId - 1 + nzz + nxx*nzz] = 1.0f;

    int aux = 0;
    int nItEikonal = 0;

    aux = (int)sqrtf(powf(sidx, 2.0f) + powf(sidy,2.0f) + powf(sidz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(sidy,2.0f) + powf(sidz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(sidz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(sidy,2.0f) + powf(nzz - sidz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(sidz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    float sqrt2 = sqrtf(2.0f);
    float sqrt3 = sqrtf(3.0f);

    # pragma acc enter data copyin(this[0:1], K[0:nPoints])
    # pragma acc enter data copyin(this[0:1], nT[0:nPoints])
    # pragma acc enter data copyin(this[0:1], nK[0:nPoints])
    # pragma acc enter data copyin(this[0:1], S[0:nPoints])
    # pragma acc enter data copyin(this[0:1], T[0:nPoints])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {  
            # pragma acc parallel loop present(S[0:nPoints],T[0:nPoints],K[0:nPoints],nT[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (nxx*nzz));         // y direction
                    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
                    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

                    if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
                    {
                        float h = dx;
                        float lowest = T[index];
                        float Tijk, T1, T2, Sref, M, N, P, Q, hs2; 

                        /* 1D operator head wave: i,j-1,k -> i,j,k (x direction) */
                        Tijk = T[index - nzz] + h * min(S[index - nzz], 
                                                                     min(S[index - 1 - nzz], 
                                                                     min(S[index - nzz - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j+1,k -> i,j,k (x direction) */
                        Tijk = T[index + nzz] + h * min(S[index], 
                                                                     min(S[index - 1], 
                                                                     min(S[index - nxx*nzz], S[index - 1 - nxx*nzz])));
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k-1 -> i,j,k (y direction) */
                        Tijk = T[index - nxx*nzz] + h * min(S[index - nxx*nzz], 
                                                                         min(S[index - nzz - nxx*nzz], 
                                                                         min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k+1 -> i,j,k (y direction) */
                        Tijk = T[index + nxx*nzz] + h * min(S[index],
                                                                         min(S[index - 1], 
                                                                         min(S[index - nzz], S[index - 1 - nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i-1,j,k -> i,j,k (z direction) */
                        Tijk = T[index - 1] + h * min(S[index - 1], 
                                                                   min(S[index - 1 - nzz], 
                                                                   min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i+1,j,k -> i,j,k (z direction) */
                        Tijk = T[index + 1] + h * min(S[index], 
                                                                   min(S[index - nzz], 
                                                                   min(S[index - nxx*nzz], S[index - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;
                    
                        /* 1D operator diffraction XZ plane */
                        
                        // i-1,j-1,k -> i,j,k
                        Tijk = T[index - 1 - nzz] + h*sqrt2*S[index - 1 - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k -> i,j,k
                        Tijk = T[index - 1 + nzz] + h*sqrt2*S[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j-1,k -> i,j,k
                        Tijk = T[index + 1 - nzz] + h*sqrt2*S[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j+1,k -> i,j,k
                        Tijk = T[index + 1 + nzz] + h*sqrt2*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction YZ plane */

                        // i-1,j,k-1 -> i,j,k
                        Tijk = T[index - 1 - nxx*nzz] + h*sqrt2*S[index - 1 - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j,k+1 -> i,j,k
                        Tijk = T[index - 1 + nxx*nzz] + h*sqrt2*S[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k-1 -> i,j,k
                        Tijk = T[index + 1 - nxx*nzz] + h*sqrt2*S[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k+1 -> i,j,k
                        Tijk = T[index + 1 + nxx*nzz] + h*sqrt2*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction XY plane */
                        
                        // i,j-1,k-1 -> i,j,k
                        Tijk = T[index - nzz - nxx*nzz] + h*sqrt2*S[index - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j-1,k+1 -> i,j,k
                        Tijk = T[index - nzz + nxx*nzz] + h*sqrt2*S[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k-1 -> i,j,k
                        Tijk = T[index + nzz - nxx*nzz] + h*sqrt2*S[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k+1 -> i,j,k
                        Tijk = T[index + nzz + nxx*nzz] + h*sqrt2*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator corner diffractions */

                        // i-1,j-1,k-1 -> i,j,k
                        Tijk = T[index - 1 - nzz - nxx*nzz] + h*sqrt3*S[index - 1 - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j-1,k+1 -> i,j,k
                        Tijk = T[index - 1 - nzz + nxx*nzz] + h*sqrt3*S[index - 1 - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k-1 -> i,j,k
                        Tijk = T[index + 1 - nzz - nxx*nzz] + h*sqrt3*S[index - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k+1 -> i,j,k
                        Tijk = T[index + 1 - nzz + nxx*nzz] + h*sqrt3*S[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k-1 -> i,j,k
                        Tijk = T[index - 1 + nzz - nxx*nzz] + h*sqrt3*S[index - 1 - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k+1 -> i,j,k
                        Tijk = T[index - 1 + nzz + nxx*nzz] + h*sqrt3*S[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k-1 -> i,j,k
                        Tijk = T[index + 1 + nzz - nxx*nzz] + h*sqrt3*S[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k+1 -> i,j,k
                        Tijk = T[index + 1 + nzz + nxx*nzz] + h*sqrt3*S[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 2D operator XZ plane: First Quadrant*/

                        Sref = S[index - 1 - nzz];

                        // i,j-1,k - i-1,j-1,k -> i,j,k
                        T1 = T[index - nzz];
                        T2 = T[index - 1 - nzz];
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
                        T2 = T[index - 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Second Quadrant*/                        

                        Sref = S[index - nzz];

                        // i,j-1,k - i+1,j-1,k -> i,j,k
                        T1 = T[index - nzz];
                        T2 = T[index + 1 - nzz];
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
                        T2 = T[index + 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Third Quadrant*/                        

                        Sref = S[index];

                        // i+1,j,k - i+1,j+1,k -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i+1,j+1,k -> i,j,k
                        T1 = T[index + nzz];
                        T2 = T[index + 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Fourth Quadrant*/                        

                        Sref = S[index - 1];

                        // i,j+1,k - i-1,j+1,k -> i,j,k
                        T1 = T[index + nzz];
                        T2 = T[index - 1 + nzz];
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
                        T2 = T[index - 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: First Quadrant */                        

                        Sref = S[index - 1 - nxx*nzz];

                        // i,j,k-1 - i-1,j,k-1 -> i,j,k
                        T1 = T[index - nxx*nzz];
                        T2 = T[index - 1 - nxx*nzz];
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
                        T2 = T[index - 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Second Quadrant */                        

                        Sref = S[index - nxx*nzz];

                        // i,j,k-1 - i+1,j,k-1 -> i,j,k
                        T1 = T[index - nxx*nzz];
                        T2 = T[index + 1 - nxx*nzz];
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
                        T2 = T[index + 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Third Quadrant*/                        

                        Sref = S[index];

                        // i+1,j,k - i+1,j,k+1 -> i,j,k
                        T1 = T[index + 1];
                        T2 = T[index + 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i+1,j,k+1 -> i,j,k
                        T1 = T[index + nxx*nzz];
                        T2 = T[index + 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Fourth Quadrant*/                        

                        Sref = S[index - 1];

                        // i,j,k+1 - i-1,j,k+1 -> i,j,k
                        T1 = T[index + nxx*nzz];
                        T2 = T[index - 1 + nxx*nzz];
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
                        T2 = T[index - 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: First Quadrant*/                        

                        Sref = S[index - nzz - nxx*nzz];

                        // i,j-1,k - i,j-1,k-1 -> i,j,k
                        T1 = T[index - nzz];
                        T2 = T[index - nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j-1,k-1 -> i,j,k
                        T1 = T[index - nxx*nzz];
                        T2 = T[index - nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Second Quadrant*/                        

                        Sref = S[index - nzz];

                        // i,j-1,k - i,j-1,k+1 -> i,j,k
                        T1 = T[index - nzz];
                        T2 = T[index - nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i,j-1,k+1 -> i,j,k
                        T1 = T[index + nxx*nzz];
                        T2 = T[index - nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Third Quadrant*/                        

                        Sref = S[index];

                        // i,j,k+1 - i,j+1,k+1 -> i,j,k
                        T1 = T[index + nxx*nzz];
                        T2 = T[index + nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i,j+1,k+1 -> i,j,k
                        T1 = T[index + nzz];
                        T2 = T[index + nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Fourth Quadrant*/                        

                        Sref = S[index - nxx*nzz];

                        // i,j+1,k - i,j+1,k-1 -> i,j,k
                        T1 = T[index + nzz];
                        T2 = T[index + nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j+1,k-1 -> i,j,k
                        T1 = T[index - nxx*nzz];
                        T2 = T[index + nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 3D operator - First octant: XY plane */

                        Sref = S[index - 1 - nzz - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k-1 */   M = T[index - 1 - nzz - nxx*nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - nzz];             
    /* i-1, j ,k-1 */   P = T[index - 1 - nxx*nzz];       
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

    /* i-1,j-1,k-1 */   M = T[index - 1 - nzz - nxx*nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - nzz];             
    /*  i ,j-1,k-1 */   P = T[index - nzz - nxx*nzz];       
    /*  i ,j-1, k  */   Q = T[index - nzz];                 

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

    /* i-1,j-1,k-1 */   M = T[index - 1 - nzz - nxx*nzz];   
    /*  i ,j-1,k-1 */   N = T[index - nzz - nxx*nzz];             
    /* i-1, j ,k-1 */   P = T[index - 1 - nxx*nzz];       
    /*  i , j ,k-1 */   Q = T[index - nxx*nzz];                 

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

                        Sref = S[index - 1 - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j+1,k-1 */   M = T[index - 1 + nzz - nxx*nzz];   
    /* i-1, j ,k-1 */   N = T[index - 1 - nxx*nzz];             
    /* i-1,j+1, k  */   P = T[index - 1 + nzz];       
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

    /* i-1,j+1,k-1 */   M = T[index - 1 + nzz - nxx*nzz];   
    /* i-1,j+1, k  */   N = T[index - 1 + nzz];             
    /*  i ,j+1,k-1 */   P = T[index + nzz - nxx*nzz];       
    /*  i ,j+1, k  */   Q = T[index + nzz];                 

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

    /* i-1,j+1,k-1 */   M = T[index - 1 + nzz - nxx*nzz];   
    /* i-1, j ,k-1 */   N = T[index - 1 - nxx*nzz];             
    /*  i ,j+1,k-1 */   P = T[index + nzz - nxx*nzz];       
    /*  i , j ,k-1 */   Q = T[index - nxx*nzz];                 

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

    /* i-1,j+1,k+1 */   M = T[index - 1 + nzz + nxx*nzz];   
    /* i-1,j+1, k  */   N = T[index - 1 + nzz];             
    /* i-1, j ,k+1 */   P = T[index - 1 + nxx*nzz];       
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

    /* i-1,j+1,k+1 */   M = T[index - 1 + nzz + nxx*nzz];   
    /*  i ,j+1,k+1 */   N = T[index + nzz + nxx*nzz];             
    /* i-1,j+1, k  */   P = T[index - 1 + nzz];       
    /*  i ,j+1, k  */   Q = T[index + nzz];                 

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

    /* i-1,j+1,k+1 */   M = T[index - 1 + nzz + nxx*nzz];   
    /* i-1, j ,k+1 */   N = T[index - 1 + nxx*nzz];             
    /*  i ,j+1,k+1 */   P = T[index + nzz + nxx*nzz];       
    /*  i , j ,k+1 */   Q = T[index + nxx*nzz];                 

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

                        Sref = S[index - 1 - nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k+1 */   M = T[index - 1 - nzz + nxx*nzz];   
    /* i-1, j ,k+1 */   N = T[index - 1 + nxx*nzz];             
    /* i-1,j-1, k  */   P = T[index - 1 - nzz];       
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

    /* i-1,j-1,k+1 */   M = T[index - 1 - nzz + nxx*nzz];   
    /* i-1,j-1, k  */   N = T[index - 1 - nzz];             
    /*  i ,j-1,k+1 */   P = T[index - nzz + nxx*nzz];       
    /*  i ,j-1, k  */   Q = T[index - nzz];                 

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

    /* i-1,j-1,k+1 */   M = T[index - 1 - nzz + nxx*nzz];   
    /*  i ,j-1,k+1 */   N = T[index - nzz + nxx*nzz];             
    /* i-1, j ,k+1 */   P = T[index - 1 + nxx*nzz];       
    /*  i , j ,k+1 */   Q = T[index + nxx*nzz];                 

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

                        Sref = S[index - nzz - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k-1 */   M = T[index + 1 - nzz - nxx*nzz];   
    /* i+1, j ,k-1 */   N = T[index + 1 - nxx*nzz];             
    /* i+1,j-1, k  */   P = T[index + 1 - nzz];       
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

    /* i+1,j-1,k-1 */   M = T[index + 1 - nzz - nxx*nzz];   
    /* i+1,j-1, k  */   N = T[index + 1 - nzz];             
    /*  i ,j-1,k-1 */   P = T[index - nzz - nxx*nzz];       
    /*  i ,j-1, k  */   Q = T[index - nzz];                 

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

    /* i+1,j-1,k-1 */   M = T[index + 1 - nzz - nxx*nzz];   
    /*  i ,j-1,k-1 */   N = T[index - nzz - nxx*nzz];             
    /* i+1, j ,k-1 */   P = T[index + 1 - nxx*nzz];       
    /*  i , j ,k-1 */   Q = T[index - nxx*nzz];                 

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

                        Sref = S[index - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j+1,k-1 */   M = T[index + 1 + nzz - nxx*nzz];   
    /* i+1,j+1, k  */   N = T[index + 1 + nzz];             
    /* i+1, j ,k-1 */   P = T[index + 1 - nxx*nzz];       
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

    /* i+1,j+1,k-1 */   M = T[index + 1 + nzz - nxx*nzz];   
    /*  i ,j+1,k-1 */   N = T[index + nzz - nxx*nzz];             
    /* i+1,j+1, k  */   P = T[index + 1 + nzz];       
    /*  i ,j+1, k  */   Q = T[index + nzz];                 

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

    /* i+1,j+1,k-1 */   M = T[index + 1 + nzz - nxx*nzz];   
    /* i+1, j ,k-1 */   N = T[index + 1 - nxx*nzz];             
    /*  i ,j+1,k-1 */   P = T[index + nzz - nxx*nzz];       
    /*  i , j ,k-1 */   Q = T[index - nxx*nzz];                 

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
                        
                        Sref = S[index - nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k+1 */   M = T[index + 1 - nzz + nxx*nzz];   
    /* i+1,j-1, k  */   N = T[index + 1 - nzz];             
    /* i+1, j ,k+1 */   P = T[index + 1 + nxx*nzz];       
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

    /* i+1,j-1,k+1 */   M = T[index + 1 - nzz + nxx*nzz];   
    /*  i ,j-1,k+1 */   N = T[index - nzz + nxx*nzz];             
    /* i+1,j-1, k  */   P = T[index + 1 - nzz];       
    /*  i ,j-1, k  */   Q = T[index - nzz];                 

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

    /* i+1,j-1,k+1 */   M = T[index + 1 - nzz + nxx*nzz];   
    /* i+1, j ,k+1 */   N = T[index + 1 + nxx*nzz];             
    /*  i ,j-1,k+1 */   P = T[index - nzz + nxx*nzz];       
    /*  i , j ,k+1 */   Q = T[index + nxx*nzz];                 

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

    /* i+1,j+1,k+1 */   M = T[index + 1 + nzz + nxx*nzz];   
    /* i+1, j ,k+1 */   N = T[index + 1 + nxx*nzz];             
    /* i+1,j+1, k  */   P = T[index + 1 + nzz];       
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

    /* i+1,j+1,k+1 */   M = T[index + 1 + nzz + nxx*nzz];   
    /* i+1,j+1, k  */   N = T[index + 1 + nzz];             
    /*  i ,j+1,k+1 */   P = T[index + nzz + nxx*nzz];       
    /*  i ,j+1, k  */   Q = T[index + nzz];                 

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

    /* i+1,j+1,k+1 */   M = T[index + 1 + nzz + nxx*nzz];   
    /*  i ,j+1,k+1 */   N = T[index + nzz + nxx*nzz];             
    /* i+1, j ,k+1 */   P = T[index + 1 + nxx*nzz];       
    /*  i , j ,k+1 */   Q = T[index + nxx*nzz];                 

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

            # pragma acc parallel loop present(nK[0:nPoints])
            for (int index = 0; index < nPoints; index++) nK[index] = 0.0f;

            # pragma acc parallel loop present(K[0:nPoints], nK[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (nxx*nzz));         // y direction
                    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
                    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

                    if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
                    {
                        nK[index - 1] = 1.0f;
                        nK[index + 1] = 1.0f;
                        nK[index - nzz] = 1.0f;
                        nK[index + nzz] = 1.0f;
                        nK[index - nxx*nzz] = 1.0f;
                        nK[index + nxx*nzz] = 1.0f;
                        nK[index + 1 - nzz] = 1.0f;
                        nK[index - 1 - nzz] = 1.0f;
                        nK[index + 1 + nzz] = 1.0f;
                        nK[index - 1 + nzz] = 1.0f;
                        nK[index + 1 + nxx*nzz] = 1.0f;
                        nK[index + 1 - nxx*nzz] = 1.0f;
                        nK[index - 1 + nxx*nzz] = 1.0f;
                        nK[index - 1 - nxx*nzz] = 1.0f;
                        nK[index - nzz - nxx*nzz] = 1.0f;
                        nK[index - nzz + nxx*nzz] = 1.0f;
                        nK[index + nzz - nxx*nzz] = 1.0f;
                        nK[index + nzz + nxx*nzz] = 1.0f;
                        nK[index + 1 + nzz + nxx*nzz] = 1.0f;
                        nK[index + 1 + nzz - nxx*nzz] = 1.0f;
                        nK[index + 1 - nzz + nxx*nzz] = 1.0f;
                        nK[index + 1 - nzz - nxx*nzz] = 1.0f;
                        nK[index - 1 - nzz - nxx*nzz] = 1.0f;
                        nK[index - 1 - nzz + nxx*nzz] = 1.0f;
                        nK[index - 1 + nzz - nxx*nzz] = 1.0f;
                        nK[index - 1 + nzz + nxx*nzz] = 1.0f;
                    }
                }
            }

            # pragma acc parallel loop present(T[0:nPoints],nT[0:nPoints],K[0:nPoints],nK[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(K[0:nPoints], this[0:1])
    # pragma acc exit data delete(nT[0:nPoints], this[0:1])
    # pragma acc exit data delete(nK[0:nPoints], this[0:1])
    # pragma acc exit data delete(S[0:nPoints], this[0:1])
    # pragma acc exit data copyout(T[0:nPoints], this[0:1])

    travel_time = eiko_m.reduce_fdm(T);
}

void Classic::destroy()
{
    delete[] S;
    delete[] T;
    delete[] K;
    delete[] nT;
    delete[] nK;
}