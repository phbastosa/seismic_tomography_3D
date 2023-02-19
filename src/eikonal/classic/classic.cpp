# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "classic.hpp"

float Classic::min(float v1, float v2) { return !(v2 < v1) ? v1 : v2; }

void Classic::prepare_volumes()
{



}

void Classic::solve()
{
    int nPoints = model.total_samples_b;

    float * K = new float[nPoints]();    
    float * nT = new float[nPoints]();    
    float * nK = new float[nPoints]();  

    int nxx = model.x_samples_b;
    int nyy = model.y_samples_b;
    int nzz = model.z_samples_b;

    int nb = model.boundary_samples;

    float dx = model.x_spacing;
    float dy = model.y_spacing;
    float dz = model.z_spacing;

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

            travel_time[sId] = dist * slowness[sId];
            nT[sId] = travel_time[sId]; 
        }
        else
        {
            nT[index] = 1e6f;
            travel_time[index] = 1e6f;
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
    # pragma acc enter data copyin(this[0:1], slowness[0:nPoints])
    # pragma acc enter data copyin(this[0:1], travel_time[0:nPoints])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {  
            # pragma acc parallel loop present(slowness[0:nPoints],travel_time[0:nPoints],K[0:nPoints],nT[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                if (K[index] == 1.0f)
                {
                    int k = (int) (index / (nxx*nzz));             // y direction
                    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
                    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

                    if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
                    {
                        float h = dx;
                        float lowest = travel_time[index];
                        float Tijk, T1, T2, Sref, M, N, P, Q, hs2; 

                        /* 1D operator head wave: i,j-1,k -> i,j,k (x direction) */
                        Tijk = travel_time[index - nzz] + h * min(slowness[index - nzz], 
                                                                     min(slowness[index - 1 - nzz], 
                                                                     min(slowness[index - nzz - nxx*nzz], slowness[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j+1,k -> i,j,k (x direction) */
                        Tijk = travel_time[index + nzz] + h * min(slowness[index], 
                                                                     min(slowness[index - 1], 
                                                                     min(slowness[index - nxx*nzz], slowness[index - 1 - nxx*nzz])));
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k-1 -> i,j,k (y direction) */
                        Tijk = travel_time[index - nxx*nzz] + h * min(slowness[index - nxx*nzz], 
                                                                         min(slowness[index - nzz - nxx*nzz], 
                                                                         min(slowness[index - 1 - nxx*nzz], slowness[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k+1 -> i,j,k (y direction) */
                        Tijk = travel_time[index + nxx*nzz] + h * min(slowness[index],
                                                                         min(slowness[index - 1], 
                                                                         min(slowness[index - nzz], slowness[index - 1 - nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i-1,j,k -> i,j,k (z direction) */
                        Tijk = travel_time[index - 1] + h * min(slowness[index - 1], 
                                                                   min(slowness[index - 1 - nzz], 
                                                                   min(slowness[index - 1 - nxx*nzz], slowness[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i+1,j,k -> i,j,k (z direction) */
                        Tijk = travel_time[index + 1] + h * min(slowness[index], 
                                                                   min(slowness[index - nzz], 
                                                                   min(slowness[index - nxx*nzz], slowness[index - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;
                    
                        /* 1D operator diffraction XZ plane */
                        
                        // i-1,j-1,k -> i,j,k
                        Tijk = travel_time[index - 1 - nzz] + h*sqrt2*slowness[index - 1 - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k -> i,j,k
                        Tijk = travel_time[index - 1 + nzz] + h*sqrt2*slowness[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j-1,k -> i,j,k
                        Tijk = travel_time[index + 1 - nzz] + h*sqrt2*slowness[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j+1,k -> i,j,k
                        Tijk = travel_time[index + 1 + nzz] + h*sqrt2*slowness[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction YZ plane */

                        // i-1,j,k-1 -> i,j,k
                        Tijk = travel_time[index - 1 - nxx*nzz] + h*sqrt2*slowness[index - 1 - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j,k+1 -> i,j,k
                        Tijk = travel_time[index - 1 + nxx*nzz] + h*sqrt2*slowness[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k-1 -> i,j,k
                        Tijk = travel_time[index + 1 - nxx*nzz] + h*sqrt2*slowness[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;
                        
                        // i+1,j,k+1 -> i,j,k
                        Tijk = travel_time[index + 1 + nxx*nzz] + h*sqrt2*slowness[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator diffraction XY plane */
                        
                        // i,j-1,k-1 -> i,j,k
                        Tijk = travel_time[index - nzz - nxx*nzz] + h*sqrt2*slowness[index - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j-1,k+1 -> i,j,k
                        Tijk = travel_time[index - nzz + nxx*nzz] + h*sqrt2*slowness[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k-1 -> i,j,k
                        Tijk = travel_time[index + nzz - nxx*nzz] + h*sqrt2*slowness[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i,j+1,k+1 -> i,j,k
                        Tijk = travel_time[index + nzz + nxx*nzz] + h*sqrt2*slowness[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator corner diffractions */

                        // i-1,j-1,k-1 -> i,j,k
                        Tijk = travel_time[index - 1 - nzz - nxx*nzz] + h*sqrt3*slowness[index - 1 - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j-1,k+1 -> i,j,k
                        Tijk = travel_time[index - 1 - nzz + nxx*nzz] + h*sqrt3*slowness[index - 1 - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k-1 -> i,j,k
                        Tijk = travel_time[index + 1 - nzz - nxx*nzz] + h*sqrt3*slowness[index - nzz - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j-1,k+1 -> i,j,k
                        Tijk = travel_time[index + 1 - nzz + nxx*nzz] + h*sqrt3*slowness[index - nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k-1 -> i,j,k
                        Tijk = travel_time[index - 1 + nzz - nxx*nzz] + h*sqrt3*slowness[index - 1 - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i-1,j+1,k+1 -> i,j,k
                        Tijk = travel_time[index - 1 + nzz + nxx*nzz] + h*sqrt3*slowness[index - 1]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k-1 -> i,j,k
                        Tijk = travel_time[index + 1 + nzz - nxx*nzz] + h*sqrt3*slowness[index - nxx*nzz]; 
                        if (Tijk < lowest) lowest = Tijk;

                        // i+1,j+1,k+1 -> i,j,k
                        Tijk = travel_time[index + 1 + nzz + nxx*nzz] + h*sqrt3*slowness[index]; 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 2D operator XZ plane: First Quadrant*/

                        Sref = slowness[index - 1 - nzz];

                        // i,j-1,k - i-1,j-1,k -> i,j,k
                        T1 = travel_time[index - nzz];
                        T2 = travel_time[index - 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j-1,k -> i,j,k
                        T1 = travel_time[index - 1];
                        T2 = travel_time[index - 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Second Quadrant*/                        

                        Sref = slowness[index - nzz];

                        // i,j-1,k - i+1,j-1,k -> i,j,k
                        T1 = travel_time[index - nzz];
                        T2 = travel_time[index + 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i+1,j,k - i+1,j-1,k -> i,j,k
                        T1 = travel_time[index + 1];
                        T2 = travel_time[index + 1 - nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Third Quadrant*/                        

                        Sref = slowness[index];

                        // i+1,j,k - i+1,j+1,k -> i,j,k
                        T1 = travel_time[index + 1];
                        T2 = travel_time[index + 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i+1,j+1,k -> i,j,k
                        T1 = travel_time[index + nzz];
                        T2 = travel_time[index + 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XZ plane: Fourth Quadrant*/                        

                        Sref = slowness[index - 1];

                        // i,j+1,k - i-1,j+1,k -> i,j,k
                        T1 = travel_time[index + nzz];
                        T2 = travel_time[index - 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j+1,k -> i,j,k
                        T1 = travel_time[index - 1];
                        T2 = travel_time[index - 1 + nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: First Quadrant */                        

                        Sref = slowness[index - 1 - nxx*nzz];

                        // i,j,k-1 - i-1,j,k-1 -> i,j,k
                        T1 = travel_time[index - nxx*nzz];
                        T2 = travel_time[index - 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j,k-1 -> i,j,k
                        T1 = travel_time[index - 1];
                        T2 = travel_time[index - 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Second Quadrant */                        

                        Sref = slowness[index - nxx*nzz];

                        // i,j,k-1 - i+1,j,k-1 -> i,j,k
                        T1 = travel_time[index - nxx*nzz];
                        T2 = travel_time[index + 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i+1,j,k - i+1,j,k-1 -> i,j,k
                        T1 = travel_time[index + 1];
                        T2 = travel_time[index + 1 - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Third Quadrant*/                        

                        Sref = slowness[index];

                        // i+1,j,k - i+1,j,k+1 -> i,j,k
                        T1 = travel_time[index + 1];
                        T2 = travel_time[index + 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i+1,j,k+1 -> i,j,k
                        T1 = travel_time[index + nxx*nzz];
                        T2 = travel_time[index + 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator YZ plane: Fourth Quadrant*/                        

                        Sref = slowness[index - 1];

                        // i,j,k+1 - i-1,j,k+1 -> i,j,k
                        T1 = travel_time[index + nxx*nzz];
                        T2 = travel_time[index - 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i-1,j,k - i-1,j,k+1 -> i,j,k
                        T1 = travel_time[index - 1];
                        T2 = travel_time[index - 1 + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: First Quadrant*/                        

                        Sref = slowness[index - nzz - nxx*nzz];

                        // i,j-1,k - i,j-1,k-1 -> i,j,k
                        T1 = travel_time[index - nzz];
                        T2 = travel_time[index - nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j-1,k-1 -> i,j,k
                        T1 = travel_time[index - nxx*nzz];
                        T2 = travel_time[index - nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Second Quadrant*/                        

                        Sref = slowness[index - nzz];

                        // i,j-1,k - i,j-1,k+1 -> i,j,k
                        T1 = travel_time[index - nzz];
                        T2 = travel_time[index - nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k+1 - i,j-1,k+1 -> i,j,k
                        T1 = travel_time[index + nxx*nzz];
                        T2 = travel_time[index - nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Third Quadrant*/                        

                        Sref = slowness[index];

                        // i,j,k+1 - i,j+1,k+1 -> i,j,k
                        T1 = travel_time[index + nxx*nzz];
                        T2 = travel_time[index + nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j+1,k - i,j+1,k+1 -> i,j,k
                        T1 = travel_time[index + nzz];
                        T2 = travel_time[index + nzz + nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 2D operator XY plane: Fourth Quadrant*/                        

                        Sref = slowness[index - nxx*nzz];

                        // i,j+1,k - i,j+1,k-1 -> i,j,k
                        T1 = travel_time[index + nzz];
                        T2 = travel_time[index + nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        // i,j,k-1 - i,j+1,k-1 -> i,j,k
                        T1 = travel_time[index - nxx*nzz];
                        T2 = travel_time[index + nzz - nxx*nzz];
                        if ((T1 - T2) > 0.0f)
                        {
                            if ((T1 - T2) < h*Sref/sqrt2)
                            {
                                Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                                if (Tijk < lowest) lowest = Tijk;
                            }
                        }

                        /* 3D operator - First octant: XY plane */

                        Sref = slowness[index - 1 - nzz - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k-1 */   M = travel_time[index - 1 - nzz - nxx*nzz];   
    /* i-1,j-1, k  */   N = travel_time[index - 1 - nzz];             
    /* i-1, j ,k-1 */   P = travel_time[index - 1 - nxx*nzz];       
    /* i-1, j , k  */   Q = travel_time[index - 1];                 

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

    /* i-1,j-1,k-1 */   M = travel_time[index - 1 - nzz - nxx*nzz];   
    /* i-1,j-1, k  */   N = travel_time[index - 1 - nzz];             
    /*  i ,j-1,k-1 */   P = travel_time[index - nzz - nxx*nzz];       
    /*  i ,j-1, k  */   Q = travel_time[index - nzz];                 

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

    /* i-1,j-1,k-1 */   M = travel_time[index - 1 - nzz - nxx*nzz];   
    /*  i ,j-1,k-1 */   N = travel_time[index - nzz - nxx*nzz];             
    /* i-1, j ,k-1 */   P = travel_time[index - 1 - nxx*nzz];       
    /*  i , j ,k-1 */   Q = travel_time[index - nxx*nzz];                 

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

                        Sref = slowness[index - 1 - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j+1,k-1 */   M = travel_time[index - 1 + nzz - nxx*nzz];   
    /* i-1, j ,k-1 */   N = travel_time[index - 1 - nxx*nzz];             
    /* i-1,j+1, k  */   P = travel_time[index - 1 + nzz];       
    /* i-1, j , k  */   Q = travel_time[index - 1];                 

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

    /* i-1,j+1,k-1 */   M = travel_time[index - 1 + nzz - nxx*nzz];   
    /* i-1,j+1, k  */   N = travel_time[index - 1 + nzz];             
    /*  i ,j+1,k-1 */   P = travel_time[index + nzz - nxx*nzz];       
    /*  i ,j+1, k  */   Q = travel_time[index + nzz];                 

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

    /* i-1,j+1,k-1 */   M = travel_time[index - 1 + nzz - nxx*nzz];   
    /* i-1, j ,k-1 */   N = travel_time[index - 1 - nxx*nzz];             
    /*  i ,j+1,k-1 */   P = travel_time[index + nzz - nxx*nzz];       
    /*  i , j ,k-1 */   Q = travel_time[index - nxx*nzz];                 

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

                        Sref = slowness[index - 1];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j+1,k+1 */   M = travel_time[index - 1 + nzz + nxx*nzz];   
    /* i-1,j+1, k  */   N = travel_time[index - 1 + nzz];             
    /* i-1, j ,k+1 */   P = travel_time[index - 1 + nxx*nzz];       
    /* i-1, j , k  */   Q = travel_time[index - 1];                 

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

    /* i-1,j+1,k+1 */   M = travel_time[index - 1 + nzz + nxx*nzz];   
    /*  i ,j+1,k+1 */   N = travel_time[index + nzz + nxx*nzz];             
    /* i-1,j+1, k  */   P = travel_time[index - 1 + nzz];       
    /*  i ,j+1, k  */   Q = travel_time[index + nzz];                 

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

    /* i-1,j+1,k+1 */   M = travel_time[index - 1 + nzz + nxx*nzz];   
    /* i-1, j ,k+1 */   N = travel_time[index - 1 + nxx*nzz];             
    /*  i ,j+1,k+1 */   P = travel_time[index + nzz + nxx*nzz];       
    /*  i , j ,k+1 */   Q = travel_time[index + nxx*nzz];                 

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

                        Sref = slowness[index - 1 - nzz];
                        hs2 = h*h*Sref*Sref;

    /* i-1,j-1,k+1 */   M = travel_time[index - 1 - nzz + nxx*nzz];   
    /* i-1, j ,k+1 */   N = travel_time[index - 1 + nxx*nzz];             
    /* i-1,j-1, k  */   P = travel_time[index - 1 - nzz];       
    /* i-1, j , k  */   Q = travel_time[index - 1];                 

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

    /* i-1,j-1,k+1 */   M = travel_time[index - 1 - nzz + nxx*nzz];   
    /* i-1,j-1, k  */   N = travel_time[index - 1 - nzz];             
    /*  i ,j-1,k+1 */   P = travel_time[index - nzz + nxx*nzz];       
    /*  i ,j-1, k  */   Q = travel_time[index - nzz];                 

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

    /* i-1,j-1,k+1 */   M = travel_time[index - 1 - nzz + nxx*nzz];   
    /*  i ,j-1,k+1 */   N = travel_time[index - nzz + nxx*nzz];             
    /* i-1, j ,k+1 */   P = travel_time[index - 1 + nxx*nzz];       
    /*  i , j ,k+1 */   Q = travel_time[index + nxx*nzz];                 

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

                        Sref = slowness[index - nzz - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k-1 */   M = travel_time[index + 1 - nzz - nxx*nzz];   
    /* i+1, j ,k-1 */   N = travel_time[index + 1 - nxx*nzz];             
    /* i+1,j-1, k  */   P = travel_time[index + 1 - nzz];       
    /* i+1, j , k  */   Q = travel_time[index + 1];                 

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

    /* i+1,j-1,k-1 */   M = travel_time[index + 1 - nzz - nxx*nzz];   
    /* i+1,j-1, k  */   N = travel_time[index + 1 - nzz];             
    /*  i ,j-1,k-1 */   P = travel_time[index - nzz - nxx*nzz];       
    /*  i ,j-1, k  */   Q = travel_time[index - nzz];                 

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

    /* i+1,j-1,k-1 */   M = travel_time[index + 1 - nzz - nxx*nzz];   
    /*  i ,j-1,k-1 */   N = travel_time[index - nzz - nxx*nzz];             
    /* i+1, j ,k-1 */   P = travel_time[index + 1 - nxx*nzz];       
    /*  i , j ,k-1 */   Q = travel_time[index - nxx*nzz];                 

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

                        Sref = slowness[index - nxx*nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j+1,k-1 */   M = travel_time[index + 1 + nzz - nxx*nzz];   
    /* i+1,j+1, k  */   N = travel_time[index + 1 + nzz];             
    /* i+1, j ,k-1 */   P = travel_time[index + 1 - nxx*nzz];       
    /* i+1, j , k  */   Q = travel_time[index + 1];                 

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

    /* i+1,j+1,k-1 */   M = travel_time[index + 1 + nzz - nxx*nzz];   
    /*  i ,j+1,k-1 */   N = travel_time[index + nzz - nxx*nzz];             
    /* i+1,j+1, k  */   P = travel_time[index + 1 + nzz];       
    /*  i ,j+1, k  */   Q = travel_time[index + nzz];                 

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

    /* i+1,j+1,k-1 */   M = travel_time[index + 1 + nzz - nxx*nzz];   
    /* i+1, j ,k-1 */   N = travel_time[index + 1 - nxx*nzz];             
    /*  i ,j+1,k-1 */   P = travel_time[index + nzz - nxx*nzz];       
    /*  i , j ,k-1 */   Q = travel_time[index - nxx*nzz];                 

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
                        
                        Sref = slowness[index - nzz];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j-1,k+1 */   M = travel_time[index + 1 - nzz + nxx*nzz];   
    /* i+1,j-1, k  */   N = travel_time[index + 1 - nzz];             
    /* i+1, j ,k+1 */   P = travel_time[index + 1 + nxx*nzz];       
    /* i+1, j , k  */   Q = travel_time[index + 1];                 

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

    /* i+1,j-1,k+1 */   M = travel_time[index + 1 - nzz + nxx*nzz];   
    /*  i ,j-1,k+1 */   N = travel_time[index - nzz + nxx*nzz];             
    /* i+1,j-1, k  */   P = travel_time[index + 1 - nzz];       
    /*  i ,j-1, k  */   Q = travel_time[index - nzz];                 

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

    /* i+1,j-1,k+1 */   M = travel_time[index + 1 - nzz + nxx*nzz];   
    /* i+1, j ,k+1 */   N = travel_time[index + 1 + nxx*nzz];             
    /*  i ,j-1,k+1 */   P = travel_time[index - nzz + nxx*nzz];       
    /*  i , j ,k+1 */   Q = travel_time[index + nxx*nzz];                 

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

                        Sref = slowness[index];
                        hs2 = h*h*Sref*Sref;

    /* i+1,j+1,k+1 */   M = travel_time[index + 1 + nzz + nxx*nzz];   
    /* i+1, j ,k+1 */   N = travel_time[index + 1 + nxx*nzz];             
    /* i+1,j+1, k  */   P = travel_time[index + 1 + nzz];       
    /* i+1, j , k  */   Q = travel_time[index + 1];                 

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

    /* i+1,j+1,k+1 */   M = travel_time[index + 1 + nzz + nxx*nzz];   
    /* i+1,j+1, k  */   N = travel_time[index + 1 + nzz];             
    /*  i ,j+1,k+1 */   P = travel_time[index + nzz + nxx*nzz];       
    /*  i ,j+1, k  */   Q = travel_time[index + nzz];                 

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

    /* i+1,j+1,k+1 */   M = travel_time[index + 1 + nzz + nxx*nzz];   
    /*  i ,j+1,k+1 */   N = travel_time[index + nzz + nxx*nzz];             
    /* i+1, j ,k+1 */   P = travel_time[index + 1 + nxx*nzz];       
    /*  i , j ,k+1 */   Q = travel_time[index + nxx*nzz];                 

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
                        if (lowest == travel_time[index]) K[index] = 0.0f;

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
                    int k = (int) (index / (nxx*nzz));             // y direction
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

            # pragma acc parallel loop present(travel_time[0:nPoints],nT[0:nPoints],K[0:nPoints],nK[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                travel_time[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(K[0:nPoints], this[0:1])
    # pragma acc exit data delete(nT[0:nPoints], this[0:1])
    # pragma acc exit data delete(nK[0:nPoints], this[0:1])
    # pragma acc exit data copyout(slowness[0:nPoints], this[0:1])
    # pragma acc exit data copyout(travel_time[0:nPoints], this[0:1])

    delete[] K;
    delete[] nT;
    delete[] nK;
}