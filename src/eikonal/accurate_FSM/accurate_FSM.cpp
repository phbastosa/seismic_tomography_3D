# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "accurate_FSM.hpp"

void Accurate_FSM::prepare_volumes()
{
    S = eiko_m.expand_fdm(slowness);

    T = new float[eiko_m.total_samples_b]();
}

void Accurate_FSM::inner_sweep()
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = i - sgnvz; 
    int j1 = j - sgnvx; 
    int k1 = k - sgnvy;

    // Get local times of surrounding points
    float tv = T[(i - sgntz) + j*nzz + k*nxx*nzz];
    float te = T[i + (j - sgntx)*nzz + k*nxx*nzz];
    float tn = T[i + j*nzz + (k - sgnty)*nxx*nzz];
    float tev = T[(i - sgntz) + (j - sgntx)*nzz + k*nxx*nzz];
    float ten = T[i + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];
    float tnv = T[(i - sgntz) + j*nzz + (k - sgnty)*nxx*nzz];
    float tnve = T[(i - sgntz) + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];     

    int ijk = i + j*nzz + k*nxx*nzz;

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dz * std::min(S[i1 + std::max(j-1,1)*nzz   + std::max(k-1,1)*nxx*nzz], 
                     std::min(S[i1 + std::max(j-1,1)*nzz   + std::min(k,nyy-1)*nxx*nzz], 
                     std::min(S[i1 + std::min(j,nxx-1)*nzz + std::max(k-1,1)*nxx*nzz],
                              S[i1 + std::min(j,nxx-1)*nzz + std::min(k,nyy-1)*nxx*nzz]))); 

    // X direction
    t1D2 = te + dx * std::min(S[std::max(i-1,1)   + j1*nzz + std::max(k-1,1)*nxx*nzz], 
                     std::min(S[std::min(i,nzz-1) + j1*nzz + std::max(k-1,1)*nxx*nzz],
                     std::min(S[std::max(i-1,1)   + j1*nzz + std::min(k,nyy-1)*nxx*nzz], 
                              S[std::min(i,nzz-1) + j1*nzz + std::min(k,nyy-1)*nxx*nzz])));

    // Y direction
    t1D3 = tn + dy * std::min(S[std::max(i-1,1)   + std::max(j-1,1)*nzz   + k1*nxx*nzz], 
                     std::min(S[std::max(i-1,1)   + std::min(j,nxx-1)*nzz + k1*nxx*nzz],
                     std::min(S[std::min(i,nzz-1) + std::max(j-1,1)*nzz   + k1*nxx*nzz], 
                              S[std::min(i,nzz-1) + std::min(j,nxx-1)*nzz + k1*nxx*nzz])));

    t1D = std::min(t1D1, std::min(t1D2, t1D3));

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[i1 + j1*nzz + std::max(k-1,1)*nxx*nzz], S[i1 + j1*nzz + std::min(k, nyy-1)*nxx*nzz]);
    
    if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[i1 + std::max(j-1,1)*nzz + k1*nxx*nzz], S[i1 + std::min(j,nxx-1)*nzz + k1*nxx*nzz]);

    if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*dz2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[std::max(i-1,1) + j1*nzz + k1*nxx*nzz],S[std::min(i,nzz-1) + j1*nzz + k1*nxx*nzz]);

    if((te < tn + dy*Sref) && (tn < te + dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*dx2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
    }

    t2D = std::min(t2D1, std::min(t2D2, t2D3));

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[i1 + j1*nzz + k1*nxx*nzz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (std::min(t1D,t2D) > std::max(tv,std::max(te, tn)))
    {
        t2 = 9.0f*Sref*Sref*dsum;
        
        t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
            
            t3D = (t1 + sqrtf(t2 - t3)) / dsum;
        }
    }
   
    T[ijk] = std::min(T[ijk], std::min(t1D, std::min(t2D, t3D)));
}

void Accurate_FSM::init_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Accurate_FSM::full_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North 
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Accurate_FSM::solve()
{
    nxx = eiko_m.x_samples_b;
    nyy = eiko_m.y_samples_b;
    nzz = eiko_m.z_samples_b;

    dx = eiko_m.x_spacing;
    dy = eiko_m.y_spacing;
    dz = eiko_m.z_spacing;

    float sx = geometry[shots_type]->shots.x[shot_id];
    float sy = geometry[shots_type]->shots.y[shot_id];
    float sz = geometry[shots_type]->shots.z[shot_id];

    int nb = 1;

    sidx = (int)(sx / dx) + nb; 
    sidy = (int)(sy / dy) + nb;
    sidz = (int)(sz / dz) + nb;
    
    int sId = sidz + sidx*nzz + sidy*nxx*nzz;     

    for (int index = 0; index < eiko_m.total_samples_b; index++)
        T[index] = 1e6f;

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));

    T[sId + nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1 + nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 + nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1 - nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    
    T[sId + 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    
    T[sId + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));

    T[sId - 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));

    dzi = 1.0f / dz;
    dxi = 1.0f / dx;
    dyi = 1.0f / dy;

    dz2i = 1.0f / (dz*dz);
    dx2i = 1.0f / (dx*dx);
    dy2i = 1.0f / (dy*dy);

    dz2dx2 = dz2i * dx2i;
    dz2dy2 = dz2i * dy2i;
    dx2dy2 = dx2i * dy2i;

    dsum = dz2i + dx2i + dy2i;

    init_sweep();
    full_sweep();

    travel_time = eiko_m.reduce_fdm(T);
}
