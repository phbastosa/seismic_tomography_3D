# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <cmath>
# include <string>
# include <algorithm>

# include "utils.hpp"
# include "model.hpp"
# include "geometry.hpp"

/*  */
typedef struct                    
{
    int sgntz; int sgntx;   // 
    int sgnty; int sgnvy;   //
    int sgnvz; int sgnvx;   //

    int i, j, k;            //

    float dzi, dxi, dyi;    //
    float dz2i, dx2i, dy2i; //
    float dz2dx2, dz2dy2;   //
    float dx2dy2, dsum;     //
    
} FSM;
    
/* */
void writeTravelTimes(float * T, int nx, int ny, int nz, int sId, std::string folder)
{
    writeBinaryFloat(folder + "eikonal_nz" + std::to_string(nz) + "_nx" + std::to_string(nx) + "_ny" + std::to_string(ny) + "_shot_" + std::to_string(sId+1) + ".bin", T, nx*ny*nz);
}
    
/* */
void writeFirstArrivals(float * T, Position *nodes, int nx, int ny, int nz, float dx, float dy, float dz, int sId, std::string folder)
{
    float * firstArrivals = new float[nodes->n]();
        
    for (int r = 0; r < nodes->n; r++)
    {
        float x = nodes->x[r];
        float y = nodes->y[r];
        float z = nodes->z[r];

        int rIdx = (int) (nodes->x[r] / dx);
        int rIdy = (int) (nodes->y[r] / dy);
        int rIdz = (int) (nodes->z[r] / dz);

        int rid = rIdz + rIdx*nz + rIdy*nx*nz;

        float x0 = rIdx*dx;
        float y0 = rIdy*dy;
        float z0 = rIdz*dz;

        float x1 = rIdx*dx + dx;
        float y1 = rIdy*dy + dy;
        float z1 = rIdz*dz + dz;

        float c000 = T[rid];
        float c001 = T[rid + 1];
        float c100 = T[rid + nz]; 
        float c010 = T[rid + nx*nz]; 
        float c101 = T[rid + 1 + nz]; 
        float c011 = T[rid + 1 + nx*nz]; 
        float c110 = T[rid + nz + nx*nz]; 
        float c111 = T[rid + 1 + nz + nx*nz];

        firstArrivals[r] = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);        
    }

    writeBinaryFloat(folder + "times_nr" + std::to_string(nodes->n) + "_shot_" + std::to_string(sId+1) + ".bin", firstArrivals, nodes->n);

    delete[] firstArrivals;
}

/* */
void innerSweep(FSM fsm, float * T, float * S, int nx, int ny, int nz, float dx, float dy, float dz)
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = fsm.i - fsm.sgnvz; 
    int j1 = fsm.j - fsm.sgnvx; 
    int k1 = fsm.k - fsm.sgnvy;

    // Get local times of surrounding points
    float tv = T[(fsm.i - fsm.sgntz) + fsm.j*nz + fsm.k*nx*nz];
    float te = T[fsm.i + (fsm.j - fsm.sgntx)*nz + fsm.k*nx*nz];
    float tn = T[fsm.i + fsm.j*nz + (fsm.k - fsm.sgnty)*nx*nz];
    float tev = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*nz + fsm.k*nx*nz];
    float ten = T[fsm.i + (fsm.j - fsm.sgntx)*nz + (fsm.k - fsm.sgnty)*nx*nz];
    float tnv = T[(fsm.i - fsm.sgntz) + fsm.j*nz + (fsm.k - fsm.sgnty)*nx*nz];
    float tnve = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*nz + (fsm.k - fsm.sgnty)*nx*nz];     

    int ijk = fsm.i + fsm.j*nz + fsm.k*nx*nz;

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dz * min4(S[i1 + imax(fsm.j-1,1)*nz   + imax(fsm.k-1,1)*nx*nz], 
                          S[i1 + imax(fsm.j-1,1)*nz   + imin(fsm.k,ny-1)*nx*nz],
                          S[i1 + imin(fsm.j,nx-1)*nz + imax(fsm.k-1,1)*nx*nz], 
                          S[i1 + imin(fsm.j,nx-1)*nz + imin(fsm.k,ny-1)*nx*nz]);

    // X direction
    t1D2 = te + dx * min4(S[imax(fsm.i-1,1)   + j1*nz + imax(fsm.k-1,1)*nx*nz], 
                          S[imin(fsm.i,nz-1) + j1*nz + imax(fsm.k-1,1)*nx*nz],
                          S[imax(fsm.i-1,1)   + j1*nz + imin(fsm.k,ny-1)*nx*nz], 
                          S[imin(fsm.i,nz-1) + j1*nz + imin(fsm.k,ny-1)*nx*nz]);

    // Y direction
    t1D3 = tn + dy * min4(S[imax(fsm.i-1,1)   + imax(fsm.j-1,1)*nz   + k1*nx*nz], 
                          S[imax(fsm.i-1,1)   + imin(fsm.j,nx-1)*nz + k1*nx*nz],
                          S[imin(fsm.i,nz-1) + imax(fsm.j-1,1)*nz   + k1*nx*nz], 
                          S[imin(fsm.i,nz-1) + imin(fsm.j,nx-1)*nz + k1*nx*nz]);

    t1D = min3(t1D1,t1D2,t1D3);

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + j1*nz + imax(fsm.k-1,1)*nx*nz], S[i1 + j1*nz + imin(fsm.k, ny-1)*nx*nz]);
    
    if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*fsm.dz2i + ta*fsm.dx2i) + sqrtf(4.0f*Sref*Sref*(fsm.dz2i + fsm.dx2i) - fsm.dz2i*fsm.dx2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + imax(fsm.j-1,1)*nz + k1*nx*nz], S[i1 + imin(fsm.j,nx-1)*nz + k1*nx*nz]);

    if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*fsm.dz2i + tb*fsm.dy2i) + sqrtf(4.0f*Sref*Sref*(fsm.dz2i + fsm.dy2i) - fsm.dz2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[imax(fsm.i-1,1) + j1*nz + k1*nx*nz],S[imin(fsm.i,nz-1) + j1*nz + k1*nx*nz]);

    if((te < tn + dy*Sref) && (tn < te + dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*fsm.dx2i + tb*fsm.dy2i) + sqrtf(4.0f*Sref*Sref*(fsm.dx2i + fsm.dy2i) - fsm.dx2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dx2i + fsm.dy2i);
    }

    t2D = min3(t2D1,t2D2,t2D3);

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[i1 + j1*nz + k1*nx*nz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (min(t1D,t2D) > max3(tv,te,tn))
    {
        t2 = 9.0f*Sref*Sref*fsm.dsum;
        
        t3 = fsm.dz2dx2*(ta - tb)*(ta - tb) + fsm.dz2dy2*(tb - tc)*(tb - tc) + fsm.dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*fsm.dz2i + ta*fsm.dx2i + tc*fsm.dy2i;        
            
            t3D = (t1 + sqrtf(t2 - t3)) / fsm.dsum;
        }
    }
   
    T[ijk] = min4(T[ijk],t1D,t2D,t3D);
}

/* */
void initSweep(FSM fsm, float * T, float * S, int nx, int ny, int nz, float dx, float dy, float dz, int sIdx, int sIdy, int sIdz)
{
    // First sweeping: Top->Bottom; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = 1; 
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = imax(1, sIdy); fsm.k < ny; fsm.k++)
    {
        for (fsm.j = imax(1, sIdx); fsm.j < nx; fsm.j++)
        {
            for (fsm.i = imax(1, sIdz); fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = imax(1, sIdy); fsm.k < ny; fsm.k++)
    {
        for (fsm.j = imax(1, sIdx); fsm.j < nx; fsm.j++)
        {
            for (fsm.i = sIdz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = sIdy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = imax(1, sIdx); fsm.j < nx; fsm.j++)
        {
            for (fsm.i = imax(1, sIdz); fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = sIdy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = imax(1, sIdx); fsm.j < nx; fsm.j++)
        {
            for (fsm.i = sIdz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = imax(1, sIdy); fsm.k < ny; fsm.k++)
    {
        for (fsm.j = sIdx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = imax(1, sIdz); fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = imax(1, sIdy); fsm.k < ny; fsm.k++)
    {
        for (fsm.j = sIdx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = sIdz + 1; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = sIdy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = sIdx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = imax(1, sIdz); fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = sIdy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = sIdx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = sIdz + 1; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }
}

/* */
void fullSweep(FSM fsm, float * T, float * S, int nx, int ny, int nz, float dx, float dy, float dz)
{
    // First sweeping: Top->Bottom; West->East; South->North 
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = 1; 
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < ny; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < nx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < ny; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < nx; fsm.j++)
        {
            for (fsm.i = nz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = ny - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < nx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = ny - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < nx; fsm.j++)
        {
            for (fsm.i = nz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < ny; fsm.k++)
    {
        for (fsm.j = nx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < ny; fsm.k++)
    {
        for (fsm.j = nx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = nz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = ny - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = nx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < nz; fsm.i++)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = ny - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = nx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = nz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);
            }
        }
    }
}

/* */
float * noble(float * V, float sx, float sy, float sz, int nx, int ny, int nz, float dx, float dy, float dz)
{
    FSM fsm;

    int nPoints = nx*ny*nz;

    int sIdx = (int)(sx / dx);
    int sIdy = (int)(sy / dy);
    int sIdz = (int)(sz / dz);

    float * T = new float[nPoints];
    float * S = new float[nPoints];

    for (int index = 0; index < nPoints; index++)
    {
        T[index] = 1e6;
        S[index] = 1.0f / V[index];
    }

    int sId = sIdz + sIdx*nz + sIdy*nx*nz;

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));

    T[sId + nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    T[sId - nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    
    T[sId + nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    T[sId - nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    
    T[sId + 1 + nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId + 1 - nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    T[sId - 1 + nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId - 1 - nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf(sIdy*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    
    T[sId + 1 + nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId + 1 - nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId - 1 + nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    T[sId - 1 - nx*nz] = S[sId] * sqrtf(powf(sIdx*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    
    T[sId + nz + nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    T[sId + nz - nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    T[sId - nz + nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    T[sId - nz - nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf(sIdz*dz - sz, 2.0f));
    
    T[sId + 1 + nz + nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId + 1 - nz + nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId + 1 + nz - nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));
    T[sId + 1 - nz - nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz+1)*dz - sz, 2.0f));

    T[sId - 1 + nz + nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    T[sId - 1 - nz + nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy+1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    T[sId - 1 + nz - nx*nz] = S[sId] * sqrtf(powf((sIdx+1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));
    T[sId - 1 - nz - nx*nz] = S[sId] * sqrtf(powf((sIdx-1)*dx - sx, 2.0f) + powf((sIdy-1)*dy - sy, 2.0f) + powf((sIdz-1)*dz - sz, 2.0f));

    fsm.dzi = 1.0f / dz;
    fsm.dxi = 1.0f / dx;
    fsm.dyi = 1.0f / dy;
    fsm.dz2i = 1.0f / (dz*dz);
    fsm.dx2i = 1.0f / (dx*dx);
    fsm.dy2i = 1.0f / (dy*dy);
    fsm.dz2dx2 = fsm.dz2i * fsm.dx2i;
    fsm.dz2dy2 = fsm.dz2i * fsm.dy2i;
    fsm.dx2dy2 = fsm.dx2i * fsm.dy2i;
    fsm.dsum = fsm.dz2i + fsm.dx2i + fsm.dy2i;

    initSweep(fsm,T,S,nx,ny,nz,dx,dy,dz,sIdx,sIdy,sIdz);
    
    fullSweep(fsm,T,S,nx,ny,nz,dx,dy,dz);

    delete[] S;

    return T;
}

/* */
float * podvin(float * V, float sx, float sy, float sz, int nx, int ny, int nz, float dx, float dy, float dz)
{    
    int nb = 2;

    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    int nPoints = nxx*nyy*nzz;

    float * T = new float[nPoints]();
    float * S = new float[nPoints]();
    float * K = new float[nPoints]();
    float * nT = new float[nPoints]();
    float * nK = new float[nPoints]();

    int sIdx = (int)(sx / dx) + nb;
    int sIdy = (int)(sy / dy) + nb;
    int sIdz = (int)(sz / dz) + nb;

    int sId = sIdz + sIdx*nzz + sIdy*nxx*nzz;

    V = expandModel(V,nx,ny,nz,nb);

    for (int index = 0; index < nPoints; index++)
    {
        if (index == sId)
        {
            float sxGrid = floorf(sx / dx) * dx;
            float syGrid = floorf(sy / dy) * dy;
            float szGrid = floorf(sz / dz) * dz;

            T[index] = S[index] * sqrtf(powf(sxGrid*dx - sx,2.0f) + powf(syGrid*dy - sy,2.0f) + powf(szGrid*dz - sz,2.0f));
            nT[index] = T[index];
        } 
        else
        {
            T[index] = 1e6;
            nT[index] = 1e6;
        }       

        K[index] = 0.0f;
        nK[index] = 0.0f;

        S[index] = 1.0f / V[index];
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

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(sIdy,2.0f) + powf(sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(sIdy,2.0f) + powf(sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(sIdy,2.0f) + powf(nzz - sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    float sqrt2 = sqrtf(2.0f);
    float sqrt3 = sqrtf(3.0f);

    # pragma acc enter data copyin(S[0:nPoints])
    # pragma acc enter data copyin(K[0:nPoints])
    # pragma acc enter data copyin(nT[0:nPoints])
    # pragma acc enter data copyin(nK[0:nPoints])
    # pragma acc enter data copyin(T[0:nPoints])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {  
            # pragma acc parallel loop present(S[0:nPoints],T[0:nPoints],K[0:nPoints],nT[0:nPoints])
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
                        float lowest = T[index];
                        float Tijk, T1, T2, Sref, M, N, P, Q, hs2; 

                        /* 1D operator head wave: i,j-1,k -> i,j,k (x direction) */
                        Tijk = T[index - nzz] + h*min(S[index - nzz], min(S[index - 1 - nzz], min(S[index - nzz - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j+1,k -> i,j,k (x direction) */
                        Tijk = T[index + nzz] + h*min(S[index], min(S[index - 1], min(S[index - nxx*nzz], S[index - 1 - nxx*nzz])));
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k-1 -> i,j,k (y direction) */
                        Tijk = T[index - nxx*nzz] + h*min(S[index - nxx*nzz], min(S[index - nzz - nxx*nzz], min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i,j,k+1 -> i,j,k (y direction) */
                        Tijk = T[index + nxx*nzz] + h*min(S[index],min(S[index - 1], min(S[index - nzz], S[index - 1 - nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i-1,j,k -> i,j,k (z direction) */
                        Tijk = T[index - 1] + h*min(S[index - 1], min(S[index - 1 - nzz], min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
                        if (Tijk < lowest) lowest = Tijk;

                        /* 1D operator head wave: i+1,j,k -> i,j,k (z direction) */
                        Tijk = T[index + 1] + h*min(S[index], min(S[index - nzz], min(S[index - nxx*nzz], S[index - nzz - nxx*nzz]))); 
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

            # pragma acc parallel loop present(T[0:nPoints],nT[0:nPoints],K[0:nPoints],nK[0:nPoints])
            for (int index = 0; index < nPoints; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(S[0:nPoints])
    # pragma acc exit data delete(K[0:nPoints])
    # pragma acc exit data delete(nT[0:nPoints])
    # pragma acc exit data delete(nK[0:nPoints])
    # pragma acc exit data copyout(T[0:nPoints])

    delete[] S;
    delete[] K;
    delete[] nK;
    delete[] nT;
 
    return reduceModel(T,nx,ny,nz,nb);
}

/* */ 
float * jeong(float * V, float sx, float sy, float sz, int nx, int ny, int nz, float dx, float dy, float dz)
{
    int nb = 2;

    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    int nPoints = nxx*nyy*nzz;

    float * T = new float[nPoints]();
    float * S = new float[nPoints]();
    float * K = new float[nPoints]();
    float * nT = new float[nPoints]();
    float * nK = new float[nPoints]();

    int sIdx = (int)(sx / dx) + nb;
    int sIdy = (int)(sy / dy) + nb;
    int sIdz = (int)(sz / dz) + nb;

    int sId = sIdz + sIdx*nzz + sIdy*nxx*nzz;

    V = expandModel(V,nx,ny,nz,nb);

    for (int index = 0; index < nPoints; index++)
    {
        if (index == sId)
        {
            float sxGrid = floorf(sx / dx) * dx;
            float syGrid = floorf(sy / dy) * dy;
            float szGrid = floorf(sz / dz) * dz;

            T[index] = S[index] * sqrtf(powf(sxGrid*dx - sx,2.0f) + powf(syGrid*dy - sy,2.0f) + powf(szGrid*dz - sz,2.0f));
            nT[index] = T[index];
        } 
        else
        {
            T[index] = 1e6;
            nT[index] = 1e6;
        }       

        K[index] = 0.0f;
        nK[index] = 0.0f;

        S[index] = 1.0f / V[index];
    }

    K[sId - 1] = 1.0f;
    K[sId + 1] = 1.0f;
    K[sId - nzz] = 1.0f;
    K[sId + nzz] = 1.0f;
    K[sId - nxx*nzz] = 1.0f;
    K[sId + nxx*nzz] = 1.0f;      
 
    int aux = 0;
    int nItEikonal = 0;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(sIdy,2.0f) + powf(sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(sIdy,2.0f) + powf(sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(sIdy,2.0f) + powf(nzz - sIdz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - sIdx,2.0f) + powf(nyy - sIdy,2.0f) + powf(nzz - sIdz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    # pragma acc enter data copyin(T[0:nPoints])
    # pragma acc enter data copyin(S[0:nPoints])
    # pragma acc enter data copyin(K[0:nPoints])
    # pragma acc enter data copyin(nT[0:nPoints])
    # pragma acc enter data copyin(nK[0:nPoints])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {
            # pragma acc parallel loop present(S[0:nPoints],T[0:nPoints],K[0:nPoints],nT[0:nPoints])
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
                        float a, b, c, tmp, Tijk;
                        float tlag = - 0.4f * h * S[sId];

                        a = min(T[index - nzz],T[index + nzz]);                 // Tx min        
                        b = min(T[index - nxx*nzz],T[index + nxx*nzz]); // Ty min        
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

            # pragma acc parallel loop present(nK[0:nPoints])
            for (int index = 0; index < nPoints; index++) 
                nK[index] = 0.0f;

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
    # pragma acc exit data delete(S[0:nPoints])
    # pragma acc exit data delete(K[0:nPoints])
    # pragma acc exit data delete(nT[0:nPoints])
    # pragma acc exit data delete(nK[0:nPoints])
    # pragma acc exit data copyout(T[0:nPoints])

    delete[] S;
    delete[] K;
    delete[] nK;
    delete[] nT;

    return reduceModel(T,nx,ny,nz,nb);
}

/* */
float * eikonalComputing(float * V, int nx, int ny, int nz, float dx, float dy, float dz, float sx, float sy, float sz, int type)
{
    switch (type)
    {
    case 0: // Podvin & Lecomte (1991)
        return podvin(V,sx,sy,sz,nx,ny,nz,dx,dy,dz);

    case 1: // Jeong & Witaker (2008)
        return jeong(V,sx,sy,sz,nx,ny,nz,dx,dy,dz);

    case 2: // Noble, Gesret and Belayouni (2014)
        return noble(V,sx,sy,sz,nx,ny,nz,dx,dy,dz);

    default:// Error message
        throw std::invalid_argument("Error: Invalid integer! Parameter " + std::to_string(type) + " is not a valid eikonal type.");        
    }   
}

# endif
