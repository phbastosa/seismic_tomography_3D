# ifndef EIKONAL3D_H_DEFINED
# define EIKONAL3D_H_DEFINED

float * importBinaryFloat(int nPoints, char filename[])
{
    float * array = new float[nPoints]; 

    FILE * read = fopen((const char *) filename,"rb");
    fread(array,sizeof(float),nPoints,read);
    fclose(read);

    return array;
}

void exportBinaryFloat(float * array, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(array, sizeof(float), nPoints, write);
    fclose(write);
}

void readParameters(int * nx,int * ny,int * nz,float * dx,float * dy,float * dz,char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",ny); fscanf(arq,"%i",nz); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dy); fscanf(arq,"%f",dz); 
    } 
    fclose(arq);
} 

float * extractContent3D(float * arrayExpanded, int nx, int ny, int nz, int nabc)
{
    float * array = new float[nx*ny*nz];

    for (int index = 0; index < nx*ny*nz; index++)  
    {
        int kk = index / (nx*ny);
        int jj = index % nx;
        int ii = index % (nx*ny) / nx;

        if ((kk >= 0) && (kk < nz) && (jj >= 0) && (jj < nx) && (ii >= 0) && (ii < ny))
        {
            array[index] = arrayExpanded[(kk+nabc)*(nx+2*nabc)*(ny+2*nabc) + (ii+nabc)*(nx+2*nabc) + (jj+nabc)];
        }
    }

    return array;
}

float * expandModel3D(float * model,int nx, int ny, int nz, int nabc)
{
    int nxx = nx + 2*nabc;
    int nyy = ny + 2*nabc;
    int nzz = nz + 2*nabc;

    float * modelExp = new float[nxx*nyy*nzz];

    // Centering
    for (int z = nabc; z < nzz-nabc; z++)
    {
        for (int y = nabc; y < nyy-nabc; y++)
        {
            for (int x = nabc; x < nxx-nabc; x++)
            {
                modelExp[z*nxx*nyy + y*nxx + x] = model[(z-nabc)*nx*ny + (y-nabc)*nx + x-nabc];
            }
        }
    }

    // Z direction
    for (int z = 0; z < nabc; z++)
    {
        for (int y = nabc; y < nyy-nabc; y++)
        {
            for (int x = nabc; x < nxx-nabc; x++)
            {
                modelExp[z*nxx*nyy + y*nxx + x] = model[0*nx*ny + (y-nabc)*nx + x-nabc];
                modelExp[(nzz-z-1)*nxx*nyy + y*nxx + x] = model[(nz-1)*nx*ny + (y-nabc)*nx + x-nabc];
            }
        }
    }

    // X direction
    for (int x = 0; x < nabc; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = nabc; y < nyy-nabc; y++)
            {
                modelExp[z*nxx*nyy + y*nxx + x] = modelExp[z*nxx*nyy + y*nxx + nabc];
                modelExp[z*nxx*nyy + y*nxx + nxx-x-1] = modelExp[z*nxx*nyy + y*nxx + nxx-nabc-1];
            }
        }
    }

    // Y direction
    for (int y = 0; y < nabc; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                modelExp[z*nxx*nyy + y*nxx + x] = modelExp[z*nxx*nyy + nabc*nxx + x];
                modelExp[z*nxx*nyy + (nyy-y-1)*nxx + x] = modelExp[z*nxx*nyy + (nyy-nabc-1)*nxx + x];
            }
        }
    }

    return modelExp;
}

float * switch_to_ZXY(float * array,int nx,int ny,int nz)
{
    float * newArray = new float[nx*ny*nz];

    // Desire indexes:
    // i -> Top to Bottom  = z direction
    // j -> West to East   = x direction
    // k -> North to South = y direction

    for (int k = 0; k < ny; k++)
    {
        for (int j = 0; j < nx; j++)
        {
            for (int i = 0; i < nz; i++)
            {
                newArray[j + i*nx + k*nx*nz] = array[j + k*nx + i*nx*ny];
            }
        }
    }

    return newArray;
}

float * switch_to_YXZ(float * array,int nx,int ny,int nz) 
{
    float * newArray = new float[nx*ny*nz];

    // Desire indexes:
    // k -> Top to Bottom  = y direction
    // j -> West to East   = x direction
    // i -> North to South = z direction

    for (int k = 0; k < ny; k++)
    {
        for (int j = 0; j < nx; j++)
        {
            for (int i = 0; i < nz; i++)
            {
                newArray[j + k*nx + i*nx*ny] = array[j + i*nx + k*nx*nz];
            }
        }
    }

    return newArray;
}

float min(float value1, float value2)
{
    if (value1 < value2)
        return value1;
    else
        return value2;
}

float max(float value1, float value2)
{
    if (value1 > value2)
        return value1;
    else
        return value2;
}

int imin(int value1, int value2)
{
    if (value1 < value2)
        return value1;
    else
        return value2;
}

int imax(int value1, int value2)
{
    if (value1 > value2)
        return value1;
    else
        return value2;
}

float max3(float value1, float value2, float value3)
{
    float max = value1;

    if (max < value2) max = value2;
    
    if (max < value3) max = value3;

    return max;
}

float min3(float value1, float value2, float value3)
{
    float min = value1;

    if (min > value2) min = value2;
    
    if (min > value3) min = value3;

    return min;
}

float min4(float value1, float value2, float value3, float value4)
{
    float min = value1;

    if (min > value2) min = value2;
    
    if (min > value3) min = value3;

    if (min > value4) min = value4;

    return min;
}

float analytic(float So,int xi,int yi,int zi,int dx,int dy,int dz,int sx, int sy, int sz)
{
    float t = So * sqrt(powf((xi - sx)*dx,2.0f) + powf((yi - sy)*dy,2.0f) + powf((zi - sz)*dz,2.0f));

    return t;
}

void innerSweepFunction(float * T,float * S,int i,int j,int k,int nx,int ny,int nz,float dx,float dy,float dz,int sx,int sy,int sz,int sgntz,int sgntx,int sgnty,int sgnvz,int sgnvx,int sgnvy,float dzi,float dxi,float dyi,float dz2i,float dx2i,float dy2i,float dz2dx2,float dz2dy2,float dx2dy2,float dsum)
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = i - sgnvz; 
    int j1 = j - sgnvx; 
    int k1 = k - sgnvy;

    // Get local times of surrounding points
    float tv = T[j + (i-sgntz)*nx + k*nx*nz];
    float te = T[(j-sgntx) + i*nx + k*nx*nz];
    float tn = T[j + i*nx + (k-sgnty)*nx*nz];
    float tev = T[(j-sgntx) + (i-sgntz)*nx + k*nx*nz];
    float ten = T[(j-sgntx) + i*nx + (k-sgnty)*nx*nz];
    float tnv = T[j + (i-sgntz)*nx + (k-sgnty)*nx*nz];
    float tnve = T[(j-sgntx) + (i-sgntz)*nx + (k-sgnty)*nx*nz];     

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dz * min4(S[imax(j-1,1)  + i1*nx + imax(k-1,1)*nx*nz], 
                          S[imax(j-1,1)  + i1*nx + imin(k,ny-1)*nx*nz],
                          S[imin(j,nx-1) + i1*nx + imax(k-1,1)*nx*nz], 
                          S[imin(j,nx-1) + i1*nx + imin(k,ny-1)*nx*nz]);

    // X direction
    t1D2 = te + dx * min4(S[j1 + imax(i-1,1)*nx    + imax(k-1,1)*nx*nz], 
                          S[j1 + imin(i,nz-1)*nx   + imax(k-1,1)*nx*nz],
                          S[j1 + imax(i-1,1)*nx    + imin(k,ny-1)*nx*nz], 
                          S[j1 + imin(i,nz-1)*nx   + imin(k,ny-1)*nx*nz]);

    // Y direction
    t1D3 = tn + dy * min4(S[imax(j-1,1)  + imax(i-1,1)*nx + k1*nx*nz], 
                          S[imin(j,nx-1) + imax(i-1,1)*nx + k1*nx*nz],
                          S[imax(j-1,1)  + imin(i,nz-1)*nx + k1*nx*nz], 
                          S[imin(j,nx-1) + imin(i,nz-1)*nx + k1*nx*nz]);

    t1D = min3(t1D1,t1D2,t1D3);

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[j1 + i1*nx + imax(k-1,1)*nx*nz],S[j1 + i1*nx + imin(k,ny-1)*nx*nz]);
    
    if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*dz2i + ta*dx2i) + sqrt(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[imax(j-1,1) + i1*nx + k1*nx*nz],S[imin(j,nx-1) + i1*nx + k1*nx*nz]);

    if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*dz2i + tb*dy2i) + sqrt(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[j1 + imax(i-1,1)*nx + k1*nx*nz],S[j1 + imin(i,nz-1)*nx + k1*nx*nz]);

    if((te < tn + dy*Sref) && (tn < te + dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*dx2i + tb*dy2i) + sqrt(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
    }

    t2D = min3(t2D1,t2D2,t2D3);

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[j1 + i1*nx + k1*nx*nz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (min(t1D,t2D) > max3(tv,te,tn))
    {
        t2 = 9.0f*Sref*Sref*dsum;
        
        t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
            
            t3D = (t1 + sqrt(t2 - t3)) / dsum;
        }
    }
   
    T[j + i*nx + k*nx*nz] = min4(T[j + i*nx + k*nx*nz],t1D,t2D,t3D);
}

void sweep3Dinit(float * T, float * S, int nx, int ny, int nz, float dx, float dy, float dz, int sx, int sy, int sz)
{
    int sgntz; int sgntx; int sgnty;
    int sgnvz; int sgnvx; int sgnvy;

    float dzi = 1.0f / dz;
    float dxi = 1.0f / dx;
    float dyi = 1.0f / dy;
    float dz2i = 1.0f / (dz*dz);
    float dx2i = 1.0f / (dx*dx);
    float dy2i = 1.0f / (dy*dy);
    float dz2dx2 = dz2i * dx2i;
    float dz2dy2 = dz2i * dy2i;
    float dx2dy2 = dx2i * dy2i;
    float dsum = dz2i + dx2i + dy2i;

    // First sweeping: Top->Bottom; West->East; South->North
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = imax(1,sy); k < ny; k++)
    {
        for (int j = imax(1,sx); j < nx; j++)
        {
            for (int i = imax(1,sz); i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = imax(1,sy); k < ny; k++)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = imax(1,sz); i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = imax(1,sx); j < nx; j++)
        {
            for (int i = imax(1,sz); i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = imax(1,sz); i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = imax(1,sy); k < ny; k++)
    {
        for (int j = imax(1,sx); j < nx; j++)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = imax(1,sy); k < ny; k++)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = imax(1,sx); j < nx; j++)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = sy+1; k >= 0; k--)
    {
        for (int j = sx+1; j >= 0; j--)
        {
            for (int i = sz+1; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
}

void fullSweep3D(float * T, float * S, int nx, int ny, int nz, float dx, float dy, float dz, int sx, int sy, int sz)
{
    int sgntz; int sgntx; int sgnty;
    int sgnvz; int sgnvx; int sgnvy;

    float dzi = 1.0f / dz;
    float dxi = 1.0f / dx;
    float dyi = 1.0f / dy;
    float dz2i = 1.0f / (dz*dz);
    float dx2i = 1.0f / (dx*dx);
    float dy2i = 1.0f / (dy*dy);
    float dz2dx2 = dz2i * dx2i;
    float dz2dy2 = dz2i * dy2i;
    float dx2dy2 = dx2i * dy2i;
    float dsum = dz2i + dx2i + dy2i;
        
    // First sweeping: Top->Bottom; West->East; South->North 
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = 1; k < ny; k++)
    {
        for (int j = 1; j < nx; j++)
        {
            for (int i = 1; i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = 1; k < ny; k++)
    {
        for (int j = nx - 2; j >= 0; j--)
        {
            for (int i = 1; i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = ny - 2; k >= 0; k--)
    {
        for (int j = 1; j < nx; j++)
        {
            for (int i = 1; i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = ny - 2; k >= 0; k--)
    {
        for (int j = nx - 2; j >= 0; j--)
        {
            for (int i = 1; i < nz; i++)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = 1; k < ny; k++)
    {
        for (int j = 1; j < nx; j++)
        {
            for (int i = nz - 2; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = 1; k < ny; k++)
    {
        for (int j = nx - 2; j >= 0; j--)
        {
            for (int i = nz - 2; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = ny - 2; k >= 0; k--)
    {
        for (int j = 1; j < nx; j++)
        {
            for (int i = nz - 2; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    # pragma acc parallel loop present(T[0:nx*ny*nz],S[0:nx*ny*nz])
    for (int k = ny - 2; k >= 0; k--)
    {
        for (int j = nx - 2; j >= 0; j--)
        {
            for (int i = nz - 2; i >= 0; i--)
            {
                innerSweepFunction(T,S,i,j,k,nx,ny,nz,dx,dy,dz,sx,sy,sz,sgntz,sgntx,sgnty,sgnvz,sgnvx,sgnvy,dzi,dxi,dyi,dz2i,dx2i,dy2i,dz2dx2,dz2dy2,dx2dy2,dsum);
            }
        }
    }
}

float * eikonal3D(float * velocity,float xsrc, float ysrc, float zsrc,int nx, int ny, int nz, float dx,float dy, float dz)
{
    system("clear");

    int nb = 5;

    int nxx = nx + 2*nb;
    int nyy = ny + 2*nb;
    int nzz = nz + 2*nb;

    int nPoints = nxx*nyy*nzz;

    //---- Allocating variables --------------------------------------------------------------------------------------------------------

    float * S = new float[nPoints];         
    float * T = new float[nPoints];

    float * vpFull = expandModel3D(velocity,nx,ny,nz,nb);    
    float * vpEikonal = switch_to_ZXY(vpFull,nxx,nyy,nzz);

    delete[] velocity;
    delete[] vpFull;

    //----- Defining source position ---------------------------------------------------------------------------------------------------

    int sx = (int) (xsrc / dx) + nb;
    int sy = (int) (ysrc / dy) + nb;
    int sz = (int) (zsrc / dz) + nb;

    //---- Initializing timeTravels ------------------------------------------------------------------------

    int sind = sx + sz*nxx + sy*nxx*nzz;

    for (int index = 0; index < nPoints; index++) 
    {
        S[index] = 1.0f / vpEikonal[index];
        T[index] = 1e5;
    }    

    //---- Initializing points around source (FSM) ---------------------------------------------------------------------

    float So = S[sx + sz*nx + sy*nx*nz];

    T[sx + sz*nxx + sy*nxx*nzz] = analytic(So,sx,sy,sz,dx,dy,dz,xsrc,ysrc,zsrc);             
    T[sx + (sz+1)*nxx + sy*nxx*nzz] = analytic(So,sx,sy,sz+1,dx,dy,dz,xsrc,ysrc,zsrc);            
    T[(sx+1) + sz*nxx + sy*nxx*nzz] = analytic(So,sx+1,sy,sz,dx,dy,dz,xsrc,ysrc,zsrc);               
    T[sx + sz*nxx + (sy+1)*nxx*nzz] = analytic(So,sx,sy+1,sz,dx,dy,dz,xsrc,ysrc,zsrc);           
    T[(sx+1) + (sz+1)*nxx + sy*nxx*nzz] = analytic(So,sx+1,sy,sz+1,dx,dy,dz,xsrc,ysrc,zsrc);            
    T[sx + (sz+1)*nxx + (sy+1)*nxx*nzz] = analytic(So,sx,sy+1,sz+1,dx,dy,dz,xsrc,ysrc,zsrc);        
    T[(sx+1) + sz*nxx + (sy+1)*nxx*nzz] = analytic(So,sx+1,sy+1,sz,dx,dy,dz,xsrc,ysrc,zsrc);             
    T[(sx+1) + (sz+1)*nxx + (sy+1)*nxx*nzz] = analytic(So,sx+1,sy+1,sz+1,dx,dy,dz,xsrc,ysrc,zsrc);            
            
    //------- Sweeps parallelized -------------------------------------------------------- 

    sweep3Dinit(T,S,nxx,nyy,nzz,dx,dy,dz,sx,sy,sz);
    
    fullSweep3D(T,S,nxx,nyy,nzz,dx,dy,dz,sx,sy,sz);

    //--------- Exporting results ---------------------------------------------------------------------------

    float * result = switch_to_YXZ(T,nxx,nyy,nzz);

    float * travelTimes = extractContent3D(result,nx,ny,nz,nb);

    delete[] T;
    delete[] S;
    delete[] result;

    return travelTimes;
}

# endif































