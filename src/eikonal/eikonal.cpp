# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "eikonal.hpp"

float Eikonal::min(float v1, float v2) { return !(v1 > v2) ? v1 : v2; }

void Eikonal::setEikonalParameters(char * parameters)
{
    eikonalType = std::stoi(catchParameter("eikonalType", parameters));    
    
    exportTimesVolume = str2bool(catchParameter("exportTravelTimes", parameters));
    exportFirstArrivals = str2bool(catchParameter("exportFirstArrivals", parameters));
    exportRayPosition = str2bool(catchParameter("exportRayPosition", parameters));
    exportIllumination = str2bool(catchParameter("exportIllumination", parameters));

    raysFolder = catchParameter("raysFolder", parameters);
    eikonalFolder = catchParameter("eikonalFolder", parameters);
    arrivalFolder = catchParameter("arrivalFolder", parameters);
    illuminationFolder = catchParameter("illuminationFolder", parameters);

    nb = 1;
    nx = std::stoi(catchParameter("nx", parameters));
    ny = std::stoi(catchParameter("ny", parameters));
    nz = std::stoi(catchParameter("nz", parameters));
    
    initialize();

    dx = std::stof(catchParameter("dx", parameters));
    dy = std::stof(catchParameter("dy", parameters));
    dz = std::stof(catchParameter("dz", parameters));
    
    reciprocity = str2bool(catchParameter("reciprocity", parameters));
    saveGeometry = str2bool(catchParameter("saveGeometry", parameters));

    std::vector<std::string> splitted;

    shots.elevation = std::stof(catchParameter("shotsElevation", parameters));
    nodes.elevation = std::stof(catchParameter("nodesElevation", parameters));

    shots.n_xline = std::stoi(catchParameter("xShotNumber", parameters));
    shots.n_yline = std::stoi(catchParameter("yShotNumber", parameters));
    
    splitted = split(catchParameter("shotSW", parameters),',');
    set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("shotNW", parameters),',');
    set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("shotSE", parameters),',');
    set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    setGridGeometry(shots);

    nodes.n_xline = std::stoi(catchParameter("xNodeNumber", parameters));
    nodes.n_yline = std::stoi(catchParameter("yNodeNumber", parameters));
    
    splitted = split(catchParameter("nodeSW", parameters),',');
    set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("nodeNW", parameters),',');
    set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("nodeSE", parameters),',');
    set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    setGridGeometry(nodes);

    shotsPath = catchParameter("shotsPath", parameters);
    nodesPath = catchParameter("nodesPath", parameters);

    if (saveGeometry) exportPositions();
    if (reciprocity) setReciprocity();
    
    V = expand(readBinaryFloat(catchParameter("modelPath", parameters), nPoints));

    T = new float[nPointsB]();
    illumination = new float[nPointsB]();
}

void Eikonal::writeTravelTimes()
{
    if (exportTimesVolume)
    {    
        float * travelTimes = reduce(T);
        
        writeBinaryFloat(eikonalFolder + "eikonal_nz" + std::to_string(nz) + "_nx" + std::to_string(nx) + "_ny" + std::to_string(ny) + "_shot_" + std::to_string(shotId+1) + ".bin", travelTimes, nPoints);

        delete[] travelTimes;
    }
}

void Eikonal::writeIllumination()
{
    if (exportIllumination)
    {    
        float * illum = reduce(illumination);
        
        writeBinaryFloat(illuminationFolder + "illumination_nz" + std::to_string(nz) + "_nx" + std::to_string(nx) + "_ny" + std::to_string(ny) + "_shot_" + std::to_string(shotId) + ".bin", illum, nPoints);

        delete[] illum;
    }
}

void::Eikonal::writeFirstArrivals()
{
    if (exportFirstArrivals) 
    {   
        float * firstArrivals = new float[nodes.all]();
        
        for (int r = 0; r < nodes.all; r++)
        {
            float x = nodes.x[r];
            float y = nodes.y[r];
            float z = nodes.z[r];

            float x0 = floorf(x/dx)*dx;
            float y0 = floorf(y/dy)*dy;
            float z0 = floorf(z/dz)*dz;

            float x1 = floorf(x/dx)*dx + dx;
            float y1 = floorf(y/dy)*dy + dy;
            float z1 = floorf(z/dz)*dz + dz;

            int id = ((int)(z/dz) + nb) + ((int)(x/dx) + nb)*nzz + ((int)(y/dy) + nb)*nxx*nzz;

            float c000 = T[id];
            float c001 = T[id + 1];
            float c100 = T[id + nzz]; 
            float c101 = T[id + 1 + nzz]; 
            float c010 = T[id + nxx*nzz]; 
            float c011 = T[id + 1 + nxx*nzz]; 
            float c110 = T[id + nzz + nxx*nzz]; 
            float c111 = T[id + 1 + nzz + nxx*nzz];

            firstArrivals[r] = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);        
        }

        writeBinaryFloat(arrivalFolder + "times_nr" + std::to_string(nodes.all) + "_shot_" + std::to_string(shotId+1) + ".bin", firstArrivals, nodes.all);

        delete[] firstArrivals;
    }
}

void Eikonal::podvin()
{
    S = new float[nPointsB]();
    
    float * K = new float[nPointsB]();    
    float * nT = new float[nPointsB]();    
    float * nK = new float[nPointsB]();  

    shots.idx = (int)(shots.x[shotId] / dx) + nb;
    shots.idy = (int)(shots.y[shotId] / dy) + nb;
    shots.idz = (int)(shots.z[shotId] / dz) + nb;

    int sId = shots.idz + shots.idx*nzz + shots.idy*nxx*nzz; 

    for (int index = 0; index < nPointsB; index++)
    {
        S[index] = 1.0f / V[index];

        if (index == sId)
        {
            float sx = floorf(shots.x[shotId] / dx) * dx;
            float sy = floorf(shots.y[shotId] / dy) * dy;
            float sz = floorf(shots.z[shotId] / dz) * dz;

            float dist = sqrtf(powf(sx - shots.x[shotId],2.0f) + powf(sy - shots.y[shotId],2.0f) + powf(sz - shots.z[shotId],2.0f));

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

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(nzz - shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    float sqrt2 = sqrtf(2.0f);
    float sqrt3 = sqrtf(3.0f);

    # pragma acc enter data copyin(this[0:1], S[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], K[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], nT[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], nK[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], T[0:nPointsB])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {  
            # pragma acc parallel loop present(S[0:nPointsB],T[0:nPointsB],K[0:nPointsB],nT[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
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

            # pragma acc parallel loop present(nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++) nK[index] = 0.0f;

            # pragma acc parallel loop present(K[0:nPointsB], nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
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

            # pragma acc parallel loop present(T[0:nPointsB],nT[0:nPointsB],K[0:nPointsB],nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(S[0:nPointsB], this[0:1])
    # pragma acc exit data delete(K[0:nPointsB], this[0:1])
    # pragma acc exit data delete(nT[0:nPointsB], this[0:1])
    # pragma acc exit data delete(nK[0:nPointsB], this[0:1])
    # pragma acc exit data copyout(T[0:nPointsB], this[0:1])

    delete[] S;
    delete[] K;
    delete[] nT;
    delete[] nK;
}

void Eikonal::jeongFIM()
{
    S = new float[nPointsB]();

    float * K = new float[nPointsB]();    
    float * nT = new float[nPointsB]();    
    float * nK = new float[nPointsB]();  

    shots.idx = (int)(shots.x[shotId] / dx) + nb;
    shots.idy = (int)(shots.y[shotId] / dy) + nb;
    shots.idz = (int)(shots.z[shotId] / dz) + nb;

    int sId = shots.idz + shots.idx*nzz + shots.idy*nxx*nzz; 

    for (int index = 0; index < nPointsB; index++)
    {
        S[index] = 1.0f / V[index];

        if (index == sId)
        {
            float sx = floorf(shots.x[shotId] / dx) * dx;
            float sy = floorf(shots.y[shotId] / dy) * dy;
            float sz = floorf(shots.z[shotId] / dz) * dz;

            float dist = sqrtf(powf(sx - shots.x[shotId],2.0f) + powf(sy - shots.y[shotId],2.0f) + powf(sz - shots.z[shotId],2.0f));

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
    K[sId - nzz] = 1.0f;
    K[sId + nzz] = 1.0f;
    K[sId - nxx*nzz] = 1.0f;
    K[sId + nxx*nzz] = 1.0f;      
 
    int aux = 0;
    int nItEikonal = 0;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(nzz - shots.idz,2.0f)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = (int)sqrtf(powf(nxx - shots.idx,2.0f) + powf(nyy - shots.idy,2.0f) + powf(nzz - shots.idz,2.0f));
    if (aux > nItEikonal) nItEikonal = aux;

    nItEikonal += (int)(3 * nItEikonal / 2);

    # pragma acc enter data copyin(this[0:1], S[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], K[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], nT[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], nK[0:nPointsB])
    # pragma acc enter data copyin(this[0:1], T[0:nPointsB])
    {
        for (int iteration = 0; iteration < nItEikonal; iteration++)
        {
            # pragma acc parallel loop present(S[0:nPointsB],T[0:nPointsB],K[0:nPointsB],nT[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
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

            # pragma acc parallel loop present(nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++) nK[index] = 0.0f;

            # pragma acc parallel loop present(K[0:nPointsB], nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
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

            # pragma acc parallel loop present(T[0:nPointsB],nT[0:nPointsB],K[0:nPointsB],nK[0:nPointsB])
            for (int index = 0; index < nPointsB; index++)
            {
                T[index] = nT[index];
                K[index] = nK[index];
            }
        }
    }
    # pragma acc exit data delete(S[0:nPointsB], this[0:1])
    # pragma acc exit data delete(K[0:nPointsB], this[0:1])
    # pragma acc exit data delete(nT[0:nPointsB], this[0:1])
    # pragma acc exit data delete(nK[0:nPointsB], this[0:1])
    # pragma acc exit data copyout(T[0:nPointsB], this[0:1])

    delete[] S;
    delete[] K;
    delete[] nT;
    delete[] nK;
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
    float tv = T[(fsm.i - fsm.sgntz) + fsm.j*nzz + fsm.k*nxx*nzz];
    float te = T[fsm.i + (fsm.j - fsm.sgntx)*nzz + fsm.k*nxx*nzz];
    float tn = T[fsm.i + fsm.j*nzz + (fsm.k - fsm.sgnty)*nxx*nzz];
    float tev = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*nzz + fsm.k*nxx*nzz];
    float ten = T[fsm.i + (fsm.j - fsm.sgntx)*nzz + (fsm.k - fsm.sgnty)*nxx*nzz];
    float tnv = T[(fsm.i - fsm.sgntz) + fsm.j*nzz + (fsm.k - fsm.sgnty)*nxx*nzz];
    float tnve = T[(fsm.i - fsm.sgntz) + (fsm.j - fsm.sgntx)*nzz + (fsm.k - fsm.sgnty)*nxx*nzz];     

    int ijk = fsm.i + fsm.j*nzz + fsm.k*nxx*nzz;

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dz * min4(S[i1 + imax(fsm.j-1,1)*nzz   + imax(fsm.k-1,1)*nxx*nzz], 
                          S[i1 + imax(fsm.j-1,1)*nzz   + imin(fsm.k,nyy-1)*nxx*nzz],
                          S[i1 + imin(fsm.j,nxx-1)*nzz + imax(fsm.k-1,1)*nxx*nzz], 
                          S[i1 + imin(fsm.j,nxx-1)*nzz + imin(fsm.k,nyy-1)*nxx*nzz]);

    // X direction
    t1D2 = te + dx * min4(S[imax(fsm.i-1,1)   + j1*nzz + imax(fsm.k-1,1)*nxx*nzz], 
                          S[imin(fsm.i,nzz-1) + j1*nzz + imax(fsm.k-1,1)*nxx*nzz],
                          S[imax(fsm.i-1,1)   + j1*nzz + imin(fsm.k,nyy-1)*nxx*nzz], 
                          S[imin(fsm.i,nzz-1) + j1*nzz + imin(fsm.k,nyy-1)*nxx*nzz]);

    // Y direction
    t1D3 = tn + dy * min4(S[imax(fsm.i-1,1)   + imax(fsm.j-1,1)*nzz   + k1*nxx*nzz], 
                          S[imax(fsm.i-1,1)   + imin(fsm.j,nxx-1)*nzz + k1*nxx*nzz],
                          S[imin(fsm.i,nzz-1) + imax(fsm.j-1,1)*nzz   + k1*nxx*nzz], 
                          S[imin(fsm.i,nzz-1) + imin(fsm.j,nxx-1)*nzz + k1*nxx*nzz]);

    t1D = min3(t1D1,t1D2,t1D3);

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + j1*nzz + imax(fsm.k-1,1)*nxx*nzz], S[i1 + j1*nzz + imin(fsm.k, nyy-1)*nxx*nzz]);
    
    if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*fsm.dz2i + ta*fsm.dx2i) + sqrtf(4.0f*Sref*Sref*(fsm.dz2i + fsm.dx2i) - fsm.dz2i*fsm.dx2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + imax(fsm.j-1,1)*nzz + k1*nxx*nzz], S[i1 + imin(fsm.j,nxx-1)*nzz + k1*nxx*nzz]);

    if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*fsm.dz2i + tb*fsm.dy2i) + sqrtf(4.0f*Sref*Sref*(fsm.dz2i + fsm.dy2i) - fsm.dz2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dz2i + fsm.dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[imax(fsm.i-1,1) + j1*nzz + k1*nxx*nzz],S[imin(fsm.i,nzz-1) + j1*nzz + k1*nxx*nzz]);

    if((te < tn + dy*Sref) && (tn < te + dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*fsm.dx2i + tb*fsm.dy2i) + sqrtf(4.0f*Sref*Sref*(fsm.dx2i + fsm.dy2i) - fsm.dx2i*fsm.dy2i*(ta - tb)*(ta - tb))) / (fsm.dx2i + fsm.dy2i);
    }

    t2D = min3(t2D1,t2D2,t2D3);

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[i1 + j1*nzz + k1*nxx*nzz];

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

void Eikonal::initSweep()
{
    // First sweeping: Top->Bottom; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = 1; 
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = imax(1, shots.idy); fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = imax(1, shots.idx); fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = imax(1, shots.idz); fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = imax(1, shots.idy); fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = imax(1, shots.idx); fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = shots.idz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = imax(1, shots.idx); fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = imax(1, shots.idz); fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = imax(1, shots.idx); fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = shots.idz + 1; fsm.i >= 0 ; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = imax(1, shots.idy); fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = imax(1, shots.idz); fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = imax(1, shots.idy); fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = shots.idz + 1; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = imax(1, shots.idz); fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = shots.idy + 1; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = shots.idx + 1; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = shots.idz + 1; fsm.i >= 0; fsm.i--)
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

    for (fsm.k = 1; fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = 1; fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = 1; fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    fsm.sgntz = -1; fsm.sgntx = 1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx = 1; fsm.sgnvy =  0;

    for (fsm.k = nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = 1; fsm.j < nxx; fsm.j++)
        {
            for (fsm.i = nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = 1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy = 1;

    for (fsm.k = 1; fsm.k < nyy; fsm.k++)
    {
        for (fsm.j = nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    fsm.sgntz = 1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz = 1; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = 1; fsm.i < nzz; fsm.i++)
            {
                innerSweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    fsm.sgntz = -1; fsm.sgntx = -1; fsm.sgnty = -1;
    fsm.sgnvz =  0; fsm.sgnvx =  0; fsm.sgnvy =  0;

    for (fsm.k = nyy - 2; fsm.k >= 0; fsm.k--)
    {
        for (fsm.j = nxx - 2; fsm.j >= 0; fsm.j--)
        {
            for (fsm.i = nzz - 2; fsm.i >= 0; fsm.i--)
            {
                innerSweep();
            }
        }
    }
}

void Eikonal::nobleFSM()
{
    S = new float[nPointsB];
    
    shots.idx = (int)(shots.x[shotId] / dx);
    shots.idy = (int)(shots.y[shotId] / dy);
    shots.idz = (int)(shots.z[shotId] / dz);

    int sId = (shots.idz + nb) + (shots.idx + nb)*nzz + (shots.idy + nb)*nxx*nzz;     

    for (int index = 0; index < nPointsB; index++)
    {
        T[index] = 1e6f;
        S[index] = 1.0f / V[index];
    }

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));

    T[sId + nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    T[sId - nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    
    T[sId + nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    T[sId - nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    
    T[sId + 1 + nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId + 1 - nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 + nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 - nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf(shots.idy*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    
    T[sId + 1 + nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId + 1 - nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 + nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 - nxx*nzz] = S[sId] * sqrtf(powf(shots.idx*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    
    T[sId + nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    T[sId + nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    T[sId - nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    T[sId - nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf(shots.idz*dz - shots.z[shotId], 2.0f));
    
    T[sId + 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId + 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId + 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));
    T[sId + 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz+1)*dz - shots.z[shotId], 2.0f));

    T[sId - 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy+1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx+1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));
    T[sId - 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf((shots.idx-1)*dx - shots.x[shotId], 2.0f) + powf((shots.idy-1)*dy - shots.y[shotId], 2.0f) + powf((shots.idz-1)*dz - shots.z[shotId], 2.0f));

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

    shots.idx += nb;
    shots.idy += nb;
    shots.idz += nb;

    initSweep();
    fullSweep();

    delete[] S;
}

void Eikonal::eikonalComputing()
{
    switch (eikonalType)
    {
    case 0:
        podvin();
        break;

    case 1:
        jeongFIM();
        break;

    case 2:
        nobleFSM();
        break;
    
    default:
        nobleFSM();
        break;
    }
    
    rayTracing();
    writeTravelTimes();
    writeFirstArrivals();
}

void Eikonal::rayTracing()
{
    if (exportRayPosition || exportIllumination)
    {
        int sIdz = (int)(shots.z[shotId] / dz) + nb;
        int sIdx = (int)(shots.x[shotId] / dx) + nb;
        int sIdy = (int)(shots.y[shotId] / dy) + nb;

        int sId = sIdz + sIdx*nzz + sIdy*nxx*nzz;     
        
        std::vector<float> xRay, yRay, zRay, iRay;

        int im, jm, km, rId;

        float rayStep = 0.2f * (dx + dy + dz) / 3.0f;

        for (int rayId = 0; rayId < nodes.all; rayId++)
        {
            float zi = nodes.z[rayId];
            float xi = nodes.x[rayId];
            float yi = nodes.y[rayId];

            im = (int)(zi / dz) + nb; 
            jm = (int)(xi / dx) + nb; 
            km = (int)(yi / dy) + nb; 

            rId = im + jm*nzz + km*nxx*nzz;

            if (exportIllumination) illumination[rId] += rayStep;

            if (exportRayPosition)
            {
                xRay.push_back(xi);
                yRay.push_back(yi);
                zRay.push_back(zi);
                iRay.push_back((float) rayId);
            }

            while (true)
            {
                int i = (int)(zi / dz) + nb;
                int j = (int)(xi / dx) + nb;
                int k = (int)(yi / dy) + nb;

                float dTz = (T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*dz);    
                float dTx = (T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*dx);    
                float dTy = (T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*dy);

                float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

                zi -= rayStep*dTz / norm; // z ray position update   
                xi -= rayStep*dTx / norm; // x ray position update   
                yi -= rayStep*dTy / norm; // y ray position update   

                im = (int)(zi / dz) + nb; 
                jm = (int)(xi / dx) + nb; 
                km = (int)(yi / dy) + nb; 

                rId = im + jm*nzz + km*nxx*nzz;

                if (rId == sId)
                {
                    if (exportRayPosition)
                    {
                        xRay.push_back(shots.x[shotId]);
                        yRay.push_back(shots.y[shotId]);
                        zRay.push_back(shots.z[shotId]);
                        iRay.push_back((float) rayId);
                    }
                
                    if (exportIllumination)
                    {
                        float finalDist = sqrtf(powf(shots.x[shotId] - xi, 2.0f) + powf(shots.y[shotId] - yi, 2.0f) + powf(shots.z[shotId] - zi, 2.0f));        
                        illumination[rId] += finalDist;                
                    }            

                    break;
                }

                if (exportIllumination) illumination[rId] += rayStep;

                if (exportRayPosition)
                {
                    xRay.push_back(xi);
                    yRay.push_back(yi);
                    zRay.push_back(zi);
                    iRay.push_back((float) rayId);
                }
            }
        }

        if (exportRayPosition)
        {
            std::ofstream fout;

            size_t allRayPoints = iRay.size();

            std::string xRayPath = raysFolder + "xRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
            std::string yRayPath = raysFolder + "yRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
            std::string zRayPath = raysFolder + "zRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
            std::string iRayPath = raysFolder + "iRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";

            fout.open(xRayPath, std::ios::out);
            fout.write((char*)&xRay[0], xRay.size() * sizeof(xRay));
            std::cout<<"Ray x position file was written with name " + xRayPath<<std::endl;  
            fout.close();

            fout.open(yRayPath, std::ios::out);
            fout.write((char*)&yRay[0], yRay.size() * sizeof(yRay));
            std::cout<<"Ray y position file was written with name " + yRayPath<<std::endl;  
            fout.close();

            fout.open(zRayPath, std::ios::out);
            fout.write((char*)&zRay[0], zRay.size() * sizeof(zRay));
            std::cout<<"Ray z position file was written with name " + zRayPath<<std::endl;  
            fout.close();

            fout.open(iRayPath, std::ios::out);
            fout.write((char*)&iRay[0], iRay.size() * sizeof(iRay));
            std::cout<<"Ray node index file was written with name " + iRayPath<<std::endl;  
            fout.close();
        }

        std::vector<float>().swap(iRay);
        std::vector<float>().swap(xRay);
        std::vector<float>().swap(yRay);
        std::vector<float>().swap(zRay);
    }
}

