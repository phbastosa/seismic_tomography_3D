# include <cmath>

# include "../../essentials/inout/inout.hpp"
# include "../../essentials/model/model.hpp"
# include "../../essentials/geometry/geometry.hpp"

# include "eikonal.hpp"

void Eikonal3D::setup()
{
    std::string line;
    std::ifstream parameters(parametersFile);
    if (parameters.is_open())
    {
        while (getline(parameters, line))
        {
            std::cout<<line<<"\n";
        }
        parameters.close();
    }        
    else 
    {
        std::cout<<"Unable to open a file!"<<std::endl;
    }


}


float Eikonal3D::min(float v1, float v2)
{
    if (v1 < v2)
        return v1;
    else
        return v2;    
}

float Eikonal3D::min4(float v1, float v2, float v3, float v4)
{
    float min = v1;
    
    if (min > v2) min = v2;
    if (min > v3) min = v3;
    if (min > v4) min = v4;

    return min;
}

void Eikonal3D::allocateVolumes()
{
    T = new float[m3D.nPointsB]();    
    S = new float[m3D.nPointsB]();    
    K = new float[m3D.nPointsB]();    
    nT = new float[m3D.nPointsB]();    
    nK = new float[m3D.nPointsB]();    
}

void Eikonal3D::deleteVolumes()
{
    delete[] T;
    delete[] S;
    delete[] K;
    delete[] nT;
    delete[] nK;
}

void Eikonal3D::podvin3D()
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

    aux = sqrt(sIdx*sIdx + sIdy*sIdy + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
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
                        float Tijk, T1, T2, Sref, M, N, P, Q;
                        float hs2 = h*h*S[index]*S[index];

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

                        Sref = S[index - 1 - m3D.nzz];

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

                        Sref = S[index - m3D.nzz];

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

                        Sref = S[index];

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

                        Sref = S[index - 1];

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

                        Sref = S[index - 1 - m3D.nxx*m3D.nzz];

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

                        Sref = S[index - m3D.nxx*m3D.nzz];

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

                        Sref = S[index];

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

                        Sref = S[index - 1];

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

                        Sref = S[index - m3D.nzz - m3D.nxx*m3D.nzz];

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

                        Sref = S[index - m3D.nzz];

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

                        Sref = S[index];

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

                        Sref = S[index - m3D.nxx*m3D.nzz];

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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Second octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Third octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fourth octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fifth octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Sixth octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Seventh octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Eighth octant: XY plane */
 
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        /* 3D operator - Fourth octant: XZ plane */

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
                            Tijk = N + P - M + sqrt(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                            if (Tijk < lowest) lowest = Tijk;
                        }   

                        // QNP -> R    
                        if ((N <= Q) && (P <= Q) && 
                           ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                            if (Tijk < lowest) lowest = Tijk;
                        }

                        // NMQ -> R
                        if ((N-M >= 0) && (N-M <= Q-N) && 
                            (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                            if (Tijk < lowest) lowest = Tijk;
                        }        

                        // PMQ -> R
                        if ((P-M >= 0) && (P-M <= Q-P) && 
                            (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
                        {
                            Tijk = Q + sqrt(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
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

void Eikonal3D::fim3D()
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

    aux = sqrt(sIdx*sIdx + sIdy*sIdy + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz)); 
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt(sIdx*sIdx + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + sIdy*sIdy + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + sIdz*sIdz);
    if (aux > nItEikonal) nItEikonal = aux;

    aux = sqrt((m3D.nxx - sIdx)*(m3D.nxx - sIdx) + (m3D.nyy - sIdy)*(m3D.nyy - sIdy) + (m3D.nzz - sIdz)*(m3D.nzz - sIdz));
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
                                tmp = 0.5f * (b + c + sqrt(2.0f*h*h*S[index]*S[index] - (b - c)*(b - c)));           

                                if (tmp > b) Tijk = tmp;

                                if (Tijk > a)
                                {
                                    tmp = (a + b + c)/3.0f + sqrt(2.0f*(a*(b - a) + b*(c - b) + c*(a - c)) + 3.0f*h*h*S[index]*S[index])/3.0f;

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

void Eikonal3D::writeTravelTimes()
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
        
        io.writeBinaryFloat(eikonalPath + "eikonal_nz" + io.toString(m3D.nz) + "_nx" + io.toString(m3D.nx) + "_ny" + io.toString(m3D.ny) + "_shot_" + io.toString(shotId+1) + ".bin", TT, m3D.nPoints);

        delete[] TT;
    }
}

void::Eikonal3D::writeFirstArrivals()
{
    if (exportFirstArrivals) 
    {
        float * firstArrivals = new float[g3D.nr]();
        
        for (int r = 0; r < g3D.nr; r++)
        {
            float x = g3D.nodes->x[r];
            float y = g3D.nodes->y[r];
            float z = g3D.nodes->z[r];

            float x0 = floorf(g3D.nodes->x[r] / m3D.dx) * m3D.dx;
            float y0 = floorf(g3D.nodes->y[r] / m3D.dy) * m3D.dy;
            float z0 = floorf(g3D.nodes->z[r] / m3D.dz) * m3D.dz;

            float x1 = floorf(g3D.nodes->x[r] / m3D.dx) * m3D.dx + m3D.dx;
            float y1 = floorf(g3D.nodes->y[r] / m3D.dy) * m3D.dy + m3D.dy;
            float z1 = floorf(g3D.nodes->z[r] / m3D.dz) * m3D.dz + m3D.dz;

            int xi = (int)(g3D.nodes->x[r] / m3D.dx) + m3D.nb;    
            int yi = (int)(g3D.nodes->y[r] / m3D.dy) + m3D.nb;    
            int zi = (int)(g3D.nodes->z[r] / m3D.dz) + m3D.nb;    

            int indT = zi + xi*m3D.nzz + yi*m3D.nxx*m3D.nzz;

            float c000 = T[indT];
            float c001 = T[indT + 1];
            float c100 = T[indT + m3D.nzz]; 
            float c101 = T[indT + 1 + m3D.nzz]; 
            float c010 = T[indT + m3D.nxx*m3D.nzz]; 
            float c011 = T[indT + 1 + m3D.nxx*m3D.nzz]; 
            float c110 = T[indT + m3D.nzz + m3D.nxx*m3D.nzz]; 
            float c111 = T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz];

            float xd = (x - x0) / (x1 - x0);
            float yd = (y - y0) / (y1 - y0);
            float zd = (z - z0) / (z1 - z0);

            float c00 = c000*(1 - xd) + c100*xd;    
            float c01 = c001*(1 - xd) + c101*xd;    
            float c10 = c010*(1 - xd) + c110*xd;    
            float c11 = c011*(1 - xd) + c111*xd;    

            float c0 = c00*(1 - yd) + c10*yd;
            float c1 = c01*(1 - yd) + c11*yd;

            firstArrivals[r] = c0*(1 - zd) + c1*zd;
        }

        // InOut::writeBinaryFloat(arrivalsPath + "times_nr" + InOut::toString(g3D.nr) + "_shot_" + InOut::toString(shotId+1) + ".bin", firstArrivals, g3D.nr);

        delete[] firstArrivals;
    }
}
