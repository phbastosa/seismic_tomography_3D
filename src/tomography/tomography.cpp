# include <cmath>
# include <limits>
# include <string>
# include <iomanip>
# include <algorithm>

# include "tomography.hpp"

# include "../eikonal/eikonal.hpp"

Tomography::Tomography(char **argv)
{
    parametersFile = argv[1];

    setup();

    dobsPath = io.catchParameter("dobsFolder", parametersFile);
    dcalPath = io.catchParameter("dcalFolder", parametersFile);
    gradPath = io.catchParameter("gradientFolder", parametersFile);

    maxIteration = std::stoi(io.catchParameter("maxIteration", parametersFile));
    tomoTolerance = std::stof(io.catchParameter("tomoTolerance", parametersFile));
    lambda = std::stof(io.catchParameter("regParam", parametersFile));

    mTomo.nx = std::stoi(io.catchParameter("nxTomo", parametersFile));
    mTomo.ny = std::stoi(io.catchParameter("nyTomo", parametersFile));
    mTomo.nz = std::stoi(io.catchParameter("nzTomo", parametersFile));

    mTomo.dx = std::stof(io.catchParameter("dxTomo", parametersFile));
    mTomo.dy = std::stof(io.catchParameter("dyTomo", parametersFile));
    mTomo.dz = std::stof(io.catchParameter("dzTomo", parametersFile));

    mTomo.nPoints = mTomo.nx * mTomo.ny * mTomo.nz;

    resPath = io.catchParameter("convergencyFolder", parametersFile);
    estModels = io.catchParameter("estimatedModelsFolder", parametersFile); 
    msv = std::stof(io.catchParameter("maxVariation", parametersFile)); 

    xMask = std::stof(io.catchParameter("xMask", parametersFile));
    yMask = std::stof(io.catchParameter("yMask", parametersFile));
    zMaskUp = std::stof(io.catchParameter("zMaskUp", parametersFile));
    zMaskDown = std::stof(io.catchParameter("zMaskDown", parametersFile));

    smoothing = utils.str2bool(io.catchParameter("smoothing", parametersFile));

    generate_dobs = utils.str2bool(io.catchParameter("generateDobs", parametersFile));

    if (generate_dobs)
    {
        arrivalsPath = dobsPath;

        m3D.vpPath = io.catchParameter("trueModelPath", parametersFile);

        m3D.readAndExpandVP();

        for (int shot = 0; shot < g3D.ns; shot++)
        {
            std::cout<<"Computing observed data shot "<<shot+1<<" of "<<g3D.ns<<"\n";

            shotId = shot;

            if (eikonalType == 1) 
                podvin();
            else 
                jeongFIM();
        }
    }

    m3D.vpPath = io.catchParameter("initModelPath", parametersFile);

    m3D.readAndExpandVP();

    iteration = 0;

    residuo.reserve(maxIteration);

    dobs = new float[g3D.ns * g3D.nr]();
    dcal = new float[g3D.ns * g3D.nr]();

    gradient = new float[mTomo.nPoints]();
    slowness = new float[mTomo.nPoints]();
}

void Tomography::infoMessage()
{
    std::cout << std::fixed << std::setprecision(1);

    system("clear");
    std::cout<<"3D first arrival tomography for OBN geometry\n\n";
    
    std::cout<<"Total x model length = "<<(m3D.nx-1)*m3D.dx<<" m\n";
    std::cout<<"Total Y model length = "<<(m3D.ny-1)*m3D.dy<<" m\n";
    std::cout<<"Total Z model length = "<<(m3D.nz-1)*m3D.dz<<" m\n\n";
    
    if (reciprocity)
        std::cout<<"Reciprocity = True\n\n";
    else
        std::cout<<"Shots reciprocity = False\n\n";

    if (iteration == maxIteration)
        std::cout<<"------- Checking final residuo ------------\n\n";
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" ------------\n\n";

    std::cout<<"Shot "<<shotId+1<<" of "<<g3D.ns<<" at position: (x = "<<g3D.shots->x[shotId]<<" m, y = "<<g3D.shots->y[shotId]<<" m, z = "<<g3D.shots->z[shotId]<<" m)\n\n";

    if (iteration > 1) std::cout<<"\nPrevious iteration residuous: "<<residuo.back()<<"\n";
} 

void Tomography::importDobs()
{   
    float * data = new float[g3D.nr]();    
    
    int ptr = 0; 
    for (int shot = 0; shot < g3D.ns; shot++)
    {
        io.readBinaryFloat(dobsPath + "times_nr" + std::to_string(g3D.nr) + "_shot_" + std::to_string(shot+1) + ".bin", data, g3D.nr);

        for (int d = ptr; d < ptr + g3D.nr; d++) dobs[d] = data[d - ptr];

        ptr += g3D.nr;
    }

    delete[] data;
}

void Tomography::setInitialModel()
{
    for (int k = 0; k < mTomo.ny; k++)
    {
        for (int j = 0; j < mTomo.nx; j++)
        {
            for (int i = 0; i < mTomo.nz; i++)
            {
                int im = (int) ((float)(i)*mTomo.dz/m3D.dz) + m3D.nb;
                int jm = (int) ((float)(j)*mTomo.dx/m3D.dx) + m3D.nb;
                int km = (int) ((float)(k)*mTomo.dy/m3D.dy) + m3D.nb;
                
                slowness[i + j*mTomo.nz + k*mTomo.nx*mTomo.nz] = 1.0f / m3D.vp[im + jm*m3D.nzz + km*m3D.nxx*m3D.nzz];    
            }
        }
    }
}

void Tomography::importDcal()
{
    float * data = new float[g3D.nr]();    
    
    int ptr = 0; 
    for (int shot = 0; shot < g3D.ns; shot++)
    {
        io.readBinaryFloat(dcalPath + "times_nr" + std::to_string(g3D.nr) + "_shot_" + std::to_string(shot+1) + ".bin", data, g3D.nr);

        for (int d = ptr; d < ptr + g3D.nr; d++) 
            dcal[d] = data[d - ptr];

        ptr += g3D.nr;
    }

    delete[] data;
}

void Tomography::fwdModeling()
{
    arrivalsPath = dcalPath;

    for (int shot = 0; shot < g3D.ns; shot++)
    {
        shotId = shot;

        if (eikonalType) 
            podvin();
        else 
            jeongFIM();

        gradientRayTracing();
    
        infoMessage();    
    }
}

void Tomography::gradientRayTracing()
{   
    int sIdz = (int)(g3D.shots->z[shotId] / mTomo.dz);
    int sIdx = (int)(g3D.shots->x[shotId] / mTomo.dx);
    int sIdy = (int)(g3D.shots->y[shotId] / mTomo.dy);

    int sId = sIdz + sIdx*mTomo.nz + sIdy*mTomo.nx*mTomo.nz;     

    int maxRayLength = 100000;
    float rayStep = 0.2f * m3D.dx;
    
    for (int rayId = 0; rayId < g3D.nr; rayId++)
    {
        float zi = g3D.nodes->z[rayId];
        float xi = g3D.nodes->x[rayId];
        float yi = g3D.nodes->y[rayId];

        std::vector < int > rayInd;

        for (int ds = 0; ds < maxRayLength; ds++)
        {
            int i = (int)(zi / m3D.dz) + m3D.nb;
            int j = (int)(xi / m3D.dx) + m3D.nb;
            int k = (int)(yi / m3D.dy) + m3D.nb;

            float dTz = (T[(i+1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz] - T[(i-1) + j*m3D.nzz + k*m3D.nxx*m3D.nzz]) / (2.0f*m3D.dz);    
            float dTx = (T[i + (j+1)*m3D.nzz + k*m3D.nxx*m3D.nzz] - T[i + (j-1)*m3D.nzz + k*m3D.nxx*m3D.nzz]) / (2.0f*m3D.dx);    
            float dTy = (T[i + j*m3D.nzz + (k+1)*m3D.nxx*m3D.nzz] - T[i + j*m3D.nzz + (k-1)*m3D.nxx*m3D.nzz]) / (2.0f*m3D.dy);

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            zi -= rayStep*dTz / norm; // z ray position atualization    
            xi -= rayStep*dTx / norm; // x ray position atualization   
            yi -= rayStep*dTy / norm; // y ray position atualization   

            int jm = (int)(xi / mTomo.dx); 
            int km = (int)(yi / mTomo.dy); 
            int im = (int)(zi / mTomo.dz); 

            rayInd.push_back(im + jm*mTomo.nz + km*mTomo.nx*mTomo.nz);

            if (rayInd.back() == sId) break;
        }

        float distance = rayStep;    
        float finalDist = sqrtf(powf(zi - g3D.shots->z[shotId],2.0f) + powf(xi - g3D.shots->x[shotId],2.0f) + powf(yi - g3D.shots->y[shotId],2.0f));

        std::sort(rayInd.begin(), rayInd.end());

        int current = rayInd[0];

        for (int index = 0; index < rayInd.size(); index++)
        {
            if (rayInd[index] == current)
            {
                distance += rayStep;
            }
            else
            {
                vM.push_back(distance);
                jM.push_back(current);
                iM.push_back(rayId + shotId * g3D.nr);

                if (current == sId) vM.back() = finalDist;

                distance = rayStep;
                current = rayInd[index];    
            }
        }

        if (current == sId)
        {
            vM.push_back(finalDist);
            jM.push_back(sId);
            iM.push_back(rayId + shotId * g3D.nr);
        }
        else 
        {
            vM.push_back(distance);
            jM.push_back(current);
            iM.push_back(rayId + shotId * g3D.nr);
        }

        std::vector < int >().swap(rayInd);
    }
}

void Tomography::makeGradient()
{
    for (int index = 0; index < vM.size(); index++)
        gradient[jM[index]] += vM[index] * (dobs[iM[index]] - dcal[iM[index]]);

    io.writeBinaryFloat(gradPath + "gradient_iteration_" + std::to_string(iteration) + ".bin", gradient, mTomo.nPoints);    
}

bool Tomography::converged()
{
    float r = 0.0f;
    for (int i = 0; i < g3D.nr * g3D.ns; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.emplace_back(sqrt(r));
    
    if ((residuo.back() < tomoTolerance) || (iteration >= maxIteration))
    {
        std::cout<<"\nFinal residuo: "<<sqrtf(r)<<std::endl;
        return true;
    }
    else
    {
        iteration += 1;
        return false;
    }
}

void Tomography::cgls_Berriman()
{
    std::cout<<"Solving linear system...\n\n";

    // G matrix construction

    int nnz = vM.size();

    int nM = mTomo.nPoints;                 // Model dimension
    int nD = g3D.ns * g3D.nr;               // Data dimension

    int * iG = new int[nnz + nM]();
    int * jG = new int[nnz + nM]();
    float * vG = new float[nnz + nM]();

    for (int index = 0; index < nnz; index++)
    {
        iG[index] = iM[index];
        jG[index] = jM[index];
        vG[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    // Berriman regularization 

    float * L = new float[nD]();
    float * C = new float[nM]();
    float * B = new float[nD + nM]();

    for (int index = 0; index < nnz; index++)
    {
        L[iG[index]] += vG[index];
        C[jG[index]] += vG[index];
    }    

    for (int index = 0; index < nD; index++) 
        L[index] = sqrtf(1.0f / L[index]);

    for (int index = 0; index < nnz; index++) 
        vG[index] *= L[iG[index]]; 

    for (int index = 0; index < nD; index++) 
        B[index] = L[index] * (dobs[index] - dcal[index]);

    for (int index = 0; index < nM; index++) 
    {
        iG[nnz + index] = nD + index;
        jG[nnz + index] = index;
        vG[nnz + index] = lambda * C[index];
    }

    // Conjugate gradient least square

    nD += nM;
    nnz += nM;

    int maxIt = 1000;
    int cgTol = 1e-6;

    float * x = utils.sparse_cgls(iG, jG, vG, B, nD, nM, nnz, maxIt, cgTol);

    int xcut = (int) (xMask / mTomo.dx);
    int ycut = (int) (yMask / mTomo.dy);
    int zcutUp = (int) (zMaskUp / mTomo.dz);
    int zcutDown = (int) (zMaskDown / mTomo.dz);

    // Tomography slowness parameters atualization    
    for (int index = 0; index < nM; index++) 
    {
        int k = (int) (index / (mTomo.nx*mTomo.nz));               // y direction
        int j = (int) (index - k*mTomo.nx*mTomo.nz) / mTomo.nz;    // x direction
        int i = (int) (index - j*mTomo.nz - k*mTomo.nx*mTomo.nz);  // z direction        

        if ((i >= zcutUp) && (i < mTomo.nz - zcutDown) && (j >= xcut) && (j < mTomo.nx - xcut - 1) && (k >= ycut) && (k < mTomo.ny - ycut - 1))
        {
            if (fabs(x[index]) < msv) slowness[index] += x[index];
        }
    }

    delete[] iG; delete[] jG; delete[] vG; delete[] B; delete[] x;
}

void Tomography::cgls_zoTikhonov()
{


}

void Tomography::cgls_foTikhonov()
{


}

void Tomography::cgls_soTikhonov()
{

    
}

void Tomography::modelUpdate()
{    
    // Trilinear interpolation - wikipedia

    for (int index = 0; index < m3D.nPoints; index++)
    {
        int k = (int) (index / (m3D.nx*m3D.nz));         // y direction
        int j = (int) (index - k*m3D.nx*m3D.nz) / m3D.nz;    // x direction
        int i = (int) (index - j*m3D.nz - k*m3D.nx*m3D.nz);  // z direction

        float z = i*m3D.dz; 
        float x = j*m3D.dx; 
        float y = k*m3D.dy; 

        float x0 = floor(x / mTomo.dx) * mTomo.dx;
        float x1 = floor(x / mTomo.dx) * mTomo.dx + mTomo.dx;

        float y0 = floor(y / mTomo.dy) * mTomo.dy;
        float y1 = floor(y / mTomo.dy) * mTomo.dy + mTomo.dy;

        float z0 = floor(z / mTomo.dz) * mTomo.dz;
        float z1 = floor(z / mTomo.dz) * mTomo.dz + mTomo.dz;

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);

        int zTomo = (int) z / mTomo.dz;
        int xTomo = (int) x / mTomo.dx;
        int yTomo = (int) y / mTomo.dy;

        int indS = zTomo + xTomo*mTomo.nz + yTomo*mTomo.nx*mTomo.nz;

        float c000 = slowness[indS];                  
        float c100 = slowness[indS + mTomo.nz];
        float c001 = slowness[indS + 1];
        float c101 = slowness[indS + 1 + mTomo.nz];
        float c010 = slowness[indS + mTomo.nx*mTomo.nz];
        float c110 = slowness[indS + mTomo.nz + mTomo.nx*mTomo.nz];
        float c011 = slowness[indS + 1 + mTomo.nx*mTomo.nz];
        float c111 = slowness[indS + 1 + mTomo.nz + mTomo.nx*mTomo.nz];  

        float c00 = c000*(1 - xd) + c100*xd;
        float c01 = c001*(1 - xd) + c101*xd;
        float c10 = c010*(1 - xd) + c110*xd;
        float c11 = c011*(1 - xd) + c111*xd;

        float c0 = c00*(1 - yd) + c10*yd;
        float c1 = c01*(1 - yd) + c11*yd;

        float c = c0*(1 - zd) + c1*zd; 
        
        if (c > 0.0f) m3D.vp[(i + m3D.nb) + (j+m3D.nb)*m3D.nzz + (k+m3D.nb)*m3D.nxx*m3D.nzz] = 1.0f / c;
    }

    float * mm = new float[m3D.nPoints]();
    
    for (int index  = 0; index < m3D.nPoints; index++)
    {
        int k = (int) (index / (m3D.nx*m3D.nz));             // y direction
        int j = (int) (index - k*m3D.nx*m3D.nz) / m3D.nz;    // x direction
        int i = (int) (index - j*m3D.nz - k*m3D.nx*m3D.nz);  // z direction

        mm[i + j*m3D.nz + k*m3D.nx*m3D.nz] = m3D.vp[(i + m3D.nb) + (j + m3D.nb)*m3D.nzz + (k + m3D.nb)*m3D.nxx*m3D.nzz];
    }

    io.writeBinaryFloat(estModels + "estimatedModel_iteration_" + std::to_string(iteration) + ".bin", mm, m3D.nPoints);    

    delete[] mm;
}

void Tomography::modelSmoothing()
{
    int init = m3D.nb;
    int samples = 2*m3D.nb + 1;

    // Slowness suavization - moving average filter

    if (smoothing)
    {
        float * smoothSlowness = new float[m3D.nPointsB];

        for (int i = 0; i < m3D.nPointsB; i++) smoothSlowness[i] = 1.0f / m3D.vp[i];

        for (int i = init; i < m3D.nzz - init; i++)
        {
            for (int j = init; j < m3D.nxx - init; j++)
            {
                for (int k = init; k < m3D.nyy - init; k++)
                {
                    float xs = 0.0f; 
                    float ys = 0.0f;
                    float zs = 0.0f;
                    
                    for (int s = 0; s < samples; s++)
                    {
                        int p = s - init;

                        xs += smoothSlowness[i + (j + p)*m3D.nzz + k*m3D.nxx*m3D.nzz];
                        ys += smoothSlowness[(i + p) + j*m3D.nzz + k*m3D.nxx*m3D.nzz];
                        zs += smoothSlowness[i + j*m3D.nzz + (k + p)*m3D.nxx*m3D.nzz];
                    }        

                    xs *= 1.0f / samples;
                    ys *= 1.0f / samples;
                    zs *= 1.0f / samples;

                    smoothSlowness[i + j*m3D.nzz + k*m3D.nxx*m3D.nzz] = (xs + ys + zs) / 3.0f;
                }
            }   
        }

        for (int i = 0; i < m3D.nPointsB; i++) m3D.vp[i] = 1.0f / smoothSlowness[i];

        delete[] smoothSlowness;
    }
}

void Tomography::exportConvergency()
{
    std::ofstream resFile(resPath + "residuo.txt");
    
    for (int r = 0; r < residuo.size(); r++)
        resFile<<residuo[r]<<"\n";

    resFile.close();  
}




