# include <cmath>
# include <limits>
# include <string>
# include <iomanip>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "tomography.hpp"

Tomography::Tomography(char **argv)
{
    eikonalType = std::stoi(catchParameter("eikonalType", argv[1]));    
    exportTimesVolume = str2bool(catchParameter("exportTravelTimes", argv[1]));
    exportFirstArrivals = str2bool(catchParameter("exportFirstArrivals", argv[1]));

    nx = std::stoi(catchParameter("nx", argv[1]));
    ny = std::stoi(catchParameter("ny", argv[1]));
    nz = std::stoi(catchParameter("nz", argv[1]));
    nb = std::stoi(catchParameter("nb", argv[1]));
    
    dx = std::stof(catchParameter("dx", argv[1]));
    dy = std::stof(catchParameter("dy", argv[1]));
    dz = std::stof(catchParameter("dz", argv[1]));
    
    shotsGeometryType = std::stoi(catchParameter("shotsGeometryType", argv[1]));
    nodesGeometryType = std::stoi(catchParameter("nodesGeometryType", argv[1]));

    reciprocity = str2bool(catchParameter("reciprocity", argv[1]));
    saveGeometry = str2bool(catchParameter("saveGeometry", argv[1]));

    std::vector<std::string> splitted;

    shots.elevation = std::stof(catchParameter("sElev", argv[1]));
    nodes.elevation = std::stof(catchParameter("rElev", argv[1]));

    if (shotsGeometryType)              // Grid shots aqcuisition
    {
        shots.n_xline = std::stoi(catchParameter("nsx", argv[1]));
        shots.n_yline = std::stoi(catchParameter("nsy", argv[1]));
        
        splitted = split(catchParameter("SWs", argv[1]),',');
        set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("NWs", argv[1]),',');
        set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("SEs", argv[1]),',');
        set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

        setGridShots();
    }
    else                             // Circular shots aqcuisition
    {   
        shots.xcenter = std::stoi(catchParameter("sxc", argv[1]));
        shots.ycenter = std::stoi(catchParameter("syc", argv[1]));
        shots.circle_spacing = std::stof(catchParameter("sds", argv[1]));
        
        splitted = split(catchParameter("sOffsets", argv[1]),',');

        for (auto offset : splitted) 
            shots.offsets.push_back(std::stof(offset));

        setCircularShots();
    }

    if (nodesGeometryType)           // Grid nodes aqcuisition
    {
        nodes.n_xline = std::stoi(catchParameter("nrx", argv[1]));
        nodes.n_yline = std::stoi(catchParameter("nry", argv[1]));
        
        splitted = split(catchParameter("SWr", argv[1]),',');
        set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("NWr", argv[1]),',');
        set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("SEr", argv[1]),',');
        set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

        setGridNodes();
    }
    else                        // Circular nodes aqcuisition
    {
        nodes.xcenter = std::stoi(catchParameter("rxc", argv[1]));
        nodes.ycenter = std::stoi(catchParameter("ryc", argv[1]));
        nodes.circle_spacing = std::stof(catchParameter("rds", argv[1]));

        splitted = split(catchParameter("rOffsets", argv[1]),',');

        for (auto offset : splitted) 
            nodes.offsets.push_back(std::stof(offset));

        setCircularNodes();
    }

    shotsPath = catchParameter("shotsPositionPath", argv[1]);
    nodesPath = catchParameter("nodesPositionPath", argv[1]);

    eikonalFolder = catchParameter("travelTimesFolder", argv[1]);
    arrivalFolder = catchParameter("firstArrivalsFolder", argv[1]);

    if (saveGeometry) exportPositions();

    if (reciprocity) setReciprocity();

    dobsPath = catchParameter("dobsFolder", argv[1]);
    dcalPath = catchParameter("dcalFolder", argv[1]);
    gradPath = catchParameter("gradientFolder", argv[1]);

    maxIteration = std::stoi(catchParameter("maxIteration", argv[1]));
    tomoTolerance = std::stof(catchParameter("tomoTolerance", argv[1]));
    lambda = std::stof(catchParameter("regParam", argv[1]));

    mTomo.nx = std::stoi(catchParameter("nxTomo", argv[1]));
    mTomo.ny = std::stoi(catchParameter("nyTomo", argv[1]));
    mTomo.nz = std::stoi(catchParameter("nzTomo", argv[1]));

    mTomo.dx = std::stof(catchParameter("dxTomo", argv[1]));
    mTomo.dy = std::stof(catchParameter("dyTomo", argv[1]));
    mTomo.dz = std::stof(catchParameter("dzTomo", argv[1]));

    mTomo.nPoints = mTomo.nx * mTomo.ny * mTomo.nz;

    resPath = catchParameter("convergencyFolder", argv[1]);
    estModels = catchParameter("estimatedModelsFolder", argv[1]);

    xMask = std::stof(catchParameter("xMask", argv[1]));
    yMask = std::stof(catchParameter("yMask", argv[1]));
    zMaskUp = std::stof(catchParameter("zMaskUp", argv[1]));
    zMaskDown = std::stof(catchParameter("zMaskDown", argv[1]));

    smoothing = str2bool(catchParameter("smoothing", argv[1]));

    generate_dobs = str2bool(catchParameter("generateDobs", argv[1]));

    initialize();
    
    if (generate_dobs)
    {
        arrivalFolder = dobsPath;

        model = readBinaryFloat(catchParameter("trueModelPath", argv[1]), nPoints);
        Vp = expandModel();

        std::cout<<"Computing observed data"<<std::endl;

        T = new float[nPointsB];
        
        for (shotId = 0; shotId < shots.all; shotId++)
        {
            eikonalComputing();    
        }

        delete[] T;
    }

    model = readBinaryFloat(catchParameter("initModelPath", argv[1]), nPoints);
    Vp = expandModel();

    iteration = 0;

    residuo.reserve(maxIteration);

    dobs = new float[shots.all * nodes.all]();
    dcal = new float[shots.all * nodes.all]();

    gradient = new float[mTomo.nPoints]();
    slowness = new float[mTomo.nPoints]();
}

void Tomography::infoMessage()
{
    std::cout << std::fixed << std::setprecision(1);

    system("clear");
    std::cout<<"3D first arrival tomography for OBN geometry\n\n";
    
    std::cout<<"Total x model length = "<<(nx-1)*dx<<" m\n";
    std::cout<<"Total Y model length = "<<(ny-1)*dy<<" m\n";
    std::cout<<"Total Z model length = "<<(nz-1)*dz<<" m\n\n";
    
    if (reciprocity)
        std::cout<<"Reciprocity = True\n\n";
    else
        std::cout<<"Shots reciprocity = False\n\n";

    if (iteration == maxIteration)
        std::cout<<"------- Checking final residuo ------------\n\n";
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" ------------\n\n";

    std::cout<<"Shot "<<shotId+1<<" of "<<shots.all<<" at position: (x = "<<shots.x[shotId]<<" m, y = "<<shots.y[shotId]<<" m, z = "<<shots.z[shotId]<<" m)\n\n";

    if (iteration > 1) std::cout<<"\nPrevious iteration residuous: "<<residuo.back()<<"\n";
} 

void Tomography::importDobs()
{       
    int ptr = 0; 
    
    for (int shot = 0; shot < shots.all; shot++)
    {
        float * data = readBinaryFloat(dobsPath + "times_nr" + std::to_string(nodes.all) + "_shot_" + std::to_string(shot+1) + ".bin", nodes.all);

        for (int d = ptr; d < ptr + nodes.all; d++) dobs[d] = data[d - ptr];

        ptr += nodes.all;
        
        delete[] data;
    }
}

void Tomography::setInitialModel()
{
    for (int k = 0; k < mTomo.ny; k++)
    {
        for (int j = 0; j < mTomo.nx; j++)
        {
            for (int i = 0; i < mTomo.nz; i++)
            {
                int im = (int) ((float)(i)*mTomo.dz/dz) + nb;
                int jm = (int) ((float)(j)*mTomo.dx/dx) + nb;
                int km = (int) ((float)(k)*mTomo.dy/dy) + nb;
                
                slowness[i + j*mTomo.nz + k*mTomo.nx*mTomo.nz] = 1.0f / Vp[im + jm*nzz + km*nxx*nzz];    
            }
        }
    }
}

void Tomography::importDcal()
{
    int ptr = 0; 
 
    for (int shot = 0; shot < shots.all; shot++)
    {
        float * data = readBinaryFloat(dcalPath + "times_nr" + std::to_string(nodes.all) + "_shot_" + std::to_string(shot+1) + ".bin", nodes.all);

        for (int d = ptr; d < ptr + nodes.all; d++) dcal[d] = data[d - ptr];

        ptr += nodes.all;
     
        delete[] data;
    }
}

void Tomography::forwardModeling()
{
    arrivalFolder = dcalPath;

    T = new float[nPointsB];

    for (shotId = 0; shotId < shots.all; shotId++)
    {
        eikonalComputing();

        gradientRayTracing();
    
        infoMessage();    
    }

    delete[] T;
}

void Tomography::gradientRayTracing()
{   
    int sIdz = (int)(shots.z[shotId] / mTomo.dz);
    int sIdx = (int)(shots.x[shotId] / mTomo.dx);
    int sIdy = (int)(shots.y[shotId] / mTomo.dy);

    int sId = sIdz + sIdx*mTomo.nz + sIdy*mTomo.nx*mTomo.nz;     

    int maxRayLength = 100000;
    float rayStep = 0.2f * (dx + dy + dz) / 3.0f;
    
    for (int rayId = 0; rayId < nodes.all; rayId++)
    {
        float zi = nodes.z[rayId];
        float xi = nodes.x[rayId];
        float yi = nodes.y[rayId];

        std::vector < int > rayInd;

        for (int ds = 0; ds < maxRayLength; ds++)
        {
            int i = (int)(zi / dz) + nb;
            int j = (int)(xi / dx) + nb;
            int k = (int)(yi / dy) + nb;

            float dTz = (T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*dz);    
            float dTx = (T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*dx);    
            float dTy = (T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*dy);

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
        float finalDist = sqrtf(powf(zi - shots.z[shotId],2.0f) + powf(xi - shots.x[shotId],2.0f) + powf(yi - shots.y[shotId],2.0f));

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
                iM.push_back(rayId + shotId * nodes.all);

                if (current == sId) vM.back() = finalDist;

                distance = rayStep;
                current = rayInd[index];    
            }
        }

        if (current == sId)
        {
            vM.push_back(finalDist);
            jM.push_back(sId);
            iM.push_back(rayId + shotId * nodes.all);
        }
        else 
        {
            vM.push_back(distance);
            jM.push_back(current);
            iM.push_back(rayId + shotId * nodes.all);
        }

        std::vector < int >().swap(rayInd);
    }
}

void Tomography::makeGradient()
{
    for (int index = 0; index < vM.size(); index++)
        gradient[jM[index]] += vM[index] * (dobs[iM[index]] - dcal[iM[index]]);

    writeBinaryFloat(gradPath + "gradient_iteration_" + std::to_string(iteration) + ".bin", gradient, mTomo.nPoints);    
}

bool Tomography::converged()
{
    float r = 0.0f;
    for (int i = 0; i < nodes.all * shots.all; i++)
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
    int nD = shots.all * nodes.all;     // Data dimension

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

    float * x = sparse_cgls(iG, jG, vG, B, nD, nM, nnz, maxIt, cgTol);

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
            slowness[index] += x[index];
        }
    }

    delete[] iG; delete[] jG; delete[] vG; delete[] B; delete[] x;
}

void Tomography::cgls_zoTikhonov()
{
    std::cout<<"Solving linear system...\n\n";

    // G matrix construction

    int nnz = vM.size();

    int nM = mTomo.nPoints;             // Model dimension
    int nD = shots.all * nodes.all;     // Data dimension

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

    for (int index = 0; index < nM; index++) 
    {
        iG[nnz + index] = nD + index;
        jG[nnz + index] = index;
        vG[nnz + index] = lambda;
    }




}

void Tomography::cgls_foTikhonov()
{
    std::cout<<"Solving linear system...\n\n";

    // G matrix construction

    int nnz = vM.size();

    int nM = mTomo.nPoints;                 // Model dimension
    int nD = shots.all * nodes.all;     // Data dimension

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


}

void Tomography::cgls_soTikhonov()
{
    std::cout<<"Solving linear system...\n\n";

    // G matrix construction

    int nnz = vM.size();

    int nM = mTomo.nPoints;                 // Model dimension
    int nD = shots.all * nodes.all;     // Data dimension

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

    
}

void Tomography::modelUpdate()
{    
    // Trilinear interpolation - wikipedia

    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         // y direction
        int j = (int) (index - k*nx*nz) / nz;    // x direction
        int i = (int) (index - j*nz - k*nx*nz);  // z direction

        float x = j*dx; 
        float y = k*dy; 
        float z = i*dz; 

        float x0 = floor(x/mTomo.dx)*mTomo.dx;
        float y0 = floor(y/mTomo.dy)*mTomo.dy;
        float z0 = floor(z/mTomo.dz)*mTomo.dz;

        float x1 = floor(x/mTomo.dx)*mTomo.dx + mTomo.dx;
        float y1 = floor(y/mTomo.dy)*mTomo.dy + mTomo.dy;
        float z1 = floor(z/mTomo.dz)*mTomo.dz + mTomo.dz;

        int indS = ((int)(z/mTomo.dz)) + ((int)(x/mTomo.dx))*mTomo.nz + ((int)(y/mTomo.dy))*mTomo.nx*mTomo.nz;

        float c000 = slowness[indS];                  
        float c001 = slowness[indS + 1];
        float c100 = slowness[indS + mTomo.nz];
        float c101 = slowness[indS + 1 + mTomo.nz];
        float c010 = slowness[indS + mTomo.nx*mTomo.nz];
        float c011 = slowness[indS + 1 + mTomo.nx*mTomo.nz];
        float c110 = slowness[indS + mTomo.nz + mTomo.nx*mTomo.nz];
        float c111 = slowness[indS + 1 + mTomo.nz + mTomo.nx*mTomo.nz];  

        float s = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);

        if (s > 0.0f) Vp[(i + nb) + (j+nb)*nzz + (k+nb)*nxx*nzz] = 1.0f / s;
    }

    float * mm = new float[nPoints]();
    
    for (int index  = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));             // y direction
        int j = (int) (index - k*nx*nz) / nz;    // x direction
        int i = (int) (index - j*nz - k*nx*nz);  // z direction

        mm[i + j*nz + k*nx*nz] = Vp[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz];
    }

    writeBinaryFloat(estModels + "estimatedModel_iteration_" + std::to_string(iteration) + ".bin", mm, nPoints);    

    delete[] mm;
}

void Tomography::modelSmoothing()
{
    int init = nb;
    int samples = 2*nb + 1;

    // Slowness suavization - moving average filter

    if (smoothing)
    {
        float * smoothSlowness = new float[nPointsB];

        for (int i = 0; i < nPointsB; i++) smoothSlowness[i] = 1.0f / Vp[i];

        for (int i = init; i < nzz - init; i++)
        {
            for (int j = init; j < nxx - init; j++)
            {
                for (int k = init; k < nyy - init; k++)
                {
                    float xs = 0.0f; 
                    float ys = 0.0f;
                    float zs = 0.0f;
                    
                    for (int s = 0; s < samples; s++)
                    {
                        int p = s - init;

                        xs += smoothSlowness[i + (j + p)*nzz + k*nxx*nzz];
                        ys += smoothSlowness[(i + p) + j*nzz + k*nxx*nzz];
                        zs += smoothSlowness[i + j*nzz + (k + p)*nxx*nzz];
                    }        

                    xs *= 1.0f / samples;
                    ys *= 1.0f / samples;
                    zs *= 1.0f / samples;

                    smoothSlowness[i + j*nzz + k*nxx*nzz] = (xs + ys + zs) / 3.0f;
                }
            }   
        }

        for (int i = 0; i < nPointsB; i++) Vp[i] = 1.0f / smoothSlowness[i];

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




