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
    generate_dobs = str2bool(catchParameter("generateDobs", argv[1]));
    exportTimesVolume = str2bool(catchParameter("exportTravelTimes", argv[1]));
    optimizationMethod = std::stoi(catchParameter("optimizationMethod", argv[1]));
    exportFirstArrivals = str2bool(catchParameter("exportFirstArrivals", argv[1]));

    nb = 2;

    nx = std::stoi(catchParameter("nx", argv[1]));
    ny = std::stoi(catchParameter("ny", argv[1]));
    nz = std::stoi(catchParameter("nz", argv[1]));
    
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
    { 
        std::cout<<"------- Checking final residuo ------------\n\n";
    }
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" ------------\n\n";

    std::cout<<"Shot "<<shotId+1<<" of "<<shots.all<<" at position: (x = "<<shots.x[shotId]<<" m, y = "<<shots.y[shotId]<<" m, z = "<<shots.z[shotId]<<" m)\n\n";

    if (iteration > 0) std::cout<<"Previous iteration residuous: "<<residuo.back()<<"\n\n";
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

void Tomography::tomographyUpdate()
{
    int xcut = (int) (xMask / mTomo.dx);
    int ycut = (int) (yMask / mTomo.dy);
    int zcutUp = (int) (zMaskUp / mTomo.dz);
    int zcutDown = (int) (zMaskDown / mTomo.dz);

    // Tomography slowness update    
    for (int index = 0; index < mTomo.nPoints; index++) 
    {
        int k = (int) (index / (mTomo.nx*mTomo.nz));               // y direction
        int j = (int) (index - k*mTomo.nx*mTomo.nz) / mTomo.nz;    // x direction
        int i = (int) (index - j*mTomo.nz - k*mTomo.nx*mTomo.nz);  // z direction        

        if ((i >= zcutUp) && (i < mTomo.nz - zcutDown) && (j >= xcut) && (j < mTomo.nx - xcut - 1) && (k >= ycut) && (k < mTomo.ny - ycut - 1))
        {
            slowness[index] += deltaSlowness[index];
        }
    }    
}

void Tomography::lscg_Berriman()
{
    std::cout<<"Solving linear system with berriman regularization...\n\n";

    // G matrix construction

    sparseMatrix G;

    G.nnz = vM.size();           // Non zero elements 
    G.m = mTomo.nPoints;         // Model dimension
    G.n = shots.all*nodes.all;   // Data dimension

    G.i = new int[G.nnz];
    G.j = new int[G.nnz];
    G.v = new float[G.nnz];

    for (int index = 0; index < vM.size(); index++)
    {
        G.i[index] = iM[index];
        G.j[index] = jM[index];
        G.v[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    float * L = new float[G.n]();
    float * C = new float[G.m]();

    // Berriman regularization 

    for (int index = 0; index < G.nnz; index++)
    {
        L[G.i[index]] += G.v[index];
        C[G.j[index]] += G.v[index];
    }    

    for (int index = 0; index < G.n; index++) 
        L[index] = sqrtf(1.0f / L[index]);

    for (int index = 0; index < G.nnz; index++) 
        G.v[index] *= L[G.i[index]]; 

    // Constructing a full matrix to solve problem with sparse least square conjugate gradient

    sparseMatrix A;

    A.nnz = G.nnz + mTomo.nPoints;
    A.n = G.n + mTomo.nPoints;
    A.m = G.m;  

    A.i = new int[A.nnz];
    A.j = new int[A.nnz];
    A.v = new float[A.nnz];

    for (int index = 0; index < A.nnz; index++) 
    {
        if (index < G.nnz)
        {
            A.i[index] = G.i[index];
            A.j[index] = G.j[index];
            A.v[index] = G.v[index];
        }
        else
        {
            A.i[index] = G.n + (index - G.nnz);
            A.j[index] = index - G.nnz;
            A.v[index] = lambda * C[index - G.nnz];        
        }
    }

    // Filling the data difference array with regularization

    float * B = new float[A.n]();

    for (int index = 0; index < G.n; index++) 
        B[index] = L[index] * (dobs[index] - dcal[index]);

    delete[] G.i; delete[] G.j; delete[] G.v; delete[] C; delete L;

    // Sparse least squares conjugate gradient

    int maxIt = 1000;
    int cgTol = 1e-6;

    deltaSlowness = sparse_lscg(A, B, maxIt, cgTol);

    tomographyUpdate();

    delete[] A.i; delete[] A.j; delete[] A.v; delete[] B; delete[] deltaSlowness;
}

void Tomography::lscg_zoTikhonov()
{
    std::cout<<"Solving linear system with zero order Tikhonov regularization...\n\n";

    // G matrix construction

    sparseMatrix G;

    G.nnz = vM.size();           // Non zero elements 
    G.m = mTomo.nPoints;         // Model dimension
    G.n = shots.all*nodes.all;   // Data dimension

    G.i = new int[G.nnz];
    G.j = new int[G.nnz];
    G.v = new float[G.nnz];

    for (int index = 0; index < vM.size(); index++)
    {
        G.i[index] = iM[index];
        G.j[index] = jM[index];
        G.v[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    // Constructing a full matrix to solve problem with sparse least square conjugate gradient

    sparseMatrix A;

    A.nnz = G.nnz + mTomo.nPoints;
    A.n = G.n + mTomo.nPoints;
    A.m = G.m;  

    A.i = new int[A.nnz];
    A.j = new int[A.nnz];
    A.v = new float[A.nnz];

    for (int index = 0; index < A.nnz; index++) 
    {
        if (index < G.nnz)
        {
            A.i[index] = G.i[index];
            A.j[index] = G.j[index];
            A.v[index] = G.v[index];
        }
        else
        {
            A.i[index] = G.n + (index - G.nnz);
            A.j[index] = index - G.nnz;
            A.v[index] = lambda;        
        }
    }

    // Filling the data difference array

    float * B = new float[A.n]();

    for (int index = 0; index < G.n; index++) 
        B[index] = dobs[index] - dcal[index];

    delete[] G.i; delete[] G.j; delete[] G.v;

    // Sparse least squares conjugate gradient

    int maxIt = 1000;
    int cgTol = 1e-6;

    deltaSlowness = sparse_lscg(A, B, maxIt, cgTol);

    tomographyUpdate();

    delete[] A.i; delete[] A.j; delete[] A.v; delete[] B; delete[] deltaSlowness;
}

void Tomography::lscg_foTikhonov()
{
    std::cout<<"Solving linear system with first order Tikhonov regularization...\n\n";

    // G matrix construction

    sparseMatrix G;

    G.nnz = vM.size();           // Non zero elements 
    G.m = mTomo.nPoints;         // Model dimension
    G.n = shots.all*nodes.all;   // Data dimension

    G.i = new int[G.nnz]();
    G.j = new int[G.nnz]();
    G.v = new float[G.nnz]();

    for (int index = 0; index < vM.size(); index++)
    {
        G.i[index] = iM[index];
        G.j[index] = jM[index];
        G.v[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    // First order tikhonov regularization 

    auto L = getDerivativeMatrix(G.m, 1);

    // Constructing a full matrix to solve problem with sparse least square conjugate gradient

    sparseMatrix A;

    A.nnz = G.nnz + L.nnz;
    A.n = G.n + L.n;
    A.m = G.m;  

    A.i = new int[A.nnz]();
    A.j = new int[A.nnz]();
    A.v = new float[A.nnz]();

    for (int index = 0; index < A.nnz; index++) 
    {
        if (index < G.nnz)
        {
            A.i[index] = G.i[index];
            A.j[index] = G.j[index];
            A.v[index] = G.v[index];
        }
        else
        {
            A.i[index] = G.n + L.i[index - G.nnz];
            A.j[index] = L.j[index - G.nnz];
            A.v[index] = lambda * L.v[index - G.nnz];        
        }
    }

    // Filling the data difference array

    float * B = new float[A.n]();

    for (int index = 0; index < G.n; index++) 
        B[index] = dobs[index] - dcal[index];

    delete[] G.i; delete[] G.j; delete[] G.v;
    delete[] L.i; delete[] L.j; delete[] L.v;

    // Sparse least squares conjugate gradient

    int maxIt = 1000;
    int cgTol = 1e-6;

    deltaSlowness = sparse_lscg(A, B, maxIt, cgTol);

    tomographyUpdate();

    delete[] A.i; delete[] A.j; delete[] A.v; delete[] B; delete[] deltaSlowness;
}

void Tomography::lscg_soTikhonov()
{
    std::cout<<"Solving linear system with second order Tikhonov regularization...\n\n";

    // G matrix construction

    sparseMatrix G;

    G.nnz = vM.size();           // Non zero elements 
    G.m = mTomo.nPoints;         // Model dimension
    G.n = shots.all*nodes.all;   // Data dimension

    G.i = new int[G.nnz];
    G.j = new int[G.nnz];
    G.v = new float[G.nnz];

    for (int index = 0; index < vM.size(); index++)
    {
        G.i[index] = iM[index];
        G.j[index] = jM[index];
        G.v[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    // Second order tikhonov regularization 

    auto L = getDerivativeMatrix(G.m, 2);

    // Constructing a full matrix to solve problem with sparse least square conjugate gradient

    sparseMatrix A;

    A.nnz = G.nnz + L.nnz;
    A.n = G.n + L.n;
    A.m = G.m;  

    A.i = new int[A.nnz];
    A.j = new int[A.nnz];
    A.v = new float[A.nnz];

    for (int index = 0; index < A.nnz; index++) 
    {
        if (index < G.nnz)
        {
            A.i[index] = G.i[index];
            A.j[index] = G.j[index];
            A.v[index] = G.v[index];
        }
        else
        {
            A.i[index] = G.n + L.i[index - G.nnz];
            A.j[index] = L.j[index - G.nnz];
            A.v[index] = lambda * L.v[index - G.nnz];        
        }
    }

    // Filling the data difference array

    float * B = new float[A.n]();

    for (int index = 0; index < G.n; index++) 
        B[index] = dobs[index] - dcal[index];

    delete[] G.i; delete[] G.j; delete[] G.v;
    delete[] L.i; delete[] L.j; delete[] L.v;

    // Sparse least squares conjugate gradient

    int maxIt = 1000;
    int cgTol = 1e-6;

    deltaSlowness = sparse_lscg(A, B, maxIt, cgTol);

    tomographyUpdate();

    delete[] A.i; delete[] A.j; delete[] A.v; delete[] B; delete[] deltaSlowness;    
}

void Tomography::gradientDescent()
{
    gradient = new float[mTomo.nPoints];

    for (int index = 0; index < vM.size(); index++)
    {
        gradient[jM[index]] += vM[index] * (dobs[iM[index]] - dcal[iM[index]]);
    }

    writeBinaryFloat("gradient.bin",gradient,mTomo.nPoints);

    int samples = 11; 

    gradient = movingAverageSmoothing(gradient,mTomo.nx,mTomo.ny,mTomo.nz,samples);

    writeBinaryFloat("gradient_smoothed.bin", gradient, mTomo.nPoints);
}

void Tomography::optimization()
{
    switch (optimizationMethod)
    {
    case 0:
        lscg_Berriman();    
        break;
    
    case 1:
        lscg_zoTikhonov();
        break;
    
    case 2:
        lscg_foTikhonov();
        break;
    
    case 3:
        lscg_soTikhonov();
        break;

    case 4:
        gradientDescent();    
    }
}

void Tomography::modelUpdate()
{    
    // Trilinear interpolation - wikipedia

    int xcut = (int) (xMask / dx);
    int ycut = (int) (yMask / dy);
    int zcutUp = (int) (zMaskUp / dz);
    int zcutDown = (int) (zMaskDown / dz);

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

        if ((i >= zcutUp) && (i < nz - zcutDown - 1) && (j >= xcut) && (j < nx - xcut - 1) && (k >= ycut) && (k < ny - ycut - 1))
        {
            int idz = ((int)(z/mTomo.dz));
            int idx = ((int)(x/mTomo.dx));
            int idy = ((int)(y/mTomo.dy));

            int indS = idz + idx*mTomo.nz + idy*mTomo.nx*mTomo.nz;

            float c000 = slowness[indS];                  
            float c001 = slowness[indS + 1];
            float c100 = slowness[indS + mTomo.nz];
            float c101 = slowness[indS + 1 + mTomo.nz];
            float c010 = slowness[indS + mTomo.nx*mTomo.nz];
            float c011 = slowness[indS + 1 + mTomo.nx*mTomo.nz];
            float c110 = slowness[indS + mTomo.nz + mTomo.nx*mTomo.nz];
            float c111 = slowness[indS + 1 + mTomo.nz + mTomo.nx*mTomo.nz];  

            float s = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);
    
            Vp[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz] = 1.0f / s;
        }
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

void Tomography::exportConvergency()
{
    std::ofstream resFile(resPath + "residuo.txt");
    
    for (int r = 0; r < residuo.size(); r++)
        resFile<<residuo[r]<<"\n";

    resFile.close();
}

