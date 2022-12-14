# include <cmath>
# include <limits>
# include <string>
# include <iomanip>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "tomography.hpp"

Tomography::Tomography() { } 

void Tomography::setParameters(char * parameters)
{
    eikonalType = std::stoi(catchParameter("eikonalType", parameters));    
    exportTimesVolume = str2bool(catchParameter("exportTravelTimes", parameters));
    exportFirstArrivals = str2bool(catchParameter("exportFirstArrivals", parameters));

    nb = 2;
    nx = std::stoi(catchParameter("nx", parameters));
    ny = std::stoi(catchParameter("ny", parameters));
    nz = std::stoi(catchParameter("nz", parameters));
    
    initialize();

    dx = std::stof(catchParameter("dx", parameters));
    dy = std::stof(catchParameter("dy", parameters));
    dz = std::stof(catchParameter("dz", parameters));
    
    shotsGeometryType = std::stoi(catchParameter("shotsGeometryType", parameters));
    nodesGeometryType = std::stoi(catchParameter("nodesGeometryType", parameters));

    reciprocity = str2bool(catchParameter("reciprocity", parameters));
    saveGeometry = str2bool(catchParameter("saveGeometry", parameters));

    std::vector<std::string> splitted;

    shots.elevation = std::stof(catchParameter("shotsElevation", parameters));
    nodes.elevation = std::stof(catchParameter("nodesElevation", parameters));

    if (shotsGeometryType) // Grid shots aqcuisition
    {
        shots.n_xline = std::stoi(catchParameter("xShotNumber", parameters));
        shots.n_yline = std::stoi(catchParameter("yShotNumber", parameters));
        
        splitted = split(catchParameter("shotSW", parameters),',');
        set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("shotNW", parameters),',');
        set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("shotSE", parameters),',');
        set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

        setGridGeometry(shots);
    }
    else                   // Circular shots aqcuisition
    {   
        shots.xcenter = std::stoi(catchParameter("xShotCenter", parameters));
        shots.ycenter = std::stoi(catchParameter("yShotCenter", parameters));
        shots.circle_spacing = std::stof(catchParameter("shotSpacing", parameters));
        
        splitted = split(catchParameter("sOffsets", parameters),',');

        for (auto offset : splitted) 
            shots.offsets.push_back(std::stof(offset));

        setCircularGeometry(shots);
    }

    if (nodesGeometryType)           // Grid nodes aqcuisition
    {
        nodes.n_xline = std::stoi(catchParameter("xNodeNumber", parameters));
        nodes.n_yline = std::stoi(catchParameter("yNodeNumber", parameters));
        
        splitted = split(catchParameter("nodeSW", parameters),',');
        set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("nodeNW", parameters),',');
        set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

        splitted = split(catchParameter("nodeSE", parameters),',');
        set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

        setGridGeometry(nodes);
    }
    else                        // Circular nodes aqcuisition
    {
        nodes.xcenter = std::stoi(catchParameter("xNodeCenter", parameters));
        nodes.ycenter = std::stoi(catchParameter("yNodeCenter", parameters));
        nodes.circle_spacing = std::stof(catchParameter("nodeSpacing", parameters));

        splitted = split(catchParameter("nodeOffsets", parameters),',');

        for (auto offset : splitted) 
            nodes.offsets.push_back(std::stof(offset));

        setCircularGeometry(nodes);
    }

    shotsPath = catchParameter("shotsPositionPath", parameters);
    nodesPath = catchParameter("nodesPositionPath", parameters);

    eikonalFolder = catchParameter("eikonalVolumeFolder", parameters);

    if (saveGeometry) exportPositions();
    if (reciprocity) setReciprocity();

    dobsPath = catchParameter("dobsPath", parameters);
    dcalPath = catchParameter("dcalPath", parameters);

    lambda = std::stof(catchParameter("regParam", parameters));
    tkOrder = std::stof(catchParameter("regOrder", parameters));
    maxIteration = std::stoi(catchParameter("maxIteration", parameters));

    smooth = str2bool(catchParameter("smooth", parameters)); 
    smoothingType = std::stoi(catchParameter("smoothingType", parameters));
    filterSamples = std::stoi(catchParameter("filterSamples", parameters));
    standardDeviation = std::stof(catchParameter("standardDeviation",parameters));

    mTomo.nx = std::stoi(catchParameter("nxTomo", parameters));
    mTomo.ny = std::stoi(catchParameter("nyTomo", parameters));
    mTomo.nz = std::stoi(catchParameter("nzTomo", parameters));

    mTomo.dx = std::stof(catchParameter("dxTomo", parameters));
    mTomo.dy = std::stof(catchParameter("dyTomo", parameters));
    mTomo.dz = std::stof(catchParameter("dzTomo", parameters));

    mTomo.nPoints = mTomo.nx * mTomo.ny * mTomo.nz;

    residuoPath = catchParameter("convergencyFolder", parameters);
    estimatedPath = catchParameter("estimatedModelsFolder", parameters);

    xMask = split(catchParameter("xMask", parameters),',');
    yMask = split(catchParameter("yMask", parameters),',');
    zMask = split(catchParameter("zMask", parameters),',');

    V = expand(readBinaryFloat(catchParameter("modelPath", parameters), nPoints));

    iteration = 0;

    residuo.reserve(maxIteration);

    dobs = new float[shots.all * nodes.all]();
    dcal = new float[shots.all * nodes.all]();    
    
    model = new float [mTomo.nPoints];
}

void Tomography::infoMessage()
{
    std::cout << std::fixed << std::setprecision(1);

    system("clear");
    std::cout<<"3D first arrival tomography\n\n";
    
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
                
                model[i + j*mTomo.nz + k*mTomo.nx*mTomo.nz] = 1.0f / V[im + jm*nzz + km*nxx*nzz];    
            }
        }
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

    float rayStep = 0.2f * (dx + dy + dz) / 3.0f;
    
    for (int rayId = 0; rayId < nodes.all; rayId++)
    {
        float zi = nodes.z[rayId];
        float xi = nodes.x[rayId];
        float yi = nodes.y[rayId];

        std::vector < int > rayInd;

        while (true)
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

            int im = (int)(zi / mTomo.dz); 
            int jm = (int)(xi / mTomo.dx); 
            int km = (int)(yi / mTomo.dy); 

            rayInd.push_back(im + jm*mTomo.nz + km*mTomo.nx*mTomo.nz);

            if (rayInd.back() == sId) break;
        }
    
        float finalDist = sqrtf(powf(zi - shots.z[shotId],2.0f) + powf(xi - shots.x[shotId],2.0f) + powf(yi - shots.y[shotId],2.0f));

        std::sort(rayInd.begin(), rayInd.end());

        int current = rayInd[0];
        float distance = rayStep;

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

bool Tomography::converged()
{
    float r = 0.0f;
    for (int i = 0; i < nodes.all * shots.all; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.emplace_back(sqrtf(r));
    
    if (iteration >= maxIteration)
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

Utils::sparseMatrix Tomography::buildForwardModelingMatrix()
{
    sparseMatrix G;

    G.nnz = vM.size();           // Non zero elements 
    G.m = mTomo.nPoints;         // Model dimension
    G.n = shots.all*nodes.all;   // Data dimension

    G.init();

    for (int index = 0; index < vM.size(); index++)
    {
        G.i[index] = iM[index];
        G.j[index] = jM[index];
        G.v[index] = vM[index];
    }

    std::vector<  int  >().swap(iM);
    std::vector<  int  >().swap(jM);
    std::vector< float >().swap(vM);

    return G;
}

Utils::sparseMatrix Tomography::applyTkRegularization(sparseMatrix G, sparseMatrix L)
{
    sparseMatrix A;

    A.nnz = G.nnz + L.nnz;
    A.n = G.n + L.n;
    A.m = G.m;  

    A.init();

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

    return A;
}

void Tomography::optimization()
{
    std::cout<<"Solving linear system using Tikhonov regularization with order "+ std::to_string(tkOrder) + "\n\n";

    sparseMatrix G = buildForwardModelingMatrix();

    sparseMatrix L = getDerivativeMatrix(G.m, tkOrder);

    sparseMatrix A = applyTkRegularization(G, L);

    G.erase(); L.erase();

    // Filling the data difference array

    float * B = new float[A.n]();

    for (int index = 0; index < nodes.all*shots.all; index++) 
        B[index] = dobs[index] - dcal[index];

    // Sparse least squares conjugate gradient

    int maxIt = 1000;
    int cgTol = 1e-6;

    dm = sparse_lscg(A, B, maxIt, cgTol);

    A.erase(); delete[] B; 
}

void Tomography::modelUpdate()
{    
    int xlCut = (int)(std::stof(xMask[0]) / mTomo.dx); // x left cut 
    int xrCut = (int)(std::stof(xMask[1]) / mTomo.dx); // x right cut
    int ylCut = (int)(std::stof(yMask[0]) / mTomo.dy); // y left cut 
    int yrCut = (int)(std::stof(yMask[1]) / mTomo.dy); // y right cut
    int zuCut = (int)(std::stof(zMask[0]) / mTomo.dz); // z up cut  
    int zdCut = (int)(std::stof(zMask[1]) / mTomo.dz); // z down cut

    // Tomography slowness update    
    for (int index = 0; index < mTomo.nPoints; index++) 
    {
        int k = (int) (index / (mTomo.nx*mTomo.nz));               // y direction
        int j = (int) (index - k*mTomo.nx*mTomo.nz) / mTomo.nz;    // x direction
        int i = (int) (index - j*mTomo.nz - k*mTomo.nx*mTomo.nz);  // z direction        

        if ((i >= zuCut) && (i < mTomo.nz - zdCut) && (j >= xlCut) && (j < mTomo.nx - xrCut - 1) && (k >= ylCut) && (k < mTomo.ny - yrCut - 1))
        { 
            model[index] += dm[index];
        }
    }    

    // Trilinear interpolation - wikipedia
    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         // y direction
        int j = (int) (index - k*nx*nz) / nz;    // x direction
        int i = (int) (index - j*nz - k*nx*nz);  // z direction

        float x = j*dx; 
        float y = k*dy; 
        float z = i*dz; 

        float x0 = floorf(x/mTomo.dx)*mTomo.dx;
        float y0 = floorf(y/mTomo.dy)*mTomo.dy;
        float z0 = floorf(z/mTomo.dz)*mTomo.dz;

        float x1 = floorf(x/mTomo.dx)*mTomo.dx + mTomo.dx;
        float y1 = floorf(y/mTomo.dy)*mTomo.dy + mTomo.dy;
        float z1 = floorf(z/mTomo.dz)*mTomo.dz + mTomo.dz;

        if ((i >= 0) && (i < nz - 1) && (j >= 0) && (j < nx - 1) && (k >= 0) && (k < ny - 1))
        {
            int idz = ((int)(z/mTomo.dz));
            int idx = ((int)(x/mTomo.dx));
            int idy = ((int)(y/mTomo.dy));

            int indS = idz + idx*mTomo.nz + idy*mTomo.nx*mTomo.nz;

            float c000 = model[indS];                  
            float c001 = model[indS + 1];
            float c100 = model[indS + mTomo.nz];
            float c101 = model[indS + 1 + mTomo.nz];
            float c010 = model[indS + mTomo.nx*mTomo.nz];
            float c011 = model[indS + 1 + mTomo.nx*mTomo.nz];
            float c110 = model[indS + mTomo.nz + mTomo.nx*mTomo.nz];
            float c111 = model[indS + 1 + mTomo.nz + mTomo.nx*mTomo.nz];  

            float s = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);

            V[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz] = 1.0f / s;
        }
    }

    if (smooth) 
    {
        float * auxS = new float[nPointsB];

        // Smoothing the inverse of velocity
        for (int i = 0; i < nPointsB; i++) auxS[i] = 1.0f / V[i];

        if (smoothingType) 
            movingAverageSmoothing(auxS, nxx, nyy, nzz, filterSamples);
        else
            gaussianFilterSmoothing(auxS, nxx, nyy, nzz, standardDeviation, filterSamples);

        // Coming back to velocity    
        for (int i = 0; i < nPointsB; i++) V[i] = 1.0f / auxS[i];
    
        delete[] auxS;
    }

    float * mm = reduce(V);
    
    writeBinaryFloat(estimatedPath + "estimatedModel_iteration_" + std::to_string(iteration) + ".bin", mm, nPoints);    

    delete[] mm;
}

void Tomography::exportConvergency()
{
    std::ofstream resFile(residuoPath + "convergency.txt", std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) resFile<<residuo[r]<<"\n";

    resFile.close();
}

