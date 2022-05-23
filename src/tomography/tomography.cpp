# include <cmath>
# include <string>
# include <algorithm>

# include "tomography.hpp"

# include "../simulations/eikonal/eikonal.hpp"

Tomography3D::Tomography3D(char **argv)
{
    parametersFile = argv[1];

    setup();

    dobsPath = io.catchParameter("dobsFolder", parametersFile);
    dcalPath = io.catchParameter("dcalFolder", parametersFile);
    gradPath = io.catchParameter("gradientFolder", parametersFile);

    maxIteration = std::stoi(io.catchParameter("maxIteration", parametersFile));
    tomoTolerance = std::stof(io.catchParameter("tomoTolerance", parametersFile));

    generate_dobs = utils.str2bool(io.catchParameter("generate_dobs", parametersFile));

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
                podvin3D();
            else 
                fim3D();
        }
    }

    m3D.vpPath = io.catchParameter("initModelPath", parametersFile);

    m3D.readAndExpandVP();

    iteration = 1;

    dobs = new float[g3D.ns * g3D.nr]();
    dcal = new float[g3D.ns * g3D.nr]();

    gradient = new float[m3D.nPoints]();
}

void Tomography3D::importDobs()
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

void Tomography3D::importDcal()
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

void Tomography3D::fwdModeling()
{
    arrivalsPath = dcalPath;

    for (int shot = 0; shot < g3D.ns; shot++)
    {
        std::cout<<"Computing shot "<<shot+1<<" of "<<g3D.ns<<"\n";

        shotId = shot;

        if (eikonalType == 1) 
            podvin3D();
        else 
            fim3D();

        gradientRayTracing();
    }
}

void Tomography3D::gradientRayTracing()
{   
    int sIdz = (int)(g3D.shots->z[shotId] / m3D.dz);
    int sIdx = (int)(g3D.shots->x[shotId] / m3D.dx);
    int sIdy = (int)(g3D.shots->y[shotId] / m3D.dy);

    int sId = sIdz + sIdx*m3D.nz + sIdy*m3D.nx*m3D.nz;     

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

            int jm = (int)(xi / m3D.dx); 
            int km = (int)(yi / m3D.dy); 
            int im = (int)(zi / m3D.dz); 

            rayInd.push_back(im + jm*m3D.nz + km*m3D.nx*m3D.nz);

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

void Tomography3D::makeGradient()
{
    for (int index = 0; index < vM.size(); index++)
        gradient[jM[index]] += vM[index] * (dobs[iM[index]] - dcal[iM[index]]);

    io.writeBinaryFloat(gradPath + "gradient_iteration_" + std::to_string(iteration) + ".bin", gradient, m3D.nPoints);    
}

bool Tomography3D::converged()
{
    float r = 0.0f;
    for (int i = 0; i < g3D.nr * g3D.ns; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    if ((r < tomoTolerance) || (iteration >= maxIteration))
    {
        std::cout<<"\nFinal residuo: "<<r<<std::endl;
        return true;
    }
    else
    {
        return false;
    }
}

