# include <cmath>
# include <iostream>

# include "acoustic.hpp"

void Acoustic::setParameters()
{
    /* Importing time domain paramFile */

    nt = std::stoi(catchParameter("nt", paramFile));    
    dt = std::stof(catchParameter("dt", paramFile));
    tlag = std::stof(catchParameter("tlag", paramFile));
    fmax = std::stof(catchParameter("fmax", paramFile));

    seisLabel = catchParameter("seisLabel", paramFile);

    /* Importing model paramFile */

    nx = std::stoi(catchParameter("nx", paramFile));
    ny = std::stoi(catchParameter("ny", paramFile));
    nz = std::stoi(catchParameter("nz", paramFile));

    dx = std::stof(catchParameter("dx", paramFile));
    dy = std::stof(catchParameter("dy", paramFile));
    dz = std::stof(catchParameter("dz", paramFile));

    nb = std::stoi(catchParameter("nb", paramFile));
    
    factor = std::stof(catchParameter("factor", paramFile));

    initialize();
    sourceGenerator();
    dampingGenerator();

    /* Importing geometry paramFile */

    reciprocity = str2bool(catchParameter("reciprocity", paramFile));
    saveGeometry = str2bool(catchParameter("saveGeometry", paramFile));

    shots.elevation = std::stof(catchParameter("shotsElevation", paramFile));
    nodes.elevation = std::stof(catchParameter("nodesElevation", paramFile));

    std::vector<std::string> splitted;

    /* Setting nodes position */

    shots.n_xline = std::stoi(catchParameter("xShotNumber", paramFile));
    shots.n_yline = std::stoi(catchParameter("yShotNumber", paramFile));
    
    splitted = split(catchParameter("shotSW", paramFile),',');
    set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("shotNW", paramFile),',');
    set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("shotSE", paramFile),',');
    set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    setGridGeometry(shots);

    /* Setting nodes position */

    nodes.n_xline = std::stoi(catchParameter("xNodeNumber", paramFile));
    nodes.n_yline = std::stoi(catchParameter("yNodeNumber", paramFile));
    
    splitted = split(catchParameter("nodeSW", paramFile),',');
    set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("nodeNW", paramFile),',');
    set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = split(catchParameter("nodeSE", paramFile),',');
    set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    setGridGeometry(nodes);

    shotsPath = catchParameter("shotsPath", paramFile);
    nodesPath = catchParameter("nodesPath", paramFile);

    if (reciprocity) setReciprocity();
    if (saveGeometry) exportPositions();

    std::vector<std::string>().swap(splitted);

    /* Importing model and set wavefield volumes */

    V = expand(readBinaryFloat(catchParameter("modelPath", paramFile), nPoints));    

    U_pas = new float[nPointsB]();
    U_pre = new float[nPointsB]();
    U_fut = new float[nPointsB]();
    
    seismogram = new float[nt*nodes.all]();
}

void Acoustic::sourceGenerator()
{
    nsrc = (int) (2.0f * tlag / dt + 1);

    source = new float[nsrc]();

    float pi = 4.0f * atanf(1.0f);
    float fc = fmax / (3.0f * sqrtf(pi));

    int s = (int) (nsrc / 2);

    for (int t = -s; t < s; t++)
    {
        float aux1 = 1.0f - 2.0f * pi * powf(t*dt, 2.0f) * powf(fc, 2.0f) * powf(pi, 2.0f);
        float aux2 = expf(-pi * powf(t*dt, 2.0f)*powf(fc, 2.0f)*powf(pi, 2.0f));

        source[t + s] = aux1 * aux2;
    }
}

void Acoustic::dampingGenerator()
{
    damp1D = new float[nb]();
    damp2D = new float[nb * nb]();
    damp3D = new float[nb * nb * nb]();

    /* 1D damp construction */
    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = expf(-powf(factor * (nb - i), 2.0f)); 
    }

    /* 2D damp construction */
    for(int i = 0; i < nb; i++) 
    {
        damp2D[i + i*nb] = damp1D[i];            

        for (int j = i; j < nb; j++)
        {   
            damp2D[j + i*nb] = damp1D[i];
            damp2D[i + j*nb] = damp1D[i];
        }    
    }

    /* 3D damp construction */
    for (int i  = 0; i < nb; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            damp3D[i + j*nb + j*nb*nb] = damp2D[i + j*nb];
        }

        for(int j = 0; j < nb-1; j++)
        {
            for(int k = j+1; k < nb; k++)
            {
                damp3D[i + j*nb + k*nb*nb] = damp3D[i + j*nb + j*nb*nb];
                damp3D[i + k*nb + j*nb*nb] = damp3D[i + j*nb + j*nb*nb];
            }
        }
    }        
}

void Acoustic::setWaveField()
{
    for (int index = 0; index < nPointsB; index++)
    {
        U_pas[index] = 0.0f;
        U_pre[index] = 0.0f;
        U_fut[index] = 0.0f;
    }
}

void Acoustic::forwardModeling()
{
    for (shotId = 0; shotId < shots.all; shotId++)
    {
        shots.idx = (int)(shots.x[shotId] / dx) + nb;
        shots.idy = (int)(shots.y[shotId] / dy) + nb;
        shots.idz = (int)(shots.z[shotId] / dz) + nb;

        setWaveField();

        # pragma acc enter data copyin(this[0:1], V[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], source[0:nsrc])
        # pragma acc enter data copyin(this[0:1], damp1D[0:nb])
        # pragma acc enter data copyin(this[0:1], damp2D[0:nb*nb])
        # pragma acc enter data copyin(this[0:1], damp3D[0:nb*nb*nb])
        # pragma acc enter data copyin(this[0:1], U_pas[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], U_pre[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], U_fut[0:nPointsB])
        # pragma acc enter data copyin(this[0:1], seismogram[0:nt*nodes.all])
        for (timeStep = 0; timeStep < nt; timeStep++)
        {
            applyWavelet();
            progressMessage();
            wavePropagation();
            dampApplication();
            wavefieldUpdate();
            buildSeismogram();
        }
        # pragma acc exit data delete(V[0:nPointsB], this[0:1])
        # pragma acc exit data delete(source[0:nsrc], this[0:1])
        # pragma acc exit data delete(damp1D[0:nb], this[0:1])
        # pragma acc exit data delete(damp2D[0:nb*nb], this[0:1])
        # pragma acc exit data delete(damp3D[0:nb*nb*nb], this[0:1])
        # pragma acc exit data delete(U_pas[0:nPointsB], this[0:1])
        # pragma acc exit data delete(U_pre[0:nPointsB], this[0:1])
        # pragma acc exit data delete(U_fut[0:nPointsB], this[0:1])
        # pragma acc exit data copyout(seismogram[0:nt*nodes.all], this[0:1])
    
        exportSeismogram();
    }
}

void Acoustic::applyWavelet()
{
    int sId = shots.idz + shots.idx*nzz + shots.idy*nxx*nzz; 

    # pragma acc kernels present(U_pre[0:nPointsB], source[0:nsrc])
    {
        if (timeStep < nsrc) 
            U_pre[sId] += source[timeStep] / (dx * dy * dz);
    }        
}

void Acoustic::wavePropagation()
{
    # pragma acc parallel loop present(V[0:nPointsB], U_pas[0:nPointsB], U_pre[0:nPointsB], U_fut[0:nPointsB])
    for (int index = 0; index < nPointsB; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        if((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4)) 
        {
            float d2_Px2 = (- 9.0f*(U_pre[i + (j-4)*nzz + k*nxx*nzz] + U_pre[i + (j+4)*nzz + k*nxx*nzz])
                        +   128.0f*(U_pre[i + (j-3)*nzz + k*nxx*nzz] + U_pre[i + (j+3)*nzz + k*nxx*nzz])
                        -  1008.0f*(U_pre[i + (j-2)*nzz + k*nxx*nzz] + U_pre[i + (j+2)*nzz + k*nxx*nzz])
                        +  8064.0f*(U_pre[i + (j-1)*nzz + k*nxx*nzz] + U_pre[i + (j+1)*nzz + k*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

            float d2_Py2 = (- 9.0f*(U_pre[i + j*nzz + (k-4)*nxx*nzz] + U_pre[i + j*nzz + (k+4)*nxx*nzz])
                        +   128.0f*(U_pre[i + j*nzz + (k-3)*nxx*nzz] + U_pre[i + j*nzz + (k+3)*nxx*nzz])
                        -  1008.0f*(U_pre[i + j*nzz + (k-2)*nxx*nzz] + U_pre[i + j*nzz + (k+2)*nxx*nzz])
                        +  8064.0f*(U_pre[i + j*nzz + (k-1)*nxx*nzz] + U_pre[i + j*nzz + (k+1)*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

            float d2_Pz2 = (- 9.0f*(U_pre[(i-4) + j*nzz + k*nxx*nzz] + U_pre[(i+4) + j*nzz + k*nxx*nzz])
                        +   128.0f*(U_pre[(i-3) + j*nzz + k*nxx*nzz] + U_pre[(i+3) + j*nzz + k*nxx*nzz])
                        -  1008.0f*(U_pre[(i-2) + j*nzz + k*nxx*nzz] + U_pre[(i+2) + j*nzz + k*nxx*nzz])
                        +  8064.0f*(U_pre[(i-1) + j*nzz + k*nxx*nzz] + U_pre[(i+1) + j*nzz + k*nxx*nzz])
                        - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));

            U_fut[index] = powf(dt, 2.0f) * powf(V[index], 2.0f) * (d2_Px2 + d2_Py2 + d2_Pz2) + 2.0f*U_pre[index] - U_pas[index];
        } 
    }
}    

void Acoustic::dampApplication()
{
    # pragma acc parallel loop present(U_pre[0:nPointsB], U_fut[0:nPointsB], damp1D[0:nb], damp2D[0:nb*nb], damp3D[0:nb*nb*nb])
    for (int index = 0; index < nPointsB; index++) 
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        // Damping 1D ------------------------------------------------------------------------

        if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
        {    
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[i];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[i];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp1D[i];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp1D[i];            
        }

        if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[j];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[j];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp1D[j];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp1D[j];            
        }

        if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp1D[k];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp1D[k];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp1D[k];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp1D[k];            
        }

        // Damping 2D ------------------------------------------------------------------------

        if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];            

            U_pre[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];
            U_fut[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[j + k*nb];            
        }

        if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];            

            U_pre[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];
            U_fut[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp2D[i + k*nb];            
        }

        if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp2D[i + j*nb];            
        }

        // Damping 3D ------------------------------------------------------------------------
        if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        {
            U_pre[i + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + j*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + k*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + j*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[i + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];

            U_pre[(nzz-i-1) + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
            U_fut[(nzz-i-1) + (nxx-j-1)*nzz + (nyy-k-1)*nxx*nzz] *= damp3D[i + j*nb + k*nb*nb];
        }
    }
}

void Acoustic::wavefieldUpdate()
{
    # pragma acc parallel loop present(U_pas[0:nPointsB], U_pre[0:nPointsB], U_fut[0:nPointsB])
    for (int index = 0; index < nPointsB; index++)
    {
        U_pas[index] = U_pre[index];
        U_pre[index] = U_fut[index]; 
    }
}

void Acoustic::buildSeismogram()
{
    for (int receiver = 0; receiver < nodes.all; receiver++)
    {
        nodes.idx = (int)(nodes.x[receiver] / dx) + nb;
        nodes.idy = (int)(nodes.y[receiver] / dy) + nb;
        nodes.idz = (int)(nodes.z[receiver] / dz) + nb;

        int index = nodes.idz + nodes.idx*nzz + nodes.idy*nxx*nzz;

        seismogram[timeStep*nodes.all + receiver] = U_fut[index];
    } 
}

void Acoustic::exportSeismogram()
{
    writeBinaryFloat("seismogram.bin", seismogram, nt * nodes.all);
}

void Acoustic::progressMessage()
{
    if (timeStep % (nt/5) == 0)
    {    
        system("clear");
        std::cout<<"Running shot "<<shotId+1<<" at position: (z = "<<shots.z[shotId]<<", x = "<<shots.x[shotId]<<", y = "<<shots.z[shotId]<<")"<<std::endl;
        std::cout<<"Time step: "<<timeStep*dt<<std::endl;
    }
}
