# include <cmath>
# include <iostream>
# include <algorithm>

# include "least_squares.hpp"
# include "../../eikonal/classic/classic.hpp"
# include "../../eikonal/block_FIM/block_FIM.hpp"
# include "../../eikonal/accurate_FSM/accurate_FSM.hpp"

void Least_squares::info_message()
{
    eikonal->info_message();

    if (iteration == max_iteration)
    { 
        std::cout<<"------- Checking final residuo ------------\n\n";
    }
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" of "<<max_iteration<<" ------------\n\n";

    if (iteration > 0) std::cout<<"Previous iteration residuo: "<<residuo.back()<<"\n\n";
}

void Least_squares::set_parameters()
{
    int type = std::stoi(catch_parameter("eikonal_type", parameters));    

    Eikonal * eikonal_types[] = 
    {
        new Classic(),
        new Block_FIM(),
        new Accurate_FSM()
    };

    eikonal = eikonal_types[type];
    
    eikonal->parameters = parameters;

    eikonal->set_parameters();
    eikonal->prepare_volumes();

    n_data = eikonal->total_shots * eikonal->total_nodes;

    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));
    lambda = std::stof(catch_parameter("tk_param", parameters));
    tkOrder = std::stof(catch_parameter("tk_order", parameters));

    smooth = str2bool(catch_parameter("smooth", parameters)); 
    window = std::stoi(catch_parameter("filter_samples", parameters));
    stdv = std::stof(catch_parameter("standard_deviation",parameters));

    dx_tomo = std::stof(catch_parameter("dx_tomo", parameters));
    dy_tomo = std::stof(catch_parameter("dy_tomo", parameters));
    dz_tomo = std::stof(catch_parameter("dz_tomo", parameters));

    nz_tomo = (int)(static_cast<float>(eikonal->nz-1) * eikonal->dz / dz_tomo) + 1;    
    nx_tomo = (int)(static_cast<float>(eikonal->nx-1) * eikonal->dx / dx_tomo) + 1;    
    ny_tomo = (int)(static_cast<float>(eikonal->ny-1) * eikonal->dy / dy_tomo) + 1;   
    
    n_model = nx_tomo * ny_tomo * nz_tomo;

    obs_data_folder = catch_parameter("dobs_folder", parameters);
    cal_data_folder = catch_parameter("first_arrival_folder", parameters);
    convergency_folder = catch_parameter("convergency_folder", parameters);
    estimated_model_folder = catch_parameter("estimated_models_folder", parameters);

    iteration = 0;

    residuo.reserve(max_iteration);

    x = new float[n_model]();
    dobs = new float[n_data]();
    dcal = new float[n_data]();    
    dm = new float[eikonal->nPoints]();
    illumination = new float[eikonal->nPoints]();
}

void Least_squares::import_obs_data()
{
    int ptr = 0; 
    
    float * data = new float[eikonal->total_nodes];

    for (int shot = 0; shot < eikonal->total_shots; shot++)
    {
        read_binary_float(obs_data_folder + "obsData_" + std::to_string(eikonal->total_nodes) + "_samples_shot_" + std::to_string(shot+1) + ".bin", data, eikonal->total_nodes);

        for (int d = ptr; d < ptr + eikonal->total_nodes; d++) 
            dobs[d] = data[d - ptr];

        ptr += eikonal->total_nodes;        
    }
}

void Least_squares::forward_modeling()
{
    for (int k = 0; k < eikonal->nPoints; k++) 
        illumination[k] = 0.0f;

    for (int i = 0; i < 1; i++)
    {
        eikonal->shot_id = i;
        
        info_message();
        
        eikonal->eikonal_equation();
        eikonal->write_first_arrival();
        
        ray_tracing();
    }

    write_binary_float("illumination_iteration_" + std::to_string(iteration) + ".bin", illumination, eikonal->nPoints);
}

void Least_squares::ray_tracing()
{
    int sIdz = (int)(eikonal->shots->z[eikonal->shot_id] / dz_tomo);
    int sIdx = (int)(eikonal->shots->x[eikonal->shot_id] / dx_tomo);
    int sIdy = (int)(eikonal->shots->y[eikonal->shot_id] / dy_tomo);

    int sId = sIdz + sIdx*nz_tomo + sIdy*nx_tomo*nz_tomo;     

    float rayStep = 0.2f * (eikonal->dz) / 3.0f;
    
    for (int ray_id = 0; ray_id < eikonal->total_nodes; ray_id++)
    {
        float zi = eikonal->nodes->z[ray_id];
        float xi = eikonal->nodes->x[ray_id];
        float yi = eikonal->nodes->y[ray_id];

        std::vector < int > ray_index;

        while (true)
        {
            int i = (int)(zi / eikonal->dz);
            int j = (int)(xi / eikonal->dx);
            int k = (int)(yi / eikonal->dy);

            float dTz = (eikonal->travel_time[(i+1) + j*eikonal->nz + k*eikonal->nx*eikonal->nz] - eikonal->travel_time[(i-1) + j*eikonal->nz + k*eikonal->nx*eikonal->nz]) / (2.0f*eikonal->dz);    
            float dTx = (eikonal->travel_time[i + (j+1)*eikonal->nz + k*eikonal->nx*eikonal->nz] - eikonal->travel_time[i + (j-1)*eikonal->nz + k*eikonal->nx*eikonal->nz]) / (2.0f*eikonal->dx);    
            float dTy = (eikonal->travel_time[i + j*eikonal->nz + (k+1)*eikonal->nx*eikonal->nz] - eikonal->travel_time[i + j*eikonal->nz + (k-1)*eikonal->nx*eikonal->nz]) / (2.0f*eikonal->dy);

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            zi -= rayStep*dTz / norm; // z ray position atualization    
            xi -= rayStep*dTx / norm; // x ray position atualization   
            yi -= rayStep*dTy / norm; // y ray position atualization   

            int im = (int)(zi / dz_tomo); 
            int jm = (int)(xi / dx_tomo); 
            int km = (int)(yi / dy_tomo); 

            ray_index.push_back(im + jm*nz_tomo + km*nx_tomo*nz_tomo);
            illumination[i + j*eikonal->nz + k*eikonal->nx*eikonal->nz] += rayStep;

            if (ray_index.back() == sId) break;
        }
    
        float finalDist = sqrtf(powf(zi - eikonal->shots->z[eikonal->shot_id],2.0f) + powf(xi - eikonal->shots->x[eikonal->shot_id],2.0f) + powf(yi - eikonal->shots->y[eikonal->shot_id],2.0f));

        std::sort(ray_index.begin(), ray_index.end());

        int current = ray_index[0];
        float distance = rayStep;

        for (int index = 0; index < ray_index.size(); index++)
        {
            if (ray_index[index] == current)
            {
                distance += rayStep;
            }
            else
            {
                vG.push_back(distance);
                jG.push_back(current);
                iG.push_back(ray_id + eikonal->shot_id * eikonal->total_nodes);

                if (current == sId) vG.back() = finalDist;

                distance = rayStep;
                current = ray_index[index];    
            }
        }

        if (current == sId)
        {
            vG.push_back(finalDist);
            jG.push_back(sId);
            iG.push_back(ray_id + eikonal->shot_id * eikonal->total_nodes);
        }
        else 
        {
            vG.push_back(distance);
            jG.push_back(current);
            iG.push_back(ray_id + eikonal->shot_id * eikonal->total_nodes);
        }

        std::vector < int >().swap(ray_index);
    }
}

void Least_squares::import_cal_data()
{
    int ptr = 0; 
    
    float * data = new float[eikonal->total_nodes]();

    for (int shot = 0; shot < 1; shot++)
    {
        read_binary_float(cal_data_folder + eikonal->name + "_times_nr" + std::to_string(eikonal->total_nodes) + "_shot_" + std::to_string(shot+1) + ".bin", data, eikonal->total_nodes);

        for (int d = ptr; d < ptr + eikonal->total_nodes; d++) 
            dcal[d] = data[d - ptr];

        ptr += eikonal->total_nodes;        
    }
}

bool Least_squares::converged()
{
    float r = 0.0f;
    for (int i = 0; i < n_data; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.emplace_back(sqrtf(r));
    
    if (iteration >= max_iteration)
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

void Least_squares::optimization()
{
//     std::cout<<"Solving linear system using Tikhonov regularization with order "+ std::to_string(tkOrder) + "\n\n";

//     sparseMatrix L = Utils::getDerivativeMatrix(mTomo.nPoints, tkOrder);
    
//     int N = shots.all*eikonal->total_nodes + L.n; // Data dimension
//     int M = mTomo.nPoints;             // Model dimension
//     int NNZ = vM.size() + L.nnz;       // Non zero elements

//     int * iA = new int[NNZ]();
//     int * jA = new int[NNZ]();
//     float * vA = new float[NNZ]();

//     for (int index = 0; index < vM.size(); index++)
//     {
//         iA[index] = iM[index];
//         jA[index] = jM[index];
//         vA[index] = vM[index];
//     }

//     for (int index = NNZ - L.nnz; index < NNZ; index++) 
//     {
//         iA[index] = shots.all*eikonal->total_nodes + L.i[index - (NNZ - L.nnz)];
//         jA[index] = L.j[index - (NNZ - L.nnz)];
//         vA[index] = lambda * L.v[index - (NNZ - L.nnz)];        
//     }

//     L.erase();
//     std::vector<  int  >().swap(iM);
//     std::vector<  int  >().swap(jM);
//     std::vector< float >().swap(vM);

//     float * x = new float[M]();
//     float * B = new float[N]();

//     for (int index = 0; index < eikonal->total_nodes*shots.all; index++) 
//         B[index] = dobs[index] - dcal[index];

//     sparse_cgls_gpu(iA, jA, vA, B, x, N, M, NNZ, 10, 1e-6f);

//     std::swap(dm,x);

//     delete[] x;
//     delete[] B;
//     delete[] iA;
//     delete[] jA;
//     delete[] vA;
}

void Least_squares::model_update()
{
    // for (int index = 0; index < ; index++)
    // {
    //     int k = (int) (index / (nx*nz));         // y direction
    //     int j = (int) (index - k*nx*nz) / nz;    // x direction
    //     int i = (int) (index - j*nz - k*nx*nz);  // z direction

    //     float x = j*dx; 
    //     float y = k*dy; 
    //     float z = i*dz; 

    //     float x0 = floorf(x/mTomo.dx)*mTomo.dx;
    //     float y0 = floorf(y/mTomo.dy)*mTomo.dy;
    //     float z0 = floorf(z/mTomo.dz)*mTomo.dz;

    //     float x1 = floorf(x/mTomo.dx)*mTomo.dx + mTomo.dx;
    //     float y1 = floorf(y/mTomo.dy)*mTomo.dy + mTomo.dy;
    //     float z1 = floorf(z/mTomo.dz)*mTomo.dz + mTomo.dz;

    //     if ((i >= 0) && (i < nz - 1) && (j >= 0) && (j < nx - 1) && (k >= 0) && (k < ny - 1))
    //     {
    //         int idz = ((int)(z/mTomo.dz));
    //         int idx = ((int)(x/mTomo.dx));
    //         int idy = ((int)(y/mTomo.dy));

    //         int indS = idz + idx*mTomo.nz + idy*mTomo.nx*mTomo.nz;

    //         float c000 = dm[indS];                  
    //         float c001 = dm[indS + 1];
    //         float c100 = dm[indS + mTomo.nz];
    //         float c101 = dm[indS + 1 + mTomo.nz];
    //         float c010 = dm[indS + mTomo.nx*mTomo.nz];
    //         float c011 = dm[indS + 1 + mTomo.nx*mTomo.nz];
    //         float c110 = dm[indS + mTomo.nz + mTomo.nx*mTomo.nz];
    //         float c111 = dm[indS + 1 + mTomo.nz + mTomo.nx*mTomo.nz];  

    //         float ds_ijk = triLinearInterpolation(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);

    //         dS[(i + nb) + (j + nb)*eikonal->nz + (k + nb)*eikonal->nx*eikonal->nz] = ds_ijk;            
    //     }
    // }

    // if (smooth) 
    // {
    //     dS = gaussianFilterSmoothing(dS, eikonal->nx, nyy, eikonal->nz, standardDeviation, filterSamples);
    // }

    // for (int index = 0; index < nPointsB; index++)
    // {
    //     V[index] = 1.0f / ((1.0f / V[index]) + dS[index]);
    // }

    // float * mm = reduce(V);
    
    // writeBinaryFloat(estimatedPath + "estimatedModel_iteration_" + std::to_string(iteration) + ".bin", mm, nPoints);    

    // delete[] mm;
    // delete[] dS;
}

void Least_squares::export_convergency()
{
//     std::ofstream resFile(residuoPath + "convergency.txt", std::ios::out);
    
//     for (int r = 0; r < residuo.size(); r++) resFile<<residuo[r]<<"\n";

//     resFile.close();
}

void Least_squares::export_estimated_model()
{


}


