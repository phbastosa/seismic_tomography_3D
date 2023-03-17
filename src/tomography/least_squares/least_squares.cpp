# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "least_squares.hpp"
# include "../../eikonal/classic/classic.hpp"
# include "../../eikonal/block_FIM/block_FIM.hpp"
# include "../../eikonal/accurate_FSM/accurate_FSM.hpp"

void Least_squares::expand_fdm()
{
    int nxx = eikonal->nx + 2;
    int nyy = eikonal->ny + 2;
    int nzz = eikonal->nz + 2;

    // Centering
    for (int z = 1; z < nzz - 1; z++)
    {
        for (int y = 1; y < nyy - 1; y++)
        {
            for (int x = 1; x < nxx - 1; x++)
            {
                T[z + x*nzz + y*nxx*nzz] = eikonal->travel_time[(z - 1) + (x - 1)*eikonal->nz + (y - 1)*eikonal->nx*eikonal->nz];
            }
        }
    }

    // Z direction
    for (int z = 0; z < 1; z++)
    {
        for (int y = 1; y < nyy - 1; y++)
        {
            for (int x = 1; x < nxx - 1; x++)
            {
                T[z + x*nzz + y*nxx*nzz] = eikonal->travel_time[0 + (x - 1)*eikonal->nz + (y - 1)*eikonal->nx*eikonal->nz];
                T[(nzz - z - 1) + x*nzz + y*nxx*nzz] = eikonal->travel_time[(eikonal->nz - 1) + (x - 1)*eikonal->nz + (y - 1)*eikonal->nx*eikonal->nz];
            }
        }
    }

    // X direction
    for (int x = 0; x < 1; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = 1; y < nyy - 1; y++)
            {
                T[z + x*nzz + y*nxx*nzz] = T[z + 1*nzz + y*nxx*nzz];
                T[z + (nxx - x - 1)*nzz + y*nxx*nzz] = T[z + (nxx - 1 - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < 1; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                T[z + x*nzz + y*nxx*nzz] = T[z + x*nzz + 1*nxx*nzz];
                T[z + x*nzz + (nyy - y - 1)*nxx*nzz] = T[z + x*nzz + (nyy - 1 - 1)*nxx*nzz];
            }
        }
    }
}

void Least_squares::info_message()
{
    eikonal->info_message();

    if (iteration == max_iteration)
    { 
        std::cout<<"------- Checking final residuo ------------\n\n";
    }
    else
        std::cout<<"------- Computing iteration "<<iteration+1<<" of "<<max_iteration<<" ------------\n\n";

    if (iteration > 0) std::cout<<"Previous iteration residuo: "<<residuo[iteration-1]<<"\n\n";
}

void Least_squares::set_parameters()
{
    int type = std::stoi(catch_parameter("eikonal_type", parameters));    

    switch (type)
    {
    case 0:
        eikonal = new Classic();
        break;
    
    case 1:
        eikonal = new Block_FIM();
        break;

    case 2:
        eikonal = new Accurate_FSM();
        break;    
    
    default:
        eikonal = new Accurate_FSM();
        break;
    }

    eikonal->parameters = parameters;

    eikonal->set_parameters();

    n_data = eikonal->total_shots * eikonal->total_nodes;

    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));
    lambda = std::stof(catch_parameter("tk_param", parameters));
    tk_order = std::stof(catch_parameter("tk_order", parameters));

    smooth = str2bool(catch_parameter("smooth", parameters)); 
    window = std::stoi(catch_parameter("filter_samples", parameters));
    stdv = std::stof(catch_parameter("standard_deviation",parameters));

    dx_tomo = std::stof(catch_parameter("dx_tomo", parameters));
    dy_tomo = std::stof(catch_parameter("dy_tomo", parameters));
    dz_tomo = std::stof(catch_parameter("dz_tomo", parameters));

    nz_tomo = (int)((eikonal->nz-1) * eikonal->dz / dz_tomo) + 1;    
    nx_tomo = (int)((eikonal->nx-1) * eikonal->dx / dx_tomo) + 1;    
    ny_tomo = (int)((eikonal->ny-1) * eikonal->dy / dy_tomo) + 1;   
    
    n_model = nx_tomo * ny_tomo * nz_tomo;

    obs_data_folder = catch_parameter("dobs_folder", parameters);
    cal_data_folder = catch_parameter("first_arrival_folder", parameters);
    convergency_folder = catch_parameter("convergency_folder", parameters);
    estimated_model_folder = catch_parameter("estimated_models_folder", parameters);

    iteration = 0;

    dobs = new float[n_data]();
    dcal = new float[n_data]();    
    xdm = new float[n_model]();
    dm = new float[eikonal->nPoints]();
    illumination = new float[eikonal->nPoints]();
}

void Least_squares::import_obs_data()
{
    int ptr = 0; 
    
    float * data = new float[eikonal->total_nodes]();

    for (int shot = 0; shot < eikonal->total_shots; shot++)
    {
        read_binary_float(obs_data_folder + "obsData_" + std::to_string(eikonal->total_nodes) + "_samples_shot_" + std::to_string(shot+1) + ".bin", data, eikonal->total_nodes);

        for (int d = ptr; d < ptr + eikonal->total_nodes; d++) 
            dobs[d] = data[d - ptr];

        ptr += eikonal->total_nodes;        
    }

    delete[] data;
}

void Least_squares::forward_modeling()
{
    for (int k = 0; k < eikonal->nPoints; k++) 
        illumination[k] = 0.0f;

    eikonal->prepare_volumes();

    T = new float[(eikonal->nx + 2)*(eikonal->ny+2)*(eikonal->nz+2)]();

    for (int i = 0; i < eikonal->total_shots; i++)
    {
        eikonal->shot_id = i;
        
        info_message();
        
        eikonal->eikonal_equation();
        eikonal->write_first_arrival();
        
        expand_fdm();
        ray_tracing();
    }

    delete[] T; 

    write_binary_float("../outputs/illumination/" + eikonal->name + "_illumination_iteration_" + std::to_string(iteration) + ".bin", illumination, eikonal->nPoints);

    eikonal->destroy_volumes();
}

void Least_squares::ray_tracing()
{
    int nxx = eikonal->nx + 2;
    int nzz = eikonal->nz + 2;

    int sIdz = (int)(eikonal->shots->z[eikonal->shot_id] / dz_tomo);
    int sIdx = (int)(eikonal->shots->x[eikonal->shot_id] / dx_tomo);
    int sIdy = (int)(eikonal->shots->y[eikonal->shot_id] / dy_tomo);

    int sId = sIdz + sIdx*nz_tomo + sIdy*nx_tomo*nz_tomo;     

    float rayStep = 0.2f * eikonal->dz;

    for (int ray_id = 0; ray_id < eikonal->total_nodes; ray_id++)
    {
        float zi = eikonal->nodes->z[ray_id];
        float xi = eikonal->nodes->x[ray_id];
        float yi = eikonal->nodes->y[ray_id];

        std::vector < int > ray_index;

        while (true)
        {
            int i = (int)(zi / eikonal->dz) + 1;
            int j = (int)(xi / eikonal->dx) + 1;
            int k = (int)(yi / eikonal->dy) + 1;

            float dTz = (T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*eikonal->dz);    
            float dTx = (T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*eikonal->dx);    
            float dTy = (T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*eikonal->dy);

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

    for (int shot = 0; shot < eikonal->total_shots; shot++)
    {
        read_binary_float(cal_data_folder + eikonal->name + "_times_nr" + std::to_string(eikonal->total_nodes) + "_shot_" + std::to_string(shot+1) + ".bin", data, eikonal->total_nodes);

        for (int d = ptr; d < ptr + eikonal->total_nodes; d++) 
            dcal[d] = data[d - ptr];

        ptr += eikonal->total_nodes;        
    }

    delete[] data;
}

bool Least_squares::converged()
{
    float r = 0.0f;
    for (int i = 0; i < n_data; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(r));
    
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
    std::cout<<"Solving linear system using Tikhonov regularization with order "+ std::to_string(tk_order) + "\n\n";

    M = n_model;                                  
    N = n_data + n_model - tk_order;                    
    NNZ = vG.size() + (tk_order + 1) * (n_model - tk_order);

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    
    for (int index = 0; index < n_data; index++) 
        B[index] = dobs[index] - dcal[index];

    for (int index = 0; index < vG.size(); index++)
    {
        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index];
    }

    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);

    tk_reg_matrix();
    sparse_cgls();

    delete[] B;
    delete[] iA;
    delete[] jA;
    delete[] vA;
}

void Least_squares::tk_reg_matrix()
{
    int elements = tk_order + 1;
		
    int n = n_model - tk_order;
    int nnz = elements * n;	
    
    int * iL = new int[nnz]();
    int * jL = new int[nnz]();
    float * vL = new float[nnz]();

    if (tk_order == 0)
	{
		for (int index = 0; index < nnz; index++)
		{
			iL[index] = index;
			jL[index] = index;
			vL[index] = 1.0f;
		}
	} 
    else
    {
        int * df = new int[elements]();	
        int * df1 = new int[elements + 1]();
        int * df2 = new int[elements + 1]();
        
        df[0] = -1; df[1] = 1;
        
        for (int index = 1; index < tk_order; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                df2[k] = df[k];
                df1[k + 1] = df[k];

                df[k] = df1[k] - df2[k]; 
            }		 
        }
        
        for (int index = 0; index < n; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                iL[elements*index + k] = index;	
                jL[elements*index + k] = index + k;
                vL[elements*index + k] = df[k];
            }	
        }

        delete[] df;
        delete[] df1;
        delete[] df2;
    }

    for (int index = NNZ - nnz; index < NNZ; index++) 
    {
        iA[index] = n_data + iL[index - (NNZ - nnz)];
        jA[index] = jL[index - (NNZ - nnz)];
        vA[index] = lambda * vL[index - (NNZ - nnz)];        
    }

    delete[] iL;
    delete[] jL;
    delete[] vL;
}

void Least_squares::sparse_cgls()
{
    float a, b, qTq, rTr, rd;
    int cg_max_iteration = 10;

    float * s = new float[N]();
    float * q = new float[N]();
    float * r = new float[M]();
    float * p = new float[M]();

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int i = 0; i < N; i++) 
        s[i] = B[i]; 

    // r = G' * s    
    for (int i = 0; i < NNZ; i++) 
        r[jA[i]] += vA[i] * s[iA[i]];        

    // p = r and x = 0;
    for (int i = 0; i < M; i++) 
    {
        p[i] = r[i]; 
        xdm[i] = 0.0f;
    }

    // q = G * p
    for (int i = 0; i < NNZ; i++) 
        q[iA[i]] += vA[i] * p[jA[i]];        

    for (int i = 0; i < cg_max_iteration; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < N; k++)           // q inner product
            qTq += q[k] * q[k];               // qTq = q' * q

        rTr = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r 

        a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

        for (int k = 0; k < M; k++)           // model atualization
            xdm[k] += a * p[k];               // x = x + a * p

        for (int k = 0; k < N; k++)           // s atualization  
            s[k] -= a * q[k];                 // s = s - a * q 

        rd = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product for division 
            rd += r[k] * r[k];                // rd = r' * r

        for (int k = 0; k < M; k++)           // Zeroing r 
            r[k] = 0.0f;                      // r = 0, for multiplication
        
        for (int k = 0; k < NNZ; k++)         // r atualization 
            r[jA[k]] += vA[k] * s[iA[k]];     // r = G' * s    

        rTr = 0.0f;                
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r

        b = rTr / rd;                         // b = (r' * r) / rd

        for (int k = 0; k < M; k++)          
            p[k] = r[k] + b * p[k];           // p = r + b * p 

        for (int k = 0; k < N; k++) 
            q[k] = 0.0f;                      // q = 0, for multiplication

        for (int k = 0; k < NNZ; k++) 
            q[iA[k]] += vA[k] * p[jA[k]];     // q = G * p   
    }
}

void Least_squares::model_update()
{
    for (int index = 0; index < eikonal->nPoints; index++)
    {
        int k = (int) (index / (eikonal->nx*eikonal->nz));        
        int j = (int) (index - k*eikonal->nx*eikonal->nz) / eikonal->nz;    
        int i = (int) (index - j*eikonal->nz - k*eikonal->nx*eikonal->nz);  

        float x = j*eikonal->dx; 
        float y = k*eikonal->dy; 
        float z = i*eikonal->dz; 

        float x0 = floorf(x/dx_tomo)*dx_tomo;
        float y0 = floorf(y/dy_tomo)*dy_tomo;
        float z0 = floorf(z/dz_tomo)*dz_tomo;

        float x1 = floorf(x/dx_tomo)*dx_tomo + dx_tomo;
        float y1 = floorf(y/dy_tomo)*dy_tomo + dy_tomo;
        float z1 = floorf(z/dz_tomo)*dz_tomo + dz_tomo;

        dm[index] = 0.0f;

        if ((i >= 11) && (i < eikonal->nz - 2) && (j >= 15) && (j < eikonal->nx - 16) && (k >= 15) && (k < eikonal->ny - 16))
        {
            int idz = (int)(z / dz_tomo);
            int idx = (int)(x / dx_tomo);
            int idy = (int)(y / dy_tomo);

            int ind_m = (int)(idz + idx*nz_tomo + idy*nx_tomo*nz_tomo);

            float c000 = xdm[ind_m];                  
            float c001 = xdm[ind_m + 1];
            float c100 = xdm[ind_m + nz_tomo];
            float c101 = xdm[ind_m + 1 + nz_tomo];
            float c010 = xdm[ind_m + nx_tomo*nz_tomo];
            float c011 = xdm[ind_m + 1 + nx_tomo*nz_tomo];
            float c110 = xdm[ind_m + nz_tomo + nx_tomo*nz_tomo];
            float c111 = xdm[ind_m + 1 + nz_tomo + nx_tomo*nz_tomo];  

            float dm_ijk = trilinear(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);

            dm[i + j*eikonal->nz + k*eikonal->nx*eikonal->nz] = dm_ijk;            
        }
    }

    float * update = new float[eikonal->nPoints]();

    if (smooth) 
    {
        gaussian(dm, update, eikonal->nx, eikonal->ny, eikonal->nz, window, stdv);
        
        for (int index = 0; index < eikonal->nPoints; index++)
            dm[index] = update[index];
    }

    for (int index = 0; index < eikonal->nPoints; index++)
    {
        int k = (int) (index / (eikonal->nx*eikonal->nz));        
        int j = (int) (index - k*eikonal->nx*eikonal->nz) / eikonal->nz;    
        int i = (int) (index - j*eikonal->nz - k*eikonal->nx*eikonal->nz);  

        if ((i >= 11) && (i < eikonal->nz - 2) && (j >= 15) && (j < eikonal->nx - 15) && (k >= 15) && (k < eikonal->ny - 15))
        {
            eikonal->slowness[index] += dm[index];
        }

        update[index] = 1.0f / eikonal->slowness[index];
    }
    
    write_binary_float(estimated_model_folder + eikonal->name + "_model_iteration_" + std::to_string(iteration) + ".bin", update, eikonal->nPoints);    

    delete[] update;
}

void Least_squares::export_convergency()
{
    std::ofstream resFile(convergency_folder + eikonal->name + "_convergency.txt", std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) resFile<<residuo[r]<<"\n";

    resFile.close();
}


