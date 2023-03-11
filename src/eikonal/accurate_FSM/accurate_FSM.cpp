# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "accurate_FSM.hpp"

void Accurate_FSM::expand_model()
{
    // Centering
    for (int z = 1; z < nzz - 1; z++)
    {
        for (int y = 1; y < nyy - 1; y++)
        {
            for (int x = 1; x < nxx - 1; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = slowness[(z - 1) + (x - 1)*nz + (y - 1)*nx*nz];
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
                S[z + x*nzz + y*nxx*nzz] = slowness[0 + (x - 1)*nz + (y - 1)*nx*nz];
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = slowness[(nz - 1) + (x - 1)*nz + (y - 1)*nx*nz];
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
                S[z + x*nzz + y*nxx*nzz] = S[z + 1*nzz + y*nxx*nzz];
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - 1 - 1)*nzz + y*nxx*nzz];
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
                S[z + x*nzz + y*nxx*nzz] = S[z + x*nzz + 1*nxx*nzz];
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - 1 - 1)*nxx*nzz];
            }
        }
    }
}

void Accurate_FSM::reduce_model()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        travel_time[z + x*nz + y*nx*nz] = T[(z + 1) + (x + 1)*nzz + (y + 1)*nxx*nzz];
    }
}

void Accurate_FSM::set_parameters()
{
    name = "fsm";

    nx = std::stoi(catch_parameter("x_samples", parameters));
    ny = std::stoi(catch_parameter("y_samples", parameters));
    nz = std::stoi(catch_parameter("z_samples", parameters));

    nPoints = nx * ny * nz;

    dx = std::stof(catch_parameter("x_spacing", parameters));    
    dy = std::stof(catch_parameter("y_spacing", parameters));    
    dz = std::stof(catch_parameter("z_spacing", parameters));    

    slowness = new float[nPoints]();    

    std::string vp_file = catch_parameter("vp_file", parameters);
    read_binary_float(vp_file, slowness, nPoints);

    for (int index = 0; index < nPoints; index++)
        slowness[index] = 1.0f / slowness[index];

    Geometry * stype[] = 
    {
        new Regular(),
        new Circular()
        // new Streamer(),
        // new Crosswell()
    };

    Geometry * ntype[] = 
    {
        new Regular(),
        new Circular()
        // new Streamer(),
        // new Crosswell()
    };

    int shots_type = std::stoi(catch_parameter("shots_geometry_type", parameters)); 
    int nodes_type = std::stoi(catch_parameter("nodes_geometry_type", parameters)); 

    shots = stype[shots_type];
    nodes = ntype[nodes_type];

    shots->set_geometry(parameters, "shots");
    nodes->set_geometry(parameters, "nodes");

    reciprocity = str2bool(catch_parameter("reciprocity", parameters)); 

    if (reciprocity)
    {
        std::swap(shots->x, nodes->x); 
        std::swap(shots->y, nodes->y); 
        std::swap(shots->z, nodes->z); 

        std::swap(shots->all, nodes->all);
    }

    total_shots = shots->all;
    total_nodes = nodes->all;

    export_time_volume = str2bool(catch_parameter("export_time_volume", parameters));
    export_first_arrival = str2bool(catch_parameter("export_first_arrival", parameters));    

    time_volume_folder = catch_parameter("time_volume_folder", parameters);
    first_arrival_folder = catch_parameter("first_arrival_folder", parameters);

    travel_time = new float[nPoints]();
    first_arrival = new float[total_nodes]();
}

void Accurate_FSM::prepare_volumes()
{
    nxx = nx + 2;    
    nyy = ny + 2;    
    nzz = nz + 2;    
    
    nPointsB = nxx * nyy * nzz;

    S = new float[nPointsB]();
    T = new float[nPointsB]();

    expand_model();
}

void Accurate_FSM::inner_sweep()
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = i - sgnvz; 
    int j1 = j - sgnvx; 
    int k1 = k - sgnvy;

    // Get local times of surrounding points
    float tv = T[(i - sgntz) + j*nzz + k*nxx*nzz];
    float te = T[i + (j - sgntx)*nzz + k*nxx*nzz];
    float tn = T[i + j*nzz + (k - sgnty)*nxx*nzz];
    float tev = T[(i - sgntz) + (j - sgntx)*nzz + k*nxx*nzz];
    float ten = T[i + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];
    float tnv = T[(i - sgntz) + j*nzz + (k - sgnty)*nxx*nzz];
    float tnve = T[(i - sgntz) + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];     

    int ijk = i + j*nzz + k*nxx*nzz;

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dz * std::min(S[i1 + std::max(j-1,1)*nzz   + std::max(k-1,1)*nxx*nzz], 
                     std::min(S[i1 + std::max(j-1,1)*nzz   + std::min(k,nyy-1)*nxx*nzz], 
                     std::min(S[i1 + std::min(j,nxx-1)*nzz + std::max(k-1,1)*nxx*nzz],
                              S[i1 + std::min(j,nxx-1)*nzz + std::min(k,nyy-1)*nxx*nzz]))); 

    // X direction
    t1D2 = te + dx * std::min(S[std::max(i-1,1)   + j1*nzz + std::max(k-1,1)*nxx*nzz], 
                     std::min(S[std::min(i,nzz-1) + j1*nzz + std::max(k-1,1)*nxx*nzz],
                     std::min(S[std::max(i-1,1)   + j1*nzz + std::min(k,nyy-1)*nxx*nzz], 
                              S[std::min(i,nzz-1) + j1*nzz + std::min(k,nyy-1)*nxx*nzz])));

    // Y direction
    t1D3 = tn + dy * std::min(S[std::max(i-1,1)   + std::max(j-1,1)*nzz   + k1*nxx*nzz], 
                     std::min(S[std::max(i-1,1)   + std::min(j,nxx-1)*nzz + k1*nxx*nzz],
                     std::min(S[std::min(i,nzz-1) + std::max(j-1,1)*nzz   + k1*nxx*nzz], 
                              S[std::min(i,nzz-1) + std::min(j,nxx-1)*nzz + k1*nxx*nzz])));

    t1D = std::min(t1D1, std::min(t1D2, t1D3));

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[i1 + j1*nzz + std::max(k-1,1)*nxx*nzz], S[i1 + j1*nzz + std::min(k, nyy-1)*nxx*nzz]);
    
    if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[i1 + std::max(j-1,1)*nzz + k1*nxx*nzz], S[i1 + std::min(j,nxx-1)*nzz + k1*nxx*nzz]);

    if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*dz2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = std::min(S[std::max(i-1,1) + j1*nzz + k1*nxx*nzz],S[std::min(i,nzz-1) + j1*nzz + k1*nxx*nzz]);

    if((te < tn + dy*Sref) && (tn < te + dx*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*dx2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
    }

    t2D = std::min(t2D1, std::min(t2D2, t2D3));

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6;

    Sref = S[i1 + j1*nzz + k1*nxx*nzz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (std::min(t1D,t2D) > std::max(tv,std::max(te, tn)))
    {
        t2 = 9.0f*Sref*Sref*dsum;
        
        t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
            
            t3D = (t1 + sqrtf(t2 - t3)) / dsum;
        }
    }
   
    T[ijk] = std::min(T[ijk], std::min(t1D, std::min(t2D, t3D)));
}

void Accurate_FSM::init_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = std::max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = std::max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = std::max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Accurate_FSM::full_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North 
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Accurate_FSM::info_message()
{
    system("clear");

    std::cout<<"3D eikonal equation solver\n\n";
    
    std::cout<<"Total x model length = "<<(nx-1)*dx<<" m\n";
    std::cout<<"Total Y model length = "<<(ny-1)*dy<<" m\n";
    std::cout<<"Total Z model length = "<<(nz-1)*dz<<" m\n\n";
    
    if (reciprocity)
        std::cout<<"Reciprocity = True\n\n";
    else
        std::cout<<"Reciprocity = False\n\n";

    std::cout<<"Shot "<<shot_id+1<<" of "<<shots->all<<"\n";

    std::cout<<"Position (x,y,z) = ("<<shots->x[shot_id]<<", "<<shots->y[shot_id]<<", "<<shots->z[shot_id]<<") m\n\n";

    std::cout<<"Solving eikonal equation with the \033[32mNoble, Gesret and Belayouni (2014)\033[0;0m formulation\n\n";
}

void Accurate_FSM::write_time_volume()
{
    if (export_time_volume)        
        write_binary_float(time_volume_folder + name + "_eikonal_nz" + std::to_string(nz) + "_nx" + std::to_string(nx) + "_ny" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin", travel_time, nPoints);
}

void Accurate_FSM::write_first_arrival()
{
    if (export_first_arrival) 
    {   
        for (int r = 0; r < total_nodes; r++)
        {
            float x = nodes->x[r];
            float y = nodes->y[r];
            float z = nodes->z[r];

            float x0 = floorf(x / dx) * dx;
            float y0 = floorf(y / dy) * dy;
            float z0 = floorf(z / dz) * dz;

            float x1 = floorf(x / dx) * dx + dx;
            float y1 = floorf(y / dy) * dy + dy;
            float z1 = floorf(z / dz) * dz + dz;

            int id = ((int)(z / dz)) + ((int)(x / dx))*nz + ((int)(y / dy))*nx*nz;

            float c000 = travel_time[id];
            float c001 = travel_time[id + 1];
            float c100 = travel_time[id + nz]; 
            float c101 = travel_time[id + 1 + nz]; 
            float c010 = travel_time[id + nx*nz]; 
            float c011 = travel_time[id + 1 + nx*nz]; 
            float c110 = travel_time[id + nz + nx*nz]; 
            float c111 = travel_time[id + 1 + nz + nx*nz];

            first_arrival[r] = trilinear(c000,c001,c100,c101,c010,c011,c110,c111,x0,x1,y0,y1,z0,z1,x,y,z);        
        }

        write_binary_float(first_arrival_folder + name + "_times_nr" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot_id+1) + ".bin", first_arrival, total_nodes);
    }
}

void Accurate_FSM::eikonal_equation()
{
    int nb = 1;

    float sx = shots->x[shot_id];
    float sy = shots->y[shot_id];
    float sz = shots->z[shot_id];

    sidx = (int)(sx / dx) + nb; 
    sidy = (int)(sy / dy) + nb;
    sidz = (int)(sz / dz) + nb;
    
    int sId = sidz + sidx*nzz + sidy*nxx*nzz;     

    for (int index = 0; index < nPointsB; index++)
        T[index] = 1e6f;

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));

    T[sId + nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1 + nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 + nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1 - nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf((sidy-nb)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    
    T[sId + 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId - 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-nb)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    
    T[sId + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    T[sId - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf((sidz-nb)*dz - sz, 2.0f));
    
    T[sId + 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));
    T[sId + 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)+1)*dz - sz, 2.0f));

    T[sId - 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)+1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)+1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));
    T[sId - 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-nb)-1)*dx - sx, 2.0f) + powf(((sidy-nb)-1)*dy - sy, 2.0f) + powf(((sidz-nb)-1)*dz - sz, 2.0f));

    dzi = 1.0f / dz;
    dxi = 1.0f / dx;
    dyi = 1.0f / dy;

    dz2i = 1.0f / (dz*dz);
    dx2i = 1.0f / (dx*dx);
    dy2i = 1.0f / (dy*dy);

    dz2dx2 = dz2i * dx2i;
    dz2dy2 = dz2i * dy2i;
    dx2dy2 = dx2i * dy2i;

    dsum = dz2i + dx2i + dy2i;

    init_sweep();
    full_sweep();

    reduce_model();
}

void Accurate_FSM::destroy_volumes()
{
    delete[] S;
    delete[] T;
}
