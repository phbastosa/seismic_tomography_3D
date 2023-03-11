# include <cmath>
# include <fstream> 
# include <iostream>
# include <algorithm>

# include "block_FIM.hpp"

void Block_FIM::expand_model()
{
    // Centering
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = slowness[z + x*nz + y*nx*nz];
            }
        }
    }

    // Z direction
    for (int z = 0; z < padz; z++)
    {
        for (int y = 0; y < nyy - pady; y++)
        {
            for (int x = 0; x < nxx - padx; x++)
            {
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = slowness[(nz - 1) + x*nz + y*nx*nz];
            }
        }
    }

    // X direction
    for (int x = 0; x < padx; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = 0; y < nyy - pady; y++)
            {
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - padx - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < pady; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - pady - 1)*nxx*nzz];
            }
        }
    }
}

void Block_FIM::reduce_model()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        travel_time[z + x*nz + y*nx*nz] = T[z + x*nzz + y*nxx*nzz];
    }
}

void Block_FIM::set_parameters()
{
    name = "fim";

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

void Block_FIM::prepare_volumes()
{
    padx = (BLOCK_LENGTH - nx % BLOCK_LENGTH) % BLOCK_LENGTH;
    pady = (BLOCK_LENGTH - ny % BLOCK_LENGTH) % BLOCK_LENGTH;
    padz = (BLOCK_LENGTH - nz % BLOCK_LENGTH) % BLOCK_LENGTH;

    nxx = nx + padx;    
    nyy = ny + pady;    
    nzz = nz + padz;    
    
    nPointsB = nxx * nyy * nzz;

    S = new float[nPointsB]();
    T = new float[nPointsB]();

    expand_model();

	auto volsize = nPointsB;

	auto blksize = BLOCK_LENGTH * BLOCK_LENGTH * BLOCK_LENGTH;

	auto nBlkX = nxx / BLOCK_LENGTH;
	auto nBlkY = nyy / BLOCK_LENGTH;
	auto nBlkZ = nzz / BLOCK_LENGTH;
	
    auto blknum = nBlkX * nBlkY * nBlkZ;

	memoryStruct_.xdim = nxx;
	memoryStruct_.ydim = nyy;
	memoryStruct_.zdim = nzz;
    
    memoryStruct_.delta_h = dx;

	memoryStruct_.nIter = 10;
	memoryStruct_.blknum = static_cast<uint>(blknum);
    memoryStruct_.volsize = static_cast<uint>(volsize);
	memoryStruct_.blksize = static_cast<uint>(blksize);
	memoryStruct_.blklength = BLOCK_LENGTH;

	memoryStruct_.h_spd = new float[volsize]();
	memoryStruct_.h_mask = new bool[volsize]();
	memoryStruct_.h_sol = new float[volsize]();   // initial solution
	
    memoryStruct_.h_list = new uint[blknum]();    // linear list contains active block indices
	memoryStruct_.h_listed = new bool[blknum]();  // whether block is added to the list
	memoryStruct_.h_listVol = new bool[blknum](); // volume list shows active/nonactive of corresponding block
}

void Block_FIM::info_message()
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

    std::cout<<"Solving eikonal equation with the \033[32mJeong & Whitaker (2008)\033[0;0m formulation\n\n";
}

void Block_FIM::eikonal_equation()
{
    apply_model_mask();
    apply_source_time();
	
    block_FIM_solver(memoryStruct_);
	
	extract_solution();
}

void Block_FIM::write_time_volume()
{
    if (export_time_volume)        
        write_binary_float(time_volume_folder + name + "_eikonal_nz" + std::to_string(nz) + "_nx" + std::to_string(nx) + "_ny" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin", travel_time, nPoints);
}

void Block_FIM::write_first_arrival()
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

void Block_FIM::apply_model_mask()
{
	int blklength = memoryStruct_.blklength;

	// make each block to be stored contiguously in 1D memory space
	uint idx = 0;
	
	for(int zStr = 0; zStr < nzz; zStr += blklength) 
	{
		for(int yStr = 0; yStr < nyy; yStr += blklength) 
		{
			for(int xStr = 0; xStr < nxx; xStr += blklength) 
			{
				// for each block
				for(int z = zStr; z < zStr + blklength; z++) 
				{
					for(int y = yStr; y < yStr + blklength; y++) 
					{
						for(int x = xStr; x < xStr + blklength; x++) 
						{
							memoryStruct_.h_spd[idx] = S[z + x*nzz + y*nxx*nzz];
							memoryStruct_.h_mask[idx] = true;
							idx++;
						}
					}
				}
			}
		}
	}
}

void Block_FIM::apply_source_time()
{    
    int blklength = memoryStruct_.blklength;

	// make each block to be stored contiguously in 1D memory space
	uint idx = 0;
	uint blk_idx = 0;
	uint list_idx = 0;
	uint nActiveBlock = 0;

    float sx = shots->x[shot_id];
    float sy = shots->y[shot_id];
    float sz = shots->z[shot_id];

    int sidx = (int)(sx / dx);
    int sidy = (int)(sy / dy);
    int sidz = (int)(sz / dz);

    int sid = sidz + sidx*nzz + sidy*nxx*nzz;

    float t0 = S[sid] * sqrtf(powf((float)(sidx*dx) - sx, 2.0f) +
                              powf((float)(sidy*dy) - sy, 2.0f) +
                              powf((float)(sidz*dz) - sz, 2.0f));

  	for(int zStr = 0; zStr < nzz; zStr += blklength) 
	{
    	for(int yStr = 0; yStr < nyy; yStr += blklength) 
		{
      		for(int xStr = 0; xStr < nxx; xStr += blklength) 
			{
        		// for each block
        		bool isSeedBlock = false;

        		for(int z = zStr; z < zStr + blklength; z++) 
				{
          			for(int y = yStr; y < yStr + blklength; y++) 
					{
            			for(int x = xStr; x < xStr + blklength; x++) 
						{
                            memoryStruct_.h_sol[idx] = INF;

                            if (x == sidx && y == sidy && z == sidz) 
                            {
                                memoryStruct_.h_sol[idx] = t0;

                                isSeedBlock = true;
                            }
              	
							idx++;
            			}
          			}
        		}
        
        		if(isSeedBlock) 
				{          			
					memoryStruct_.h_listVol[blk_idx] = true;
          			memoryStruct_.h_listed[blk_idx] = true;
          			memoryStruct_.h_list[list_idx] = blk_idx;
          			
                    list_idx++;
          			nActiveBlock++;
        		} 
				else 
				{
          			memoryStruct_.h_listVol[blk_idx] = false;
          			memoryStruct_.h_listed[blk_idx] = false;
        		}
        		
				blk_idx++;
      		}
    	}
  	}

	memoryStruct_.nActiveBlock = nActiveBlock;
}

void Block_FIM::extract_solution()
{
	uint idx = 0;

    int blklength = memoryStruct_.blklength;

	for(int zStr = 0; zStr < nzz; zStr += blklength) 
	{
		for(int yStr = 0; yStr < nyy; yStr += blklength) 
		{
			for(int xStr = 0; xStr < nxx; xStr += blklength) 
			{
				// for each block
				for(int z = zStr; z < zStr + blklength; z++) 
				{
					for(int y = yStr; y < yStr + blklength; y++) 
					{
						for(int x = xStr; x < xStr + blklength; x++) 
						{
							T[z + x*nzz + y*nxx*nzz] = memoryStruct_.h_sol[idx];
							
							idx++;
						}
					}
				}
			}
		}
	}

    reduce_model();
}

void Block_FIM::destroy_volumes()
{
    delete[] S;
    delete[] T;
    delete[] memoryStruct_.h_spd;
    delete[] memoryStruct_.h_mask;
    delete[] memoryStruct_.h_sol;
    delete[] memoryStruct_.h_list;
    delete[] memoryStruct_.h_listed;
    delete[] memoryStruct_.h_listVol;
}