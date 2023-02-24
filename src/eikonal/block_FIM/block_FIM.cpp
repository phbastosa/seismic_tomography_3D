# include <cmath>

# include <fstream>
# include <iostream>

# include "block_FIM.hpp"

void Block_FIM::prepare_volumes()
{
    padx = (BLOCK_LENGTH - eiko_m.x_samples % BLOCK_LENGTH) % BLOCK_LENGTH;
    pady = (BLOCK_LENGTH - eiko_m.y_samples % BLOCK_LENGTH) % BLOCK_LENGTH;
    padz = (BLOCK_LENGTH - eiko_m.z_samples % BLOCK_LENGTH) % BLOCK_LENGTH;

    S = eiko_m.expand_pad(slowness, padx, pady, padz);

    T = new float[eiko_m.total_samples_b]();

	auto volsize = eiko_m.total_samples_b;

	auto blksize = BLOCK_LENGTH * BLOCK_LENGTH * BLOCK_LENGTH;

	auto nBlkX = eiko_m.x_samples_b / BLOCK_LENGTH;
	auto nBlkY = eiko_m.y_samples_b / BLOCK_LENGTH;
	auto nBlkZ = eiko_m.z_samples_b / BLOCK_LENGTH;
	
    auto blknum = nBlkX * nBlkY * nBlkZ;

	memoryStruct_.xdim = eiko_m.x_samples_b;
	memoryStruct_.ydim = eiko_m.y_samples_b;
	memoryStruct_.zdim = eiko_m.z_samples_b;
    
    memoryStruct_.delta_h = eiko_m.x_spacing;

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

void Block_FIM::solve()
{
    apply_model_mask();
    apply_source_time();
	
    block_FIM_solver(memoryStruct_);
	
	extract_solution();
}

void Block_FIM::apply_model_mask() 
{
	int nx = memoryStruct_.xdim;
	int ny = memoryStruct_.ydim;
	int nz = memoryStruct_.zdim;

	int blklength = memoryStruct_.blklength;

	// make each block to be stored contiguously in 1D memory space
	uint idx = 0;
	
	for(int zStr = 0; zStr < nz; zStr += blklength) 
	{
		for(int yStr = 0; yStr < ny; yStr += blklength) 
		{
			for(int xStr = 0; xStr < nx; xStr += blklength) 
			{
				// for each block
				for(int z = zStr; z < zStr + blklength; z++) 
				{
					for(int y = yStr; y < yStr + blklength; y++) 
					{
						for(int x = xStr; x < xStr + blklength; x++) 
						{
							memoryStruct_.h_spd[idx] = S[z + x*nz + y*nx*nz];
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
	int nx = memoryStruct_.xdim;
	int ny = memoryStruct_.ydim;
	int nz = memoryStruct_.zdim;

    int blklength = memoryStruct_.blklength;

	// make each block to be stored contiguously in 1D memory space
	uint idx = 0;
	uint blk_idx = 0;
	uint list_idx = 0;
	uint nActiveBlock = 0;

    float sx = geometry[shots_type]->shots.x[shot_id];
    float sy = geometry[shots_type]->shots.y[shot_id];
    float sz = geometry[shots_type]->shots.z[shot_id];

    int sidx = (int)(sx / eiko_m.x_spacing);
    int sidy = (int)(sy / eiko_m.y_spacing);
    int sidz = (int)(sz / eiko_m.z_spacing);

    int sid = sidz + sidx*nz + sidy*nx*nz;

    float t0 = memoryStruct_.h_spd[sid] * sqrtf(powf((float)(sidx*eiko_m.x_spacing) - sx, 2.0f) +
                                                powf((float)(sidy*eiko_m.y_spacing) - sy, 2.0f) +
                                                powf((float)(sidz*eiko_m.z_spacing) - sz, 2.0f));

  	for(int zStr = 0; zStr < nz; zStr += blklength) 
	{
    	for(int yStr = 0; yStr < ny; yStr += blklength) 
		{
      		for(int xStr = 0; xStr < nx; xStr += blklength) 
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
	int nx = memoryStruct_.xdim;
	int ny = memoryStruct_.ydim;
	int nz = memoryStruct_.zdim;

	uint idx = 0;

    int blklength = memoryStruct_.blklength;

	for(int zStr = 0; zStr < nz; zStr += blklength) 
	{
		for(int yStr = 0; yStr < ny; yStr += blklength) 
		{
			for(int xStr = 0; xStr < nx; xStr += blklength) 
			{
				// for each block
				for(int z = zStr; z < zStr + blklength; z++) 
				{
					for(int y = yStr; y < yStr + blklength; y++) 
					{
						for(int x = xStr; x < xStr + blklength; x++) 
						{
							T[z + x*nz + y*nx*nz] = memoryStruct_.h_sol[idx];
							
							idx++;
						}
					}
				}
			}
		}
	}

    travel_time = eiko_m.reduce_pad(T, padx, pady, padz);
}

void Block_FIM::destroy()
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