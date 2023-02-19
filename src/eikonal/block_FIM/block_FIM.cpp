# include <cmath>

# include <fstream>
# include <iostream>

# include "block_FIM.cuh"

void Block_FIM::solve()
{
	initialization();
	
	// apply_source();

	// block_FIM_solver(this->memoryStruct_, this->verbose_);
	
	// this->get_solution();
}

void Block_FIM::initialization() 
{
    model.new_x_samples = model.x_samples + (WARP - model.x_samples % WARP) % WARP;
    model.new_y_samples = model.y_samples + (WARP - model.y_samples % WARP) % WARP;
    model.new_z_samples = model.z_samples + (WARP - model.z_samples % WARP) % WARP;

    model.new_x_spacing = model.x_spacing * model.new_x_samples / model.x_samples;
    model.new_y_spacing = model.y_spacing * model.new_y_samples / model.y_samples;
    model.new_z_spacing = model.z_spacing * model.new_z_samples / model.z_samples;

    slowness = model.resize(slowness);




	auto volSize = nx * ny * nz;

	auto blkSize = BLOCK_LENGTH * BLOCK_LENGTH * BLOCK_LENGTH;

	auto nBlkX = nx / BLOCK_LENGTH;
	auto nBlkY = ny / BLOCK_LENGTH;
	auto nBlkZ = nz / BLOCK_LENGTH;
	
    auto blockNum = nBlkX * nBlkY * nBlkZ;

	memoryStruct_.xdim = 
	memoryStruct_.ydim = ny + (BLOCK_LENGTH - ny % BLOCK_LENGTH) % BLOCK_LENGTH;
	memoryStruct_.zdim = nz + (BLOCK_LENGTH - nz % BLOCK_LENGTH) % BLOCK_LENGTH;
    
    memoryStruct_.delta_h = model.x_spacing;

	memoryStruct_.nIter = 10;
	memoryStruct_.blknum = blockNum;
    memoryStruct_.volsize = volSize;
	memoryStruct_.blksize = blkSize;
	memoryStruct_.blklength = BLOCK_LENGTH;

	memoryStruct_.h_sol = new float[volSize];     // initial solution
	memoryStruct_.h_list = new uint[blockNum];    // linear list contains active block indices
	memoryStruct_.h_listed = new bool[blockNum];  // whether block is added to the list
	memoryStruct_.h_listVol = new bool[blockNum]; // volume list shows active/nonactive of corresponding block
	memoryStruct_.blockOrder = new int[blockNum]; 
}

// void Block_FIM::set_attribute_mask() 
// {
// 	uint volSize = memoryStruct_.volsize;

// 	int nx = memoryStruct_.xdim;
// 	int ny = memoryStruct_.ydim;
// 	int nz = memoryStruct_.zdim;

// 	int blklength = memoryStruct_.blklength;

// 	// create host memory
// 	memoryStruct_.h_spd = new float[volSize];
// 	memoryStruct_.h_mask = new bool[volSize];

// 	// copy input volume to host memory
// 	// make each block to be stored contiguously in 1D memory space
// 	uint idx = 0;
	
// 	for(int zStr = 0; zStr < nz; zStr += blklength) 
// 	{
// 		for(int yStr = 0; yStr < ny; yStr += blklength) 
// 		{
// 			for(int xStr = 0; xStr < nx; xStr += blklength) 
// 			{
// 				// for each block
// 				for(int z = zStr; z < zStr + blklength; z++) 
// 				{
// 					for(int y = yStr; y < yStr + blklength; y++) 
// 					{
// 						for(int x = xStr; x < xStr + blklength; x++) 
// 						{
// 							memoryStruct_.h_spd[idx] = slowness[z + x*nz + y*nz*ny];
// 							memoryStruct_.h_mask[idx] = true;
// 							idx++;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// }

// void Block_FIM::useSeeds() 
// {  
// 	if (this->verbose_) std::cout<<"Loading seed volume..."<<std::endl;
  
// 	uint volSize, blockNum;

// 	int nx = this->memoryStruct_.xdim;
// 	int ny = this->memoryStruct_.ydim;
// 	int nz = this->memoryStruct_.zdim;
	
//     int blklength = this->memoryStruct_.blklength;

// 	uint volSize = this->memoryStruct_.volsize;
// 	uint blockNum = this->memoryStruct_.blknum;

// 	// copy input volume to host memory
// 	// make each block to be stored contiguously in 1D memory space
// 	uint idx = 0;
// 	uint blk_idx = 0;
// 	uint list_idx = 0;
// 	uint nActiveBlock = 0;

//     float sx = geometry[shots_type]->shots.x[shot_id];
//     float sy = geometry[shots_type]->shots.y[shot_id];
//     float sz = geometry[shots_type]->shots.z[shot_id];

//     int sidx = (sx / model.x_spacing);
//     int sidy = (sy / model.x_spacing);
//     int sidz = (sz / model.x_spacing);

//     t0 = ;

//   	for(int zStr = 0; zStr < nz; zStr += blklength) 
// 	{
//     	for(int yStr = 0; yStr < ny; yStr += blklength) 
// 		{
//       		for(int xStr = 0; xStr < nx; xStr += blklength) 
// 			{
//         		// for each block
//         		bool isSeedBlock = false;

//         		for(int z = zStr; z < zStr + blklength; z++) 
// 				{
//           			for(int y = yStr; y < yStr + blklength; y++) 
// 					{
//             			for(int x = xStr; x < xStr + blklength; x++) 
// 						{
//               				this->memoryStruct_.h_sol[idx] = INF;
                         
//                             if (x == sx && y == sy && z == sz) 
//                             {
//                                 this->memoryStruct_.h_sol[idx] = 0;
                                
//                                 isSeedBlock = true;
//                             }
              	
// 							idx++;
//             			}
//           			}
//         		}
        
//         		if(isSeedBlock) 
// 				{          			
// 					this->memoryStruct_.h_listVol[blk_idx] = true;
//           			this->memoryStruct_.h_listed[blk_idx] = true;
//           			this->memoryStruct_.h_list[list_idx] = blk_idx;
          			
//                     list_idx++;
//           			nActiveBlock++;
//         		} 
// 				else 
// 				{
//           			this->memoryStruct_.h_listVol[blk_idx] = false;
//           			this->memoryStruct_.h_listed[blk_idx] = false;
//         		}
        		
// 				blk_idx++;
//       		}
//     	}
//   	}

// 	this->memoryStruct_.nActiveBlock = nActiveBlock;
// }

// void Block_FIM::get_solution() 
// {
// 	// put the data where it belongs in the grand scheme of data!
// 	this->answer_ = std::vector<std::vector<std::vector<float>>>(nx, 
// 					std::vector<std::vector<float>>(ny, 
// 					std::vector<float>(nz, 0)));

// 	for (size_t blockID = 0; blockID < this->memoryStruct_.blknum; blockID++) 
// 	{
// 		size_t baseAddr = blockID * this->memoryStruct_.blksize;
// 		size_t xgridlength = this->memoryStruct_.xdim/BLOCK_LENGTH;
// 		size_t ygridlength = this->memoryStruct_.ydim/BLOCK_LENGTH;
		
// 		// compute block index
// 		size_t bx = blockID%xgridlength;
// 		size_t tmpIdx = (blockID - bx)/xgridlength;
// 		size_t by = tmpIdx%ygridlength;
// 		size_t bz = (tmpIdx-by)/ygridlength;
		
// 		//translate back to real space
// 		for (int k = 0; k < BLOCK_LENGTH; k++) 
// 		{
// 			for (int j = 0; j < BLOCK_LENGTH; j++) 
// 			{
// 				for(int i = 0; i < BLOCK_LENGTH; i++) 
// 				{
// 					float d = this->memoryStruct_.h_sol[baseAddr + k * BLOCK_LENGTH * BLOCK_LENGTH + j * BLOCK_LENGTH + i];
				
// 					if ((i + bx * BLOCK_LENGTH) < this->width_ && (j + by * BLOCK_LENGTH) < this->height_ && (k + bz * BLOCK_LENGTH) < this->depth_) 
// 					{
// 						this->answer_[(i + bx * BLOCK_LENGTH)][(j + by * BLOCK_LENGTH)][k + bz * BLOCK_LENGTH] = d;
// 					}
// 				}
// 			}
// 		}
// 	}

//   	for(size_t k = 0; k < nz; k++) 
// 	{
//     	for(size_t j = 0; j < ny; j++) 
// 		{
//       		for(size_t i = 0; i < nx; i++) 
// 			{
//         		travel_time[k + i*nz + j*nz*ny] = this->answer_[i][j][k];  
//       		}
//     	}
//   	}

// }


