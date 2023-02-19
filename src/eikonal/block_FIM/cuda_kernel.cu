# include <cmath> 
# include <vector>
# include <cassert>
# include <iostream>

# include "cuda_kernel.cuh"

void cuda_safe_call(cudaError_t error)
{
    if(error != cudaSuccess)
    {
        std::cout<<"CUDA error! "<<error<<"\n";
        exit(EXIT_FAILURE);
    }
}

void block_FIM_solver(CUDAMEMSTRUCT &cmem, bool verbose)
{
    int deviceID; 
    cudaGetDevice(&deviceID);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceID);

	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);

	std::cout << "Device id : "<<deviceID<<", name : "<<deviceProp.name<<"\n";	
    std::cout << "Total Memory : " << totalMem / (1024 * 1024) << "MB" << "\n";
    std::cout << "Free Memory  : " << freeMem / (1024 * 1024) << "MB" << "\n";

    int xdim = cmem.xdim;
    int ydim = cmem.ydim;
    int zdim = cmem.zdim;
    float dh = cmem.delta_h;    

    uint volsize = cmem.volsize;
    uint blknum = cmem.blknum;

    int nIter = cmem.nIter;
    uint nActiveBlock = cmem.nActiveBlock; // active list

    float * d_spd;
    float * d_sol;
    float * h_sol;
    float * t_sol;

    uint * d_list;
    bool * d_listVol;

    bool * d_con;
    bool * d_mask;

    // copy so that original value should not be modified
    uint *h_list = (uint*) malloc(blknum*sizeof(uint));
    bool *h_listed = (bool*) malloc(blknum*sizeof(bool));
    bool *h_listVol = (bool*) malloc(blknum*sizeof(bool));

    // initialization
    memcpy(h_list, cmem.h_list, blknum*sizeof(uint));
    memcpy(h_listed, cmem.h_listed, blknum*sizeof(bool));
    memcpy(h_listVol, cmem.h_listVol, blknum*sizeof(bool));

  	// create host/device memory using CUDA mem functions

	cuda_safe_call(cudaMalloc((void**)&(d_spd), volsize*sizeof(float)));
	cuda_safe_call(cudaMalloc((void**)&(d_sol), volsize*sizeof(float)));
	cuda_safe_call(cudaMalloc((void**)&(t_sol), volsize*sizeof(float))); 
	cuda_safe_call(cudaMalloc((void**)&(d_con), volsize*sizeof(bool)));  
	cuda_safe_call(cudaMalloc((void**)&(d_list), blknum*sizeof(uint)));
	cuda_safe_call(cudaMalloc((void**)&(d_listVol), blknum*sizeof(bool)));
	cuda_safe_call(cudaMalloc((void**)&(d_mask), volsize*sizeof(bool)));

	cuda_safe_call(cudaMemcpy(d_spd, cmem.h_spd, volsize*sizeof(float), cudaMemcpyHostToDevice));
	cuda_safe_call(cudaMemcpy(d_mask, cmem.h_mask, volsize*sizeof(bool), cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(d_listVol, h_listVol, blknum*sizeof(bool), cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(d_sol, h_sol, volsize*sizeof(float), cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(t_sol, h_sol, volsize*sizeof(float), cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemset(d_con, 1, volsize*sizeof(bool)));

    // set dimension of block and entire grid size
    dim3 dimBlock(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH);
    dim3 dimEntireGrid(blknum);
    dim3 dimGrid(nActiveBlock);

    int nTotalIter = 0;

    std::vector<int> sourceList;
    sourceList.push_back((zdim/2)*ydim*xdim + (ydim/2)*xdim + (xdim/2));

    uint nTotalBlockProcessed = 0;

    while(nActiveBlock > 0)
    {
        assert(nActiveBlock < 4294967295);

        nTotalBlockProcessed += nActiveBlock;

        nTotalIter++;

        // 1. run solver on current active tiles
        
        dimGrid.y = (unsigned int)floorf((float)(nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceilf((float)nActiveBlock/(float)dimGrid.y);

        cuda_safe_call(cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice));
        
        run_solver<<<dimGrid,dimBlock>>>(d_spd, d_mask, d_sol, t_sol, d_con, d_list, xdim, ydim, zdim, dh, nIter, nActiveBlock);
        
        cuda_safe_call(cudaGetLastError());       
        cudaDeviceSynchronize();   

        // 2. reduction (only active tiles)

        run_reduction<<<dimGrid,dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2)>>>(d_con, d_listVol, d_list, nActiveBlock);
        
        cuda_safe_call(cudaGetLastError());
        cudaDeviceSynchronize();

        // 3. check neighbor tiles of converged tile
        // Add any active block of neighbor of converged block is inserted
        // to the list

        cuda_safe_call(cudaMemcpy(h_listVol, d_listVol, blknum*sizeof(bool), cudaMemcpyDeviceToHost));

        uint nBlkX = xdim/BLOCK_LENGTH;
        uint nBlkY = ydim/BLOCK_LENGTH;
        uint nOldActiveBlock = nActiveBlock;

        for(uint i = 0; i < nOldActiveBlock; i++)
        {
            // check 6-neighbor of current active tile
            uint currBlkIdx = h_list[i];

            if(!h_listVol[currBlkIdx]) // not active : converged
            {
                uint nb[6];
                nb[0] = (currBlkIdx < nBlkX*nBlkY) ? currBlkIdx : (currBlkIdx - nBlkX*nBlkY);  //tp
                nb[1] = ((currBlkIdx + nBlkX*nBlkY) >= blknum) ? currBlkIdx : (currBlkIdx + nBlkX*nBlkY); //bt
                nb[2] = (currBlkIdx < nBlkX) ? currBlkIdx : (currBlkIdx - nBlkX); //up
                nb[3] = ((currBlkIdx + nBlkX) >= blknum) ? currBlkIdx : (currBlkIdx + nBlkX); //dn
                nb[4] = (currBlkIdx%nBlkX == 0) ? currBlkIdx : currBlkIdx-1; //lf
                nb[5] = ((currBlkIdx+1)%nBlkX == 0) ? currBlkIdx : currBlkIdx+1; //rt

                for(int nbIdx = 0; nbIdx < 6; nbIdx++)
                {
                    uint currIdx = nb[nbIdx];

                    if(!h_listed[currIdx])
                    {
                        h_listed[currIdx] = true;
                        h_list[nActiveBlock++] = currIdx;
                    }
                }
            }
        }

        cudaDeviceSynchronize();
    
        // 4. run solver only once for neighbor blocks of converged block
        // current active list contains active blocks and neighbor blocks of
        // any converged blocks.

        // update grid dimension because nActiveBlock is changed
        dimGrid.y = (unsigned int)floor(((float)nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceil((float)nActiveBlock/(float)dimGrid.y);

        cuda_safe_call(cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice));
        
        run_check_neighbor<<< dimGrid, dimBlock >>>(d_spd, d_mask, t_sol, d_sol, d_con, d_list, xdim, ydim, zdim, dh, nOldActiveBlock, nActiveBlock);
        
        cuda_safe_call(cudaGetLastError());
        cudaDeviceSynchronize();

        // 5. reduction

        run_reduction<<<dimGrid,dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2)>>>(d_con, d_listVol, d_list, nActiveBlock);
        cuda_safe_call(cudaGetLastError());
        cudaDeviceSynchronize();

        // 6. update active list
        // read back active volume from the device and add
        // active block to active list on the host memory

        nActiveBlock = 0;
        cuda_safe_call(cudaMemcpy(h_listVol, d_listVol, blknum*sizeof(bool), cudaMemcpyDeviceToHost));

        for(uint i=0; i<blknum; i++)
        {
            if(h_listVol[i]) // true : active block (not converged)
            {
                h_listed[i] = true;
                h_list[nActiveBlock++] = i;
            }
            else
            { 
                h_listed[i] = false;
            }
        }

        cudaDeviceSynchronize();
    }

	cuda_safe_call(cudaMemcpy(cmem.h_sol, d_sol, volsize*sizeof(float), cudaMemcpyDeviceToHost));

    free(h_list);
    free(h_listed);
    free(h_listVol);

    cuda_safe_call(cudaFree(d_spd));
    cuda_safe_call(cudaFree(d_sol));
    cuda_safe_call(cudaFree(t_sol));  // temp solution for ping-pong
    cuda_safe_call(cudaFree(d_con));  // convergence volume
    cuda_safe_call(cudaFree(d_list));
    cuda_safe_call(cudaFree(d_listVol));
    cuda_safe_call(cudaFree(d_mask));
}

__device__ float get_time_eikonal(float a, float b, float c, float h, float s)
{
	float ret, tmp;

	// a > b > c
	if(a < b) { tmp = a; a = b; b = tmp; }
	if(b < c) { tmp = b; b = c; c = tmp; }
	if(a < b) { tmp = a; a = b; b = tmp; }

	ret = INF;

	if(c < INF)
	{
		ret = c + h*s;
		
		if(ret > b) 
		{	
			tmp = ((b+c) + sqrtf(2.0f*s*s*h*h - (b-c)*(b-c)))*0.5f;
		
			if(tmp > b) ret = tmp; 

			if(ret > a)	
			{				
                tmp = (a+b+c)/3.0f + sqrtf(2.0f*(a*(b-a) + b*(c-b) + c*(a-c)) + 3.0f*s*s*h*h) / 3.0f;

				if(tmp > a) ret = tmp;
			}
		}
	}

	return ret;
}

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock)
{
	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nActiveBlock)
	{
		// retrieve actual block index from the active list
		uint block_idx = list[list_idx];

		float F;
		bool isValid;
		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
		uint base_addr = block_idx*blocksize;

		uint xgridlength = xdim/BLOCK_LENGTH;
		uint ygridlength = ydim/BLOCK_LENGTH;
		uint zgridlength = zdim/BLOCK_LENGTH;

		// compute block index
		uint bx = block_idx%xgridlength;
		uint tmpIdx = (block_idx - bx)/xgridlength;
		uint by = tmpIdx%ygridlength;
		uint bz = (tmpIdx-by)/ygridlength;

		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		__shared__ float _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

		// copy global to shared memory
		dim3 idx(tx+1,ty+1,tz+1);

		SOL(idx.x,idx.y,idx.z) = sol_in[base_addr + tIdx];
		
        F = spd[base_addr + tIdx];
		
		isValid = mask[base_addr + tIdx];

		uint new_base_addr, new_tIdx;

		// 1-neighborhood values
		if(tx == 0) 
		{
			if(bx == 0) // end of the grid
			{	
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + BLOCK_LENGTH-1;
				new_base_addr = (block_idx - 1)*blocksize;	
			}

			SOL(tx,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
		}

		if(tx == BLOCK_LENGTH-1)
		{
			if(bx == xgridlength-1) // end of the grid
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1);
				new_base_addr = (block_idx + 1)*blocksize;	
			}
			SOL(tx+2,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
		}

		if(ty == 0)
		{
			if(by == 0)
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
				new_base_addr = (block_idx - xgridlength)*blocksize;
			}

			SOL(idx.x,ty,idx.z) = sol_in[new_base_addr + new_tIdx];
		}

		if(ty == BLOCK_LENGTH-1)
		{
			if(by == ygridlength-1) 
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
				new_base_addr = (block_idx + xgridlength)*blocksize;
			}

			SOL(idx.x,ty+2,idx.z) = sol_in[new_base_addr + new_tIdx];
		}

		if(tz == 0)
		{
			if(bz == 0)
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
				new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
			}

			SOL(idx.x,idx.y,tz) = sol_in[new_base_addr + new_tIdx];
		}

		if(tz == BLOCK_LENGTH-1)
		{
			if(bz == zgridlength-1) 
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
				new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
			}

			SOL(idx.x,idx.y,tz+2) = sol_in[new_base_addr + new_tIdx];
		}

		__syncthreads();

		float a,b,c,oldT,newT;

		for(int iter=0; iter<nIter; iter++)	
		{
			// compute new value
			oldT = newT = SOL(idx.x,idx.y,idx.z);

			if(isValid)
			{
				a = min(SOL(tx,idx.y,idx.z),SOL(tx+2,idx.y,idx.z));
				b = min(SOL(idx.x,ty,idx.z),SOL(idx.x,ty+2,idx.z));
				c = min(SOL(idx.x,idx.y,tz),SOL(idx.x,idx.y,tz+2));

				float tmp = (float) get_time_eikonal(a, b, c, dh, F);

				newT = min(tmp,oldT);
			}
			__syncthreads();	

			if(isValid) SOL(idx.x,idx.y,idx.z) = newT;
		}

		float residue = oldT - newT;

		// write back to global memory
		con[base_addr + tIdx] = (residue < EPS) ? true : false;
		sol_out[base_addr + tIdx] = newT;		
	}
}

__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock)
{
	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nActiveBlock)
	{
		uint block_idx = list[list_idx];

		__shared__ bool conv[BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH];

		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH/2;
		uint base_addr = block_idx*blocksize*2;
		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		conv[tIdx] = con[base_addr + tIdx];
		conv[tIdx + blocksize] = con[base_addr + tIdx + blocksize];

		__syncthreads();

		for(uint i=blocksize; i>0; i/=2)
		{
			if(tIdx < i)
			{
				bool b1, b2;
				b1 = conv[tIdx];
				b2 = conv[tIdx+i];
				conv[tIdx] = (b1 && b2) ? true : false ;
			}
			__syncthreads();
		}

        // active list is negation of tile convergence (active = not converged)

		if(tIdx == 0) listVol[block_idx] = !conv[0]; 
	}
}

__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock)
{
	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nTotalBlock)
	{
		float F;
		bool isValid;
		
        __shared__ float _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

		uint block_idx = list[list_idx];
		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
		uint base_addr = block_idx*blocksize;

		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		if(list_idx < nActiveBlock) // copy value
		{
			sol_out[base_addr + tIdx] = sol_in[base_addr + tIdx];
		} 
		else
		{
			uint xgridlength = xdim/BLOCK_LENGTH;
			uint ygridlength = ydim/BLOCK_LENGTH;
			uint zgridlength = zdim/BLOCK_LENGTH;

			// compute block index
			uint bx = block_idx%xgridlength;
			uint tmpIdx = (block_idx - bx)/xgridlength;
			uint by = tmpIdx%ygridlength;
			uint bz = (tmpIdx-by)/ygridlength;

			// copy global to shared memory
			dim3 idx(tx+1,ty+1,tz+1);
			
            _sol[idx.x][idx.y][idx.z] = sol_in[base_addr + tIdx];
			
            F = spd[base_addr + tIdx];
			
            if(F > 0) F = 1.0/F;
			
            isValid = mask[base_addr + tIdx];

			uint new_base_addr, new_tIdx;

			// 1-neighborhood values
			if(tx == 0) 
			{
				if(bx == 0) // end of the grid
				{	
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + BLOCK_LENGTH-1;
					new_base_addr = (block_idx - 1)*blocksize;	
				}
				_sol[tx][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
			}

			if(tx == BLOCK_LENGTH-1)
			{
				if(bx == xgridlength-1) // end of the grid
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1);
					new_base_addr = (block_idx + 1)*blocksize;	
				}
				_sol[tx+2][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
			}

			if(ty == 0)
			{
				if(by == 0)
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
					new_base_addr = (block_idx - xgridlength)*blocksize;
				}
				_sol[idx.x][ty][idx.z] = sol_in[new_base_addr + new_tIdx];
			}

			if(ty == BLOCK_LENGTH-1) 
			{
				if(by == ygridlength-1) 
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
					new_base_addr = (block_idx + xgridlength)*blocksize;
				}
				_sol[idx.x][ty+2][idx.z] = sol_in[new_base_addr + new_tIdx];
			}

			if(tz == 0)
			{
				if(bz == 0)
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
					new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
				}
				_sol[idx.x][idx.y][tz] = sol_in[new_base_addr + new_tIdx];
			}

			if(tz == BLOCK_LENGTH-1)
			{
				if(bz == zgridlength-1) // end of the grid
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
					new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
				}
				_sol[idx.x][idx.y][tz+2] = sol_in[new_base_addr + new_tIdx];
			}

			__syncthreads();


			float a, b, c, oldT, newT;

			// compute new value
			oldT = newT = _sol[idx.x][idx.y][idx.z];

			if(isValid)
			{
				a = min(_sol[tx][idx.y][idx.z],_sol[tx+2][idx.y][idx.z]);
				b = min(_sol[idx.x][ty][idx.z],_sol[idx.x][ty+2][idx.z]);
				c = min(_sol[idx.x][idx.y][tz],_sol[idx.x][idx.y][tz+2]);

				float tmp = (float) get_time_eikonal(a, b, c, dh, F);
				newT = min(tmp,oldT);

				sol_out[base_addr + tIdx] = newT;
			}

			// write back to global memory
			float residue = oldT - newT;
			con[base_addr + tIdx] = (residue < EPS) ? true : false;	
		}
	}
}
