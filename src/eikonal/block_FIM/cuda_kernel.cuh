# ifndef CUDA_KERNELS_CUH
# define CUDA_KERNELS_CUH

# include <cuda.h>
# include <cuda_runtime.h>

# define INF 1e6f
# define EPS 1e-6f
# define BLOCK_LENGTH 4

# define MEM(index) _mem[index]
# define SOL(i,j,k) _sol[i][j][k]
# define SPD(i,j,k) _spd[i][j][k]

typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct 
{
    uint blknum;
    uint volsize;
    uint blksize;
    uint nActiveBlock;

    int xdim;
    int ydim;
    int zdim;
    int nIter;
    int blklength;
    
    float delta_h;

    uint * h_list;
    bool * h_listed;
    bool * h_listVol;

    float * h_sol;
    float * h_spd;
    bool * h_mask;

} CUDAMEMSTRUCT;

void cuda_safe_call(cudaError_t error);
void block_FIM_solver(CUDAMEMSTRUCT &cmem);

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock);
__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock);
__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock);
__device__ float get_time_eikonal(float a, float b, float c, float h, float s);

# endif