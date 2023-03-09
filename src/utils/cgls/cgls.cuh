# ifndef CGLS_CUH
# define CGLS_CUH

/* 
    Function to compute a sparse matrix least square conjugate gradient 
    
    Solution of A'A x = A'B without generate A'A 

    iA - rows in coo format 
    iA - columns in coo format
    vA - values of sparse matrix A 
    
    B - second member of linear system  
    
    x - solution

    N - number of rows 
    M - number of columns
    NNZ - Non zero elements
    NIT - max iterations 
    TOL - tolerance 
    
    Solving problem using CPU only 
    Memory needed = 4 * (3*M + 2*N + 3*NNZ) / 1024 / 1024 MB    
*/
void sparse_cgls_cpu(int * iA, int * jA, float * vA, float * B, float * x, int N, int M, int NNZ, int NIT, float TOL);

/* 
    Function to compute a sparse matrix least square conjugate gradient 
    
    Solution of A'A x = A'B without generate A'A 

    iA - rows in coo format 
    iA - columns in coo format
    vA - values of sparse matrix A 
    
    B - second member of linear system  
    
    x - solution

    N - number of rows 
    M - number of columns
    NNZ - Non zero elements
    NIT - max iterations 
    TOL - tolerance 

    Solving problem using CUDA, cublas and cusparse on GPU
    Memmory needed = 4 * (3*M + 2*N + 2*NNZ + (N+1)) / 1024 / 1024 MB
*/
void sparse_cgls_gpu(int * iA, int * jA, float * vA, float * B, float * x, int N, int M, int NNZ, int NIT, float TOL);

# endif  

