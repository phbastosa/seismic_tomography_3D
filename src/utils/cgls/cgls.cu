# include <cmath>
# include <iostream>

# include <cuda.h>
# include <cublas_v2.h>
# include <cusparse.h>

# include "cgls.cuh"

void sparse_cgls_gpu(int * iA, int * jA, float * vA, float * B, float * x, int N, int M, int NNZ, int NIT, float TOL)
{
    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&cublasHandle);
    if (cublasStatus != CUBLAS_STATUS_SUCCESS) printf("Error cublas\n");

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);
    if (cusparseStatus != CUSPARSE_STATUS_SUCCESS) printf("Error cusparse\n");

    size_t bsize;
    void * buffer;
    float beta = 0.0f;
    float alpha = 1.0f;
    float a, b, qTq, rTr, rd;

    int * d_iA_coo; cudaMalloc((void **)&d_iA_coo, NNZ * sizeof(int));
    int * d_iA_csr; cudaMalloc((void **)&d_iA_csr,(N+1)* sizeof(int));

    cudaMemcpy(d_iA_coo, iA, NNZ * sizeof(int), cudaMemcpyHostToDevice);

    cusparseXcoo2csr(cusparseHandle, d_iA_coo, NNZ, N, d_iA_csr, CUSPARSE_INDEX_BASE_ZERO);

    cudaFree(d_iA_coo);
		
    float * d_p; cudaMalloc((void **)&d_p, M * sizeof(float)); 
    float * d_q; cudaMalloc((void **)&d_q, N * sizeof(float));  
    float * d_r; cudaMalloc((void **)&d_r, M * sizeof(float)); 
    float * d_s; cudaMalloc((void **)&d_s, N * sizeof(float)); 
    float * d_x; cudaMalloc((void **)&d_x, M * sizeof(float)); 

    float * d_vA; cudaMalloc((void **)&d_vA, NNZ * sizeof(float));     
    int * d_jA_coo; cudaMalloc((void **)&d_jA_coo, NNZ * sizeof(int)); 

    cudaMemset(d_x, 0, M * sizeof(float));    
    cudaMemset(d_p, 0, M * sizeof(float));
    cudaMemset(d_q, 0, N * sizeof(float));
    cudaMemset(d_r, 0, M * sizeof(float));
    cudaMemcpy(d_s, B, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_vA, vA, NNZ * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_jA_coo, jA, NNZ * sizeof(int), cudaMemcpyHostToDevice);

    cusparseDnVecDescr_t Dn_p;    
    cusparseDnVecDescr_t Dn_q;    
    cusparseDnVecDescr_t Dn_r;    
    cusparseDnVecDescr_t Dn_s;    
    cusparseSpMatDescr_t Sp_matA; 

    cusparseCreateCsr(&Sp_matA, N, M, NNZ, d_iA_csr, d_jA_coo, d_vA, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_p, M, d_p, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_q, N, d_q, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_r, M, d_r, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_s, N, d_s, CUDA_R_32F);

    cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, &bsize);
    cudaMalloc(&buffer, bsize);

    cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
    cudaDeviceSynchronize();    

    cublasScopy_v2(cublasHandle, M, d_r, 1, d_p, 1);

    cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, Sp_matA, Dn_p, &beta, Dn_q, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
    cudaDeviceSynchronize();    

    for (int iteration = 0; iteration < NIT; iteration++)
    {
        qTq = 0.0f;
        cublasSdot_v2(cublasHandle, N, d_q, 1, d_q, 1, &qTq);
        cudaDeviceSynchronize();    // qTq = q' * q

        rTr = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rTr);
        cudaDeviceSynchronize();    // rTr = r' * r 

        a = rTr / qTq;              // a = (r' * r) / (q' * q)
        cublasSaxpy_v2(cublasHandle, M, &a, d_p, 1, d_x, 1);
        cudaDeviceSynchronize();    // x = x + a * p

        a *= -1.0f;
        cublasSaxpy_v2(cublasHandle, N, &a, d_q, 1, d_s, 1);
        cudaDeviceSynchronize();    // s = s - a * q 

        rd = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rd);
        cudaDeviceSynchronize();    // rd = r' * r

        if (sqrtf(rd) < TOL) break; // Convergence condition 

        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
        cudaDeviceSynchronize();   // r = G' * s    

        rTr = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rTr);
        cudaDeviceSynchronize();   // rTr = r' * r 

        b = rTr / rd;              // b = (r' * r) / rd  
        cublasSscal_v2(cublasHandle, M, &b, d_p, 1);
        cudaDeviceSynchronize();   // p = b * p  

        b = 1.0f;
        cublasSaxpy_v2(cublasHandle, M, &b, d_r, 1, d_p, 1);
        cudaDeviceSynchronize();   // p += r  <---> p = r + b * p  
        
        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, Sp_matA, Dn_p, &beta, Dn_q, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
        cudaDeviceSynchronize();   // q = G * p    
    }

    cudaMemcpy(x, d_x, M * sizeof(float), cudaMemcpyDeviceToHost);

    cusparseDestroyDnVec(Dn_p);
    cusparseDestroyDnVec(Dn_q);
    cusparseDestroyDnVec(Dn_r);
    cusparseDestroyDnVec(Dn_s);
    cusparseDestroySpMat(Sp_matA);

    cudaFree(d_vA);
    cudaFree(d_iA_csr);
    cudaFree(d_jA_coo);

    cudaFree(d_x);
    cudaFree(d_p);
    cudaFree(d_q);
    cudaFree(d_r);
    cudaFree(d_s);

    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);
}

