# include <cmath>
# include <iostream>

# include <cuda.h>
# include <cublas_v2.h>
# include <cusparse.h>

# include "cgls.cuh"

void sparse_cgls_cpu(int * iA, int * jA, float * vA, float * B, float * x, int N, int M, int NNZ, int NIT, float TOL)
{
    float a, b, qTq, rTr, rd;

    float * s = new float[N]();
    float * q = new float[N]();
    float * r = new float[M]();
    float * p = new float[M]();

    // s = d - A * x
    for (int row = 0; row < N; row++) 
        s[row] = B[row]; 

    // r = A' * s    
    for (int ind = 0; ind < NNZ; ind++) 
        r[jA[ind]] += vA[ind] * s[iA[ind]];        

    // p = r and x = 0
    for (int col = 0; col < M; col++) 
    {
        x[col] = 0.0f;
        p[col] = r[col]; 
    }

    // q = A * p
    for (int ind = 0; ind < NNZ; ind++) 
        q[iA[ind]] += vA[ind] * p[jA[ind]];        

    // Iterations loop
    for (int iteration = 0; iteration < NIT; iteration++)
    {
        qTq = 0.0f;
        for (int row = 0; row < N; row++)             // q inner product
            qTq += q[row] * q[row];                   // qTq = q' * q

        rTr = 0.0f;
        for (int col = 0; col < M; col++)             // r inner product
            rTr += r[col] * r[col];                   // rTr = r' * r 

        a = rTr / qTq;                                // a = (r' * r) / (q' * q)                    

        for (int col = 0; col < M; col++)             // model atualization
            x[col] += a * p[col];                     // x = x + a * p

        for (int row = 0; row < N; row++)             // s atualization  
            s[row] -= a * q[row];                     // s = s - a * q 

        rd = 0.0f;
        for (int col = 0; col < M; col++)             // r inner product for division 
            rd += r[col] * r[col];                    // rd = r' * r

        if (sqrtf(rd) < TOL) break;                   // Convergence condition

        for (int col = 0; col < M; col++)             // Zeroing r 
            r[col] = 0.0f;                            // r = 0, for multiplication

        for (int index = 0; index < NNZ; index++)     // r atualization 
            r[jA[index]] += vA[index] * s[iA[index]]; // r = G' * s    
                
        rTr = 0.0f;    
        for (int col = 0; col < M; col++)             // r inner product
            rTr += r[col] * r[col];                   // rTr = r' * r
        
        b = rTr / rd;                                 // b = (r' * r) / rd

        for (int col = 0; col < M; col++)
            p[col] = r[col] + b * p[col];             // p = r + b * p 

        for (int row = 0; row < N; row++) 
            q[row] = 0.0f;                            // q = 0, for multiplication

        for (int index = 0; index < NNZ; index++) 
            q[iA[index]] += vA[index] * p[jA[index]]; // q = G * p           
    }

    delete[] s; delete[] q; delete[] r; delete[] p;
}

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

