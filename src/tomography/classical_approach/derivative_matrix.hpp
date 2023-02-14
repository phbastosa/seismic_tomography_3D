# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>
# include <string>

class Utils
{
public:
    /* 
    Struct to define a sparse matrix:
       n - number of rows
       m - number of cols
       nnz - number of non zero elements
       i - row index
       j - col index
       v - values 
    */ 
    typedef struct     
    {
        int * i;       // Rows indexes
        int * j;       // Cols indexes 
        float * v;     // Value

        int n;         // Rows number
        int m;         // Cols number
        int nnz;       // Non-zero elements
    
        /* Function to initialize the vectorial components of sparse matrix */
        void init()
        {
            i = new int[nnz]();  
            j = new int[nnz](); 
            v = new float[nnz](); 
        }

        /* Function to delete the vectorial components of sparse matrix */
        void erase()
        {
            delete[] i;
            delete[] j;
            delete[] v;
        }

    } sparseMatrix;

    std::string parameters;

    /* */
    sparseMatrix getDerivativeMatrix(int n, int degree);

};

# endif
