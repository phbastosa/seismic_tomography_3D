# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../eikonal/eikonal.hpp"

class Tomography : public Eikonal
{
private:
    
    int iteration;                  // It counts the current iteration of the inversion
    int maxIteration;               // It defines the amoung of iterations the inversion have in total 

    int tkOrder;                    // Tikhonov regularization order parameter 
    float lambda;                   // Regularization parameter

    bool smooth;                    // Parameter to select model smoothing per iteration
    int smoothingType;              // 0 - gaussian filter; 1 - moving average filter  
    int filterSamples;              // The amoung of samples in filter
    float standardDeviation;        // Standard deviation for gaussian filter

    float * dm;                     // Delta slowness model, used to realize inversion
    float * dobs;                   // Observed data 
    float * dcal;                   // Calculated data
    float * model;                  // Slowness model
    
    typedef struct                  // To get together the reducecd tomography model dimensions
    {  
        int nx;                     // Reduced samples in x dimension
        int ny;                     // Reduced samples in y dimension
        int nz;                     // Reduced samples in z dimension
        int nPoints;                // Total points in reduced dimension             
        float dx;                   // Sparse sample spacing in x dimension
        float dy;                   // Sparse sample spacing in y dimension
        float dz;                   // Sparse sample spacing in z dimension    

    } tomoModel;               
    
    tomoModel mTomo;                // Struct to organize reduced model parameter   

    std::string residuoPath;        // Folder and file to write the inversion convergency 
    std::string estimatedPath;      // Folder to write the estimated model 

    std::vector< int > iM;          // To store initially the rows of sparse G matrix 
    std::vector< int > jM;          // To store initially the cols of sparse G matrix
    std::vector<float> vM;          // To store initially the values of sparse G matrix

    float * B;                      // Data difference (dobs - dcal) with regularization     
    sparseMatrix A;                 // Regularized matrix to solve linear system (A dm = dT)  

    std::vector<float> residuo;     // To store the convergency at each iteration

    std::vector<std::string> xMask; // Mask to avoid boundary outliers in x direction 
    std::vector<std::string> yMask; // Mask to avoid boundary outliers in y direction
    std::vector<std::string> zMask; // Mask to avoid boundary outliers in z direction

    /*  
        Ray tracing using travel times as a function
        Stepest descent gradient used to compute the steps of ray path
        Indexes organized to store iM,jM, and vM to build inversion matrix
    */
    void gradientRayTracing();

    /* 
        Construct the A matrix (G + L) in form of sparse matrix deleting iM, jM and vM 
    */
    void buildRegularizedMatrix();   

    /*
    
    */
    void buildRegularizedData();

public:    

    std::string dobsPath;           // Observed data location with its preamble
    std::string parameters;         // 

    /* Tomograpy class constructor */
    Tomography();

    /* Importing from file and setting all tomography parameters */
    void setParameters();

    /* 
        It informs:
        What is the current iteration 
        The shot position 
        The previous resuduo at each iteration 
    */
    void infoMessage();

    /* 
        To import observed data 
        Each shot in different binary file 
    */
    void importDobs();

    /* 
        To import calculated data 
        Each shot in different binary file 
    */
    void importDcal();
    
    /* To reduce velocity model to update inversion */
    void setInitialModel();

    /* */
    void forwardModeling();

    /* */
    bool converged();

    /* */
    void optimization();

    /* */
    void modelUpdate();    
    
    /* */
    void exportConvergency();
};

# endif
