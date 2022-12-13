# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../eikonal/eikonal.hpp"

class Tomography : public Eikonal
{
private:
    
    /* 0 - L2 norm least squares (Berriman regularization)
       1 - L2 norm least squares (zero order Tikonov regularization)  
       2 - L2 norm least squares (first order Tikonov regularization)  
       3 - L2 norm least squares (second order Tikonov regularization) 
    */
    int inversionMethod;       

    int iteration;                  // It counts the current iteration of the inversion
    int maxIteration;               // It defines the amoung of iterations the inversion have in total 

    float lambda;                   // Regularization parameter
    std::vector<std::string> xMask; // Mask to avoid boundary outliers in x direction 
    std::vector<std::string> yMask; // Mask to avoid boundary outliers in y direction
    std::vector<std::string> zMask; // Mask to avoid boundary outliers in z direction

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

    std::vector<float> residuo;     // To store the convergency at each iteration

    /* */
    void gradientRayTracing();

    /* */
    void tomographyUpdate();
    
    /* */
    void lscg_Berriman();

    /* */
    void lscg_zoTikhonov();
    
    /* */
    void lscg_foTikhonov();

    /* */
    void lscg_soTikhonov();

    /* */
    void gradientDescent();

public:    
    
    std::string dobsPath;           // Observed data location with its preamble
    std::string dcalPath;           // Calculated data location with its preamble

    /* Tomograpy class constructor */
    Tomography();

    void setParameters(char * parametersFile);

    /* */
    void infoMessage();

    /* */
    void importDobs();

    /* */
    void importDcal();
    
    /* */
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
    void modelSmoothing();
    
    /* */
    void exportConvergency();
};

# endif
