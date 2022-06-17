# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../eikonal/eikonal.hpp"

class Tomography : public Eikonal
{
private:

    int iteration;                 //
    int maxIteration;              //  
    
    /* 0 - Least squares (Berriman regularization) 
       1 - ...
       2 - ...
       3 - comming soon... 
    */
    int optimizationMethod;       

    float xMask;                   // 
    float yMask;                   // 
    float lambda;                  //
    float zMaskUp;                 //
    float zMaskDown;               //
    float tomoTolerance;           //

    float * dobs;                  //
    float * dcal;                  //
    float * gradient;              //
    float * slowness;              //

    bool smoothing;                //
    bool generate_dobs;            //

    typedef struct                 //
    {  
        int nx, ny, nz, nPoints;   //             
        float dx, dy, dz;          //    

    } tomoModel;               
    
    tomoModel mTomo;               //    

    std::string resPath;           //
    std::string dobsPath;          //
    std::string dcalPath;          //
    std::string gradPath;          //
    std::string estModels;         //

    std::vector < int > iM;        // 
    std::vector < int > jM;        //
    std::vector <float> vM;        //

    std::vector <float> residuo;   //

    /* */
    void generateDobs();

    /* */
    void gradientRayTracing();

    /* */
    void makeGradient();

    /* */
    void cgls_Berriman();

    /* */
    void cgls_zoTikhonov();
    
    /* */
    void cgls_foTikhonov();

    /* */
    void cgls_soTikhonov();

public:    
    
    /* */
    Tomography(char **argv);

    /* */
    void infoMessage();

    /* */
    void importDobs();

    /* */
    void setInitialModel();

    /* */
    void importDcal();
    
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
