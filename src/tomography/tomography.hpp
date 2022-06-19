# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../eikonal/eikonal.hpp"

class Tomography : public Eikonal
{
private:

    int iteration;                 //
    int maxIteration;              //  
    
    /* 0 - L2 norm least squares (Berriman regularization)
       1 - L2 norm least squares (zero order Tikonov regularization)  
       2 - L2 norm least squares (first order Tikonov regularization)  
       3 - L2 norm least squares (second order Tikonov regularization) 
       4 - comming soon... 
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
    float * deltaSlowness;         //

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
    void tomographyUpdate();
    
    /* */
    void lscg_Berriman();

    /* */
    void lscg_zoTikhonov();
    
    /* */
    void lscg_foTikhonov();

    /* */
    void lscg_soTikhonov();

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
