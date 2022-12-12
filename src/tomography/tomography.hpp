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

    int iteration;                  //
    int maxIteration;               //  

    float lambda;                   //
    std::vector<std::string> xMask; // 
    std::vector<std::string> yMask; // 
    std::vector<std::string> zMask; //

    bool smooth;                    //
    int smoothingType;              // 
    int filterSamples;              // 
    float standardDeviation;        //

    float * dm;                     //
    float * dobs;                   //
    float * dcal;                   //
    float * model;                  //

    typedef struct                  //
    {  
        int nx, ny, nz, nPoints;    //             
        float dx, dy, dz;           //    

    } tomoModel;               
    
    tomoModel mTomo;                //    

    std::string dobsPath;           //
    std::string dcalPath;           //
    std::string residuoPath;        //
    std::string estimatedPath;      //

    std::vector< int > iM;          // 
    std::vector< int > jM;          //
    std::vector<float> vM;          //

    std::vector<float> residuo;     //

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
    
    /* */
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
