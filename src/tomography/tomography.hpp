# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../simulations/eikonal/eikonal.hpp"

class Tomography3D : public Eikonal3D
{
public:    
    
    bool generate_dobs;        

    float * dobs;
    float * dcal;
    float * gradient;

    std::string dobsPath;
    std::string dcalPath;
    std::string gradPath;

    std::vector < int > iM;
    std::vector < int > jM;
    std::vector <float> vM;

    int iteration;
    int maxIteration;
    float tomoTolerance;

    Tomography3D(char **argv);

    /* */
    void importDobs();

    /* */
    void importDcal();
    
    /* */
    void generateDobs();

    /* */
    void fwdModeling();

    /* */
    void gradientRayTracing();

    /* */
    void makeGradient();

    /* */
    bool converged();

    /* */
    void cgls_berriman_reg();

    /* */
    void cgls_tikhonov_reg();

    /* */
    


};

# endif
