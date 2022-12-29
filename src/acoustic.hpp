# ifndef ACOUSTIC_HPP
# define ACOUSTIC_HPP

# include <string>

# include "../essentials/utils.hpp"
# include "../essentials/model.hpp"
# include "../essentials/geometry.hpp"

class Acoustic : public Utils, public Model, public Geometry
{
private:
    

public:
 
    int nt;                //
    int nsrc;              //
    int shotId;            //
    int timeStep;          //

    float dt;              // 
    float fmax;            //
    float tlag;            //

    float * V;             //
    float * U_pas;         //
    float * U_pre;         //
    float * U_fut;         //
    
    float * source;        //
    
    float factor;          // 
    float * damp1D;        //
    float * damp2D;        //
    float * damp3D;        //

    float * seismogram;    //

    std::string seisLabel; // 
    std::string paramFile; //

    /* */
    void setParameters();

    /* */
    void sourceGenerator();
    
    /* */
    void dampingGenerator();
    
    /* */
    void setWaveField();

    /* */
    void progressMessage();

    /* */
    void forwardModeling();

    /* */
    void applyWavelet();

    /* */
    void wavePropagation();

    /* */
    void dampApplication();
    
    /* */
    void wavefieldUpdate();
    
    /* */
    void buildSeismogram();
    
    /* */
    void exportSeismogram();
};

# endif
