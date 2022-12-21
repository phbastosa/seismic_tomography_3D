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

    int nt;               //
    int nsrc;             //
    int shotId;           //
    int timeStep;         //

    int sIdx;             // 
    int sIdy;             //
    int sIdz;             //

    float dt;             // 
    float fcut;           //
    float tlag;           //

    float * V;            //
    float * U_pas;        //
    float * U_pre;        //
    float * U_fut;        //
    
    float * source;       //
    
    float factor;         // 
    float * damp1D;       //
    float * damp2D;       //
    float * damp3D;       //

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
