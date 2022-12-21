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
    float dt;             // 
    float fcut;           //
    float tlag;           //

    float * U_pre;        //
    float * U_pas;        //

    float * source;       //
    
    float factor;         // 
    float * damp1D;       //
    float * damp2D;       //
    float * damp3D;       //

    /* */
    void sourceGenerator();
    
    
    void dampingGenerator();
    void acousticWaveSolver();
    void dampingApplicator();
    void wavefieldUpdate();
    void getSeismogram();
    void exportSeismogram();

};

# endif
