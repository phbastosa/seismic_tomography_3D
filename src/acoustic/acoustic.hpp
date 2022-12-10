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

    int ns;
    int nt;
    float dt;
    float fcut;
    float * source;

    float * U_pre;
    float * U_pas;

    float * factor;
    float * prismDampX;
    float * prismDampY;
    float * prismDampZ;
    float * cubeDamper;

    void sourceGenerator();
    void dampingGenerator();
    void acousticWaveSolver();
    void dampingApplicator();
    void wavefieldUpdate();
    void getSeismogram();
    void exportSeismogram();

};

# endif