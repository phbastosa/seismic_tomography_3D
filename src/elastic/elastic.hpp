# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include <string>

# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

class Elastic : public Utils, public Model, public Geometry
{
private:

    float * M;                           //
    float * L;                           //                
    
    float * Vx;                          //
    float * Vy;                          // 
    float * Vz;                          //  
    float * Txx;                         //
    float * Tyy;                         //
    float * Tzz;                         //         
    float * Txz;                         //
    float * Tyz;                         //
    float * Txy;                         //

public:

    int nt;
    float dt;                            // 
    float * cerjan;                      //
    
    int nsrc;                            // 
    float * source;                      // 

    int shotId;                          // Current source index available
    int timeId;                          //

    float * seismogram;                  //
    bool exportSeismograms;              //  
    std::string seismFolder;             // 

    /* */
    void deleteVolumes();
    
    /* */
    void allocateVolumes();

    /* */
    void setWavefields();
    
    /* */
    void forwardModeling();

    /* */
    void elasticIsotropic_FD8E2T();

    /* */
    void cerjanAbsorbingCondition();

    /* */
    void getPressureSeismogram();
};

# endif
