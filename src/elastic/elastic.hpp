# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

class Elastic
{
private:

    float * M;
    float * L;
    
    float * Vx;
    float * Vy;
    float * Vz;
    float * Txx;
    float * Tyy;
    float * Tzz;
    float * Txz;
    float * Tyz;
    float * Txy; 

public:

    InOut io;
    Model m3D;
    Utils utils;
    Geometry g3D;

    std::string geomPath;                //
    std::string seismPath;               // 

    int nt;
    int nsrc;
    float dt;
    float * cerjan;
    float * source;
    float * seismogram;

    int shotId;                          // Current source index available
    int timeId;
    bool reciprocity;                    // 
    bool saveGeometry;                   //
    bool exportSeismograms;              //  

    void deleteVolumes();
    void allocateVolumes();

    void setWavefields();
    void forwardModeling();
    void elasticIsotropic_FD8E2T();
    void cerjanAbsorbingCondition();
    void getPressureSeismogram();
};

# endif
