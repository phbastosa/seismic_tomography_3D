# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

class Eikonal
{   
private:

public:    

    int shotId;                  // Current source index available
    int eikonalType;             //
    int shotsGeometryType;
    int nodesGeometryType;

    float * T;                   // Travel times volume
    float * S;                   // Slowness volume
    float * K;                   // Wavefront expansion volume
    float * nT;                  // Auxiliar travel times volume
    float * nK;                  // Auxiliar wavefromt expansion volume

    InOut io;
    Model m3D;
    Utils utils;
    Geometry g3D;

    std::string geomPath;
    std::string eikonalPath;             // Folder to write travel times volume 
    std::string arrivalsPath;            // Folder to write first arrivals
    std::string parametersFile;

    bool reciprocity;
    bool saveGeometry;
    bool exportTimesVolume;              // To set if you want to write the times volume 
    bool exportFirstArrivals;            // To set if you want to write the first arrivals

    /* Function to calculate minimum value between two inputs */
    float min(float v1, float v2);

    /* */
    void deleteVolumes();

    /* */
    void allocateVolumes();

    /* */
    void podvin();
    
    /* */ 
    void jeongFIM();


    /* */
    void nobleFSM();
    void initSweep();
    void fullSweep();
    void innerSweep(int i, int j, int k, int sx, int sy, int sz, int sgntz, int sgntx, int sgnty, int sgnvz, int sgnvx, int sgnvy, float dzi, float dxi, float dyi, float dz2i, float dx2i, float dy2i, float dz2dx2, float dz2dy2, float dx2dy2, float dsum);

    /* */
    void writeTravelTimes();
    
    /* */
    void writeFirstArrivals();
};

# endif
