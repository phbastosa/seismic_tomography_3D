# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../../essentials/inout/inout.hpp"
# include "../../essentials/utils/utils.hpp"
# include "../../essentials/model/model.hpp"
# include "../../essentials/geometry/geometry.hpp"

class Eikonal3D
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
    Utils utils;
    Model3D m3D;
    Geometry3D g3D;

    std::string geomPath;
    std::string eikonalPath;             // Folder to write travel times volume 
    std::string arrivalsPath;            // Folder to write first arrivals

    bool reciprocity;
    bool saveGeometry;
    bool exportTimesVolume;              // To set if you want to write the times volume 
    bool exportFirstArrivals;            // To set if you want to write the first arrivals

    std::string parametersFile;

    Eikonal3D(char **argv);

    /* Function to calculate minimum value between two inputs */
    float min(float v1, float v2);
    
    /* Function to calculate minimum value between four inputs */
    float min4(float v1, float v2, float v3, float v4);

    /* */
    void setup();

    /* */
    void deleteVolumes();

    /* */
    void allocateVolumes();

    /* */
    void podvin3D();
    
    /* */ 
    void fim3D();

    /* */
    void writeTravelTimes();
    
    /* */
    void writeFirstArrivals();
};

# endif
