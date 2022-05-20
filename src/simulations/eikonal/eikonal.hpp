# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../../essentials/model/model.hpp"
# include "../../essentials/geometry/geometry.hpp"

class Eikonal3D
{   
private:

public:    

    int shotId;                  // Current source index available
    
    float * T;                   // Travel times volume
    float * S;                   // 
    float * K;
    float * nT;
    float * nK;

    std::string eikonalPath;     // Folder to write travel times volume 
    std::string arrivalsPath;    // Folder to write first arrivals

    bool exportTimesVolume;      // To set if you want to write the times volume 
    bool exportFirstArrivals;    // To set if you want to write the first arrivals
    
    /* Function to calculate minimum value between two inputs */
    float min(float v1, float v2);
    
    /* Function to calculate minimum value between four inputs */
    float min4(float v1, float v2, float v3, float v4);

    /* */
    void deleteVolumes();

    /* */
    void allocateVolumes(Model3D m3D);

    /* */
    void podvin3D(Model3D m3D, Geometry3D geom);
    
    /* */
    void writeTravelTimes(Model3D m3D);
    
    /* */
    void writeFirstArrivals(Model3D m3D, Geometry3D geom);
};

# endif
