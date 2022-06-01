# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../essentials/inout/inout.hpp"
# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

class Eikonal
{   
private:
    typedef struct
    {
        int sgntz; int sgntx; int sgnty;
        int sgnvz; int sgnvx; int sgnvy;

        int i, j, k;

        float dzi, dxi, dyi;
        float dz2i, dx2i, dy2i;
        float dz2dx2, dz2dy2, dx2dy2;
        float dsum;  
        
    } FSM;
    
    FSM fsm;
    void initSweep();
    void fullSweep();
    void innerSweep();

    /* Function to calculate minimum value between two inputs */
    float min(float v1, float v2);

public:    

    int shotId;                          // Current source index available
    int eikonalType;                     //
    int shotsGeometryType;               // 
    int nodesGeometryType;               //

    float * T;                           // Travel times volume
    float * S;                           // Slowness volume
    float * K;                           // Wavefront expansion volume
    float * nT;                          // Auxiliar travel times volume
    float * nK;                          // Auxiliar wavefromt expansion volume

    InOut io;                            //
    Model m3D;                           // 
    Utils utils;                         //
    Geometry g3D;                        //

    std::string geomPath;                //
    std::string eikonalPath;             // Folder to write travel times volume 
    std::string arrivalsPath;            // Folder to write first arrivals
    std::string parametersFile;          //

    bool reciprocity;                    // 
    bool saveGeometry;                   //
    bool exportTimesVolume;              // To set if you want to write the times volume 
    bool exportFirstArrivals;            // To set if you want to write the first arrivals

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

    /* */
    void writeTravelTimes();
    
    /* */
    void writeFirstArrivals();
};

# endif
