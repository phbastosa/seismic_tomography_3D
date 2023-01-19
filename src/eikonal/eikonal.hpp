# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <string>

# include "../essentials/utils.hpp"
# include "../essentials/model.hpp"
# include "../essentials/geometry.hpp"

class Eikonal : public Utils, public Model, public Geometry
{   
private:

    typedef struct
    {
        int sgntz; int sgntx;            // 
        int sgnty; int sgnvy;            //
        int sgnvz; int sgnvx;            //

        int i, j, k;                     //

        float dzi, dxi, dyi;             //
        float dz2i, dx2i, dy2i;          //
        float dz2dx2, dz2dy2;            //
        float dx2dy2, dsum;              //
        
    } FSM;
    
    FSM fsm;                             // Compressing fast sweeping method variables

    /* */
    float min(float v1, float v2);

    /* */
    void writeTravelTimes();
    
    /* */
    void writeFirstArrivals();

    /* */
    void initSweep();
    
    /* */
    void fullSweep();

    /* */
    void innerSweep();

    /* */
    void podvin();
    
    /* */ 
    void jeongFIM();

    /* */
    void nobleFSM();

    /* */
    void rayTracing();

public:    
    
    /* 0 - podvin & Lecomte (1991); 
       1 - Jeong & Whitaker (2008); 
       2 - Noble, Gesret & Belayouni (2014); 
    */  
    int eikonalType;

    int shotId;                          // Current source index available

    float * V;                           // Velocity volume 
    float * T;                           // Travel times volume
    float * S;                           // Slowness volume 

    float * illumination;                // Illumination matrix for all shots 

    bool exportTimesVolume;              // To set if you want to write the times volume 
    bool exportFirstArrivals;            // To set if you want to write the first arrivals
    bool exportIllumination;             // To set if you want to write illumination matrix
    bool exportRayPosition;              // To set if you want to write ray positions 

    std::string raysFolder;              // Folder to write ray position points 
    std::string eikonalFolder;           // Folder to write travel times volume 
    std::string arrivalFolder;           // Folder to write first arrivals
    std::string illuminationFolder;      // Folder to write illumination volume

    /* */
    void eikonalComputing();

    /* */
    void writeIllumination();
    
    /* */
    void setEikonalParameters();
};

# endif
