# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <string>

# include "../essentials/utils/utils.hpp"
# include "../essentials/model/model.hpp"
# include "../essentials/geometry/geometry.hpp"

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
    
    float * S;                           // Slowness volume
    float * K;                           // Wavefront expansion volume
    float * nT;                          // Auxiliar travel times volume
    float * nK;                          // Auxiliar wavefromt expansion volume

    FSM fsm;                             // Compressing fast sweeping method variables

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

public:    
    
    /* 0 - podvin & Lecomte (1991); 
       1 - Jeong & Whitaker (2008); 
       2 - Noble, Gesret & Belayouni (2014); 
    */  
    int eikonalType;

    int shotId;                          // Current source index available

    float * T;                           // Travel times volume

    bool exportTimesVolume;              // To set if you want to write the times volume 
    bool exportFirstArrivals;            // To set if you want to write the first arrivals

    std::string eikonalFolder;           // Folder to write travel times volume 
    std::string arrivalFolder;           // Folder to write first arrivals

    /* */
    void eikonalComputing();
};

# endif
