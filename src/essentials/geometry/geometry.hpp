# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

class Geometry
{
private:

    typedef struct                   // Position struct to storage coordinates    
    {
        float * x;                   // Position x of geometry element
        float * y;                   // Position y of geometry element
        float * z;                   // Position z of geometry element
    
        int idx;                     // X index position of geometry element
        int idy;                     // Y index position of geometry element
        int idz;                     // Z index position of geometry element

        int n;                       // Total elements in simulation
        int nx;                      // Total geometry unit in xline
        int ny;                      // Total geometry unit in crossline

        float xc;                    // X center circle position
        float yc;                    // Y center circle position   

        float ds;                    // Distance between geometry element in circle 

        float elevation;             // Geometry element elevation (always positive number)

        std::vector<float> offsets;  // Circle ratios

    } Position;             

public:  

    Position shots;                  // Struct to storage x,y,z coordinates of shots
    Position nodes;                  // Struct to storage x,y,z coordinates of nodes

    Utils::Point SW;               // South western point
    Utils::Point NW;               // North western point
    Utils::Point SE;               // South eastern point

    std::string shotsPath;           // Location to export shots position in .txt file
    std::string nodesPath;           // Location to export nodes position in .txt file

    /* Method to set shots positions giving three points in grid */
    void setGridShots();
    
    /* Method to set nodes positions giving three points in grid */
    void setGridNodes();

    /* Method to set circular shots positions */
    void setCircularShots();

    /* Method to set circular nodes positions */
    void setCircularNodes();

    /* Switch the nodes position to shots position */
    void setReciprocity();

    /* Methood to save geometry in hard drive */
    void exportPositions();
};

# endif
