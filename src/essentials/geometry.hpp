# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include <string>
# include <vector>

class Geometry
{
private:

    typedef struct                   // Position struct to storage coordinates    
    {
        float * x;                   // Position x of geometry element
        float * y;                   // Position y of geometry element
        float * z;                   // Position z of geometry element

        int idx;                     // Position x in samples with boundary compansation
        int idy;                     // Position y in samples with boundary compansation  
        int idz;                     // Position z in samples with boundary compansation 

        int all;                     // Total elements in simulation
        int n_xline;                 // Total geometry unit in xline
        int n_yline;                 // Total geometry unit in yline

        float xcenter;               // X center circle position
        float ycenter;               // Y center circle position   
        float elevation;             // Geometry element elevation (always positive number)
        float circle_spacing;        // Distance between geometry element in circle 

        std::vector<float> offsets;  // Circle ratios

    } Position;

    typedef struct                   // Struct to storage a simple 2D point
    {
        float x;                     // x plane position
        float y;                     // y plane position

    } Point;            

    Point SW;                        // South western point
    Point NW;                        // North western point
    Point SE;                        // South eastern point

public:  

    Position shots;                  // Struct to storage x,y,z coordinates of shots
    Position nodes;                  // Struct to storage x,y,z coordinates of nodes

    std::string shotsPath;           // Location to export shots position in .txt file
    std::string nodesPath;           // Location to export nodes position in .txt file

    int shotsGeometryType;           // 0 - circular shots | 1 - grid shots
    int nodesGeometryType;           // 0 - circular nodes | 1 - grid nodes

    bool shotsTopography;            // To set if shots have specific topography
    bool nodesTopography;            // To set if nodes have specific topography

    std::string shotsTopographyPath; // Binary file of shots z position
    std::string nodesTopographyPath; // Binary file of nodes z position

    bool reciprocity;                // To set reciprocity 
    bool saveGeometry;               // To save geometry

    /* Function to calculate a linear spaced position */
    std::vector<float> linspace(float xi, float xf, int n);        

    /* South-western element in geometry */
    void set_SW(float x, float y);

    /* North-western element in geometry */
    void set_NW(float x, float y);
    
    /* South-estern element in geometry */
    void set_SE(float x, float y);

    /* Method to set shots positions giving three points in grid */
    void setGridGeometry(Position &obj);
    
    /* Method to set nodes positions giving three points in grid */
    void setCircularGeometry(Position &obj);

    /* Switch the nodes position to shots position */
    void setReciprocity();

    /* Method to save geometry in hard drive */
    void exportPositions();
};

# endif
