# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/file_manager/file_manager.hpp"

class Geometry
{
protected:

    File_manager fm;
  
    std::string geometry_folder;    
    std::string topography_file;
    
    std::vector<std::string> splitted;
    
    bool topography;
    float elevation;

    void set_topography();
    void get_shot_parameters();
    void get_node_parameters();

    virtual void build_geometry() = 0;

public:  

    float * x;
    float * y;
    float * z;

    int all;
      
    virtual void set_shot_parameters(std::string file) = 0;
    virtual void set_node_parameters(std::string file) = 0;
    
    void import_positions(std::string file);
    void export_positions(std::string file);
};

# endif



