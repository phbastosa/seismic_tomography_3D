# ifndef LEAST_SQUARES_HPP
# define LEAST_SQUARES_HPP

# include "../tomography.hpp"

class Least_squares : public Tomography
{
private:
    
    int tkOrder;   
    float lambda; 

    float * x;             // A x = B
    float * illumination;

    float dx_tomo;
    float dy_tomo;
    float dz_tomo;

    float nx_tomo;  
    float ny_tomo;  
    float nz_tomo;  

    std::vector< int > iG; 
    std::vector< int > jG;
    std::vector<float> vG;

    void info_message();
    void ray_tracing();
    void sparse_cgls();
    void resizing();


public:    

    void set_parameters();
    void import_obs_data();
    
    void forward_modeling();
    void import_cal_data();

    bool converged();

    void optimization();
    void model_update();

    void export_convergency();
    void export_estimated_model();    
};

# endif
