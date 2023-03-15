# ifndef LEAST_SQUARES_HPP
# define LEAST_SQUARES_HPP

# include "../tomography.hpp"

class Least_squares : public Tomography
{
private:
    
    int tk_order;   
    float lambda; 

    float * T;
    float * illumination;

    float dx_tomo;
    float dy_tomo;
    float dz_tomo;

    int nx_tomo;  
    int ny_tomo;  
    int nz_tomo;  

    std::vector< int > iG; 
    std::vector< int > jG;
    std::vector<float> vG;

    int * iA;
    int * jA;
    float * vA;
    float * B;
    float * xdm;             // A x = B

    int M, N, NNZ;

    void info_message();
    void ray_tracing();

    void expand_fdm();
    void tk_reg_matrix();
    void sparse_cgls();

public:    

    void set_parameters();
    void import_obs_data();
    
    void forward_modeling();
    void import_cal_data();

    bool converged();

    void optimization();
    void model_update();

    void export_convergency();
};

# endif
