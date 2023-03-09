# ifndef LEAST_SQUARES_HPP
# define LEAST_SQUARES_HPP

# include "../tomography.hpp"

class Least_squares : public Tomography
{
private:
    
    int tkOrder;   
    float lambda; 

    std::vector< int > iM;          // To store initially the rows of sparse G matrix 
    std::vector< int > jM;          // To store initially the cols of sparse G matrix
    std::vector<float> vM;          // To store initially the values of sparse G matrix

    void info_message();

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
