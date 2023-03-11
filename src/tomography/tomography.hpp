# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include <string>
# include <vector>

# include "../eikonal/eikonal.hpp"

class Tomography 
{
protected:
    
    int iteration;                  
    int max_iteration;             

    int window;        
    float stdv;        
    bool smooth;       

    float * dm;      
    float * dobs;     
    float * dcal;    

    int n_data;
    int n_model;

    std::string obs_data_folder;
    std::string cal_data_folder;
    std::string convergency_folder; 
    std::string estimated_model_folder; 

    std::vector<float> residuo; 

    Eikonal * eikonal;

    virtual void info_message() = 0;

public:    

    std::string parameters;

    virtual void set_parameters() = 0;
    virtual void import_obs_data() = 0;
    
    virtual void forward_modeling() = 0;
    virtual void import_cal_data() = 0;

    virtual bool converged() = 0;

    virtual void optimization() = 0;
    virtual void model_update() = 0;

    virtual void export_convergency() = 0;
    virtual void export_estimated_model() = 0;    
};

# endif
