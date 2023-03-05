# ifndef CLASSIC_HPP
# define CLASSIC_HPP

# include "../eikonal.hpp"

class Classic : public Eikonal
{
private:

    float * S;
    float * T;
    float * K;
    float * nT;
    float * nK;

    void expand_model();
    void reduce_model();

public:

    void set_parameters();
    void prepare_volumes();

    void info_message();
    void eikonal_equation();
    void write_time_volume();
    void write_first_arrival();

    void destroy_volumes();    
};

# endif