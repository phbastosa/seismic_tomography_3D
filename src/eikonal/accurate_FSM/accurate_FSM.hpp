# ifndef ACCURATE_FSM_HPP
# define ACCURATE_FSM_HPP

# include "../eikonal.hpp"

class Accurate_FSM : public Eikonal
{
private:

    int i, j, k;

    int sgntz, sgntx, sgnty; 
    int sgnvz, sgnvx, sgnvy;

    float dsum;
    float dxi, dyi, dzi;
    float dx2i, dy2i, dz2i; 
    float dx2dy2, dz2dx2, dz2dy2;

    int sidx, sidy, sidz;

    float * S;
    float * T;

    void init_sweep();
    void full_sweep();
    void inner_sweep();        

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