# ifndef EIKONAL_MODEL_HPP
# define EIKONAL_MODEL_HPP

# include "../model.hpp"

class Eikonal_model : public Model
{
public:

    float * expand_fdm(float * volume);
    float * reduce_fdm(float * volume);
    float * expand_pad(float * volume, int padx, int pady, int padz);
    float * reduce_pad(float * volume, int padx, int pady, int padz);
};


# endif