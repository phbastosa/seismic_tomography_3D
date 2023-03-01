# ifndef EIKONAL_MODEL_HPP
# define EIKONAL_MODEL_HPP

# include "../model.hpp"

class Eikonal_model : public Model
{
public:

    void expand_fdm(float * input, float * output);
    void reduce_fdm(float * input, float * output);
    void expand_pad(float * input, float * output, int padx, int pady, int padz);
    void reduce_pad(float * input, float * output, int padx, int pady, int padz);
};


# endif