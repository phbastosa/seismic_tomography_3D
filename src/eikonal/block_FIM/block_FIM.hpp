# ifndef BLOCK_FIM_CUH
# define BLOCK_FIM_CUH

# include <array>
# include <vector>
# include <string>

# include "../eikonal.hpp"
# include "cuda_kernel.cuh"

class Block_FIM : public Eikonal
{
private:

    int padx, pady, padz;

    float * S;
    float * T;

    void expand_model();
    void reduce_model();
    void apply_model_mask();
    void extract_solution();
    void apply_source_time();

    CUDAMEMSTRUCT memoryStruct_;

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
