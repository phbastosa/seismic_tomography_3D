# ifndef BLOCK_FIM_CUH
# define BLOCK_FIM_CUH

# include <array>
# include <vector>
# include <string>

# include "../eikonal.hpp"
# include "cuda_kernel.cuh"

class Block_FIM : public Eikonal
{
public:
    
    void solve();
    void prepare_volumes();
    void destroy();
    
private:

    int padx, pady, padz;

    float * S;
    float * T;

    void apply_model_mask();
    void extract_solution();
    void apply_source_time();

    float min(float v1, float v2);

    std::vector<std::vector<std::vector<float>>> speeds_;
    std::vector<std::vector<std::vector<float>>> answer_;

    CUDAMEMSTRUCT memoryStruct_;
};

# endif
