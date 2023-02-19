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

private:

    int WARP = 32;

    void initialization();
    void apply_source_time();
    void apply_model_mask();
    void extract_solution();

    std::vector<std::vector<std::vector<float>>> speeds_;
    std::vector<std::vector<std::vector<float>>> answer_;

    CUDAMEMSTRUCT memoryStruct_;
};

# endif
