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

private:

    int WARP = 32;
    int padx, pady, padz;

    float * S;
    float * T;
    float * K;
    float * nT;
    float * nK;

    void original_FIM();
    void open_acc_FIM();

    void apply_model_mask();
    void extract_solution();
    void apply_source_time();

    std::vector<std::vector<std::vector<float>>> speeds_;
    std::vector<std::vector<std::vector<float>>> answer_;

    CUDAMEMSTRUCT memoryStruct_;
};

# endif
