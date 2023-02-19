# ifndef ACCURATE_FSM_HPP
# define ACCURATE_FSM_HPP

# include "../eikonal.hpp"

class Accurate_FSM : public Eikonal
{
private:

    void init_sweep();
    void full_sweep();
    void inner_sweep();        

    int i, j, k;

    int sgntz, sgntx, sgnty; 
    int sgnvz, sgnvx, sgnvy;

    int nb;

    float dsum;
    float dxi, dyi, dzi;
    float dx2i, dy2i, dz2i; 
    float dx2dy2, dz2dx2, dz2dy2;

    int nxx, nyy, nzz;
    int sidx, sidy, sidz;

    float dx, dy, dz;

public:

    void solve();
    void prepare_volumes();
};

# endif