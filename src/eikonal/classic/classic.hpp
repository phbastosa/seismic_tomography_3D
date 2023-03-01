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

public:

    void solve();
    void prepare_volumes();
    void destroy();    
};

# endif