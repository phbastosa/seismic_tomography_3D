# ifndef CLASSIC_HPP
# define CLASSIC_HPP

# include "../eikonal.hpp"

class Classic : public Eikonal
{
private:
    
    float min(float v1, float v2);

public:

    void solve();
    void prepare_volumes();    
};

# endif