# ifndef INOUT_HPP
# define INOUT_HPP

# include <string>
# include <sstream>
# include <fstream>
# include <iostream>

class InOut
{
private:

public:    
    static void readBinaryFloat(std::string path, float *array, int n);
    static void writeBinaryFloat(std::string path, float *array, int n);

    static std::string catchParameter(std::string target, std::string file);
};

# endif
