# include <string>
# include <sstream>
# include <fstream>
# include <iostream>

# include "inout.hpp"

std::string InOut::toString(int n)
{
    std::ostringstream s;
    s << n;
    return s.str();
}

void InOut::readBinaryFloat(std::string path, float *array, int n)
{
    std::ifstream file(path, std::ios::binary);
    file.read((char *) array, n * sizeof(float));
    file.close();    
}

void InOut::writeBinaryFloat(std::string path, float *array, int n)
{
    std::ofstream file(path, std::ios::binary);
    file.write((char *) array, n * sizeof(float));
    file.close();
}
