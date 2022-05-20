# include "inout.hpp"

int main(int argc, char **argv)
{
    int n = 400;

    float * ricker = new float[n]();

    auto * io = new InOut();

    std::string filePath = "ricker.bin";

    InOut::readBinaryFloat(filePath,ricker,n);

    std::string newFilePath = "ricker_" + InOut::toString(n) + "_samples.bin";

    InOut::writeBinaryFloat(newFilePath,ricker,n);

    std::cout<<"\nTest concluded successfully!"<<std::endl;

    return 0;
}
