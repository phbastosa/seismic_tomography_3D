# include <string>
# include <chrono>
# include <iostream>

# include "least_squares/least_squares.hpp"

int main(int argc, char **argv)
{
    auto ti = std::chrono::system_clock::now();

    Tomography * tomography = new Least_squares();

    tomography->parameters = std::string(argv[1]);

    tomography->set_parameters();









    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds;

    elapsed_seconds = tf - ti;

    std::cout<<"\nTomography run time: "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}
