#include <iostream>
#include "Charge.hpp"



int main()
{
    double freq = 300e6;
    double omega = 2*PI*freq;
    double k0 = omega * std::sqrt(MU0 * EPS0);
    double x0 = 3.0;
    double y0 = 3.0;
    double z0 = 3.0;
    double Q = 1.0;
    Charge q1(Q, x0, y0, z0);

    double x = 1.0;
    double y = 1.0;
    double z = 1.0;
    double r = q1.r(x,y,z);

    std::cout << r << std::endl;
    std::cerr << "this is error" << std::endl;
    std::cout << Greens(k0, 0.0) << std::endl;


    return 0;
}
