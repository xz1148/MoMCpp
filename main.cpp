#include <iostream>
#include "Charge.hpp"
#include "PrintArray.hpp"
#include "Current.hpp"


int main()
{
    double freq = 300e6;
    double omega = 2*PI*freq;
    double k0 = omega * std::sqrt(MU0 * EPS0);
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;
    double Q = 1.0;
    Charge q1(Q, k0, x0, y0, z0);
    std::cout << q1.GetQ() << std::endl;
    double x = 6.0;
    double y = 2.0;
    double z = 4.0;


    std::complex<double> output[3];
    q1.GetE(output, x, y, z);
    PrintArray< std::complex<double> >(output, 3);
    std::cout << q1.GetR(x, y, z) << std::endl;

    double ux = 1.0;
    double uy = 1.0;
    double uz = 3.0;
    double I = 1.0;
    Current J(I, k0, x0, y0, z0, ux, uy, uz);  //give all the properties to current density
    std::complex<double> E[3];
    J.GetE(E, x, y, z);
    std::cout << J.GetR(x, y, z) << std::endl;
    PrintArray< std::complex<double> > (E, 3);
    std::cout << sizeof(J) << std::endl;
    std::cout << sizeof(q1) << std::endl;
    std::cout << sizeof(E) << std::endl;




//
//    std::cout << output[0] << std::endl;
//    std::cout << output[1] << std::endl;
//    std::cout << output[2] << std::endl;

    return 0;
}
