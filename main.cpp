#include <iostream>
#include "Charge.hpp"



int main()
{


    Charge q1;
    Charge q2(2.0);
    Charge q3(3.0, 3.0, 3.0, 3.0);
    double unit_r[3];
    q3.unit_r(unit_r, 1.0, 1.0, 1.0);
    std::cout << unit_r[0] << unit_r[1] << unit_r[2] << std::endl;
    std::cout << 1/std::sqrt(3) << std::endl;
    std::cout << c0 << std::endl;
    std::cout << q3.r(1.0,1.0,1.0) << std::endl;
    std::cout << std::sqrt(12) << std::endl;

    return 0;
}
