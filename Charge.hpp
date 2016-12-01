#ifndef CHARGE_HPP_INCLUDED
#define CHARGE_HPP_INCLUDED
#include <cmath>
#include "EMConst.hpp"

std::complex<double> Greens(double k0, double r);
void GradGreensr(std::complex<double> dGdr_r[3], double k0, double r);


class Charge
{
private:
    double mQ;  //number of charges
    double mx0;
    double my0;
    double mz0;  // the position of the charge
public:
    Charge();
    Charge(double Q);
    Charge(double Q, double x, double y, double z);
    void unit_r(double r[3], double x, double y, double z);
    // this function calculates the unit vector pointing from (mx0, my0, mz0) to (x, y, z) and it returns an array
    double r(double x, double y, double z);
    //the function calculates the distance form charge to position (x,y,z)
    double GetQ();
    double* returnQ();
};

#endif // CHARGE_HPP_INCLUDED
