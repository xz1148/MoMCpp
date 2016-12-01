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
    double mk0;  // the wave number in free space, which is related to the frequency, and also it is a real number

public:
    Charge();
    Charge(double Q);
    Charge(double Q, double x, double y, double z);
    Charge(double Q, double k, double x, double y, double z);


    void GetUnit_R(double r[3], double x, double y, double z);
    // this function calculates the unit vector pointing from (mx0, my0, mz0) to (x, y, z) and it returns an array
    double GetR(double x, double y, double z);
    //the function calculates the distance form charge to position (x,y,z)
    void Grad_Greens_r(std::complex<double> output[3], double x, double y, double z);
    // this calculates the gradient of Greens function as a complex vector
    void GetE(std::complex<double> output[3], double x, double y, double z);
    // this function calculates the E field as a complex vector at ( x, y, z )
    double GetQ();
    double* returnQ();
};

#endif // CHARGE_HPP_INCLUDED
