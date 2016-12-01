#include "Current.hpp"
#include "Greens.hpp"
#include <cmath>
#include <iostream>
#include "EMConst.hpp"


Current::Current(double I, double k0,
                 double x0, double y0, double z0,
                 double vx, double vy, double vz)
{
    mI = I;
    mk0 = k0;
    mx0 = x0;
    my0 = y0;
    mz0 = z0;
    double mag_v = CalculateR(0,0,0,vx, vy, vz);
    if (mag_v == 0)
    {
        std::cout << "Warning: the direction of current density is zero vector." << std::endl;
        std::cout << "Please reset the current direction." << std::endl;
        mux = 0.0;
        muy = 0.0;
        muz = 0.0;
    }
    else
    {
        mux = vx/mag_v;
        muy = vy/mag_v;
        muz = vz/mag_v;
    }

};

double Current::GetR(double x, double y, double z)
{
    double mag_r = CalculateR(mx0, my0, mz0, x, y, z);
    return mag_r;
}

std::complex<double> Current::Greens(double k0, double r)
{
    if (r==0.0)
    {
        std::cerr << "Warning: the field is estimating at singular point" << std::endl;
    }
    std::complex<double> alpha = -1.0i*k0*r;
    std::complex<double> Greens_numerator = -std::exp(alpha);
    std::complex<double> Greens_denominator = 4.0*PI*r;
    std::complex<double> Greens_total = Greens_numerator / Greens_denominator;
    return Greens_total;
};


void Current::A(std::complex<double> output[3], double x, double y, double z)
{
    double right_side = -MU0*mI;  //right side of PDE
    double r = GetR(x, y, z);
    std::complex<double> G = Greens(mk0, r);
    output[0] = right_side * G * mux;
    output[1] = right_side * G * muy;
    output[2] = right_side * G * muz;
}


void Current::GetE(std::complex<double> output[3], double x, double y, double z)
{
    std::complex<double> A_vec[3];
    A(A_vec, x, y, z);
    double omega = mk0 / std::sqrt(EPS0 * MU0);
    std::complex<double> pure_imag(0.0,1.0);
    output[0] = -pure_imag*omega*A_vec[0];
    output[1] = -pure_imag*omega*A_vec[1];
    output[2] = -pure_imag*omega*A_vec[2];
}
