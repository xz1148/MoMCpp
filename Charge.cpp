#include "Charge.hpp"
#include "Greens.hpp"
#include <iostream>

Charge::Charge()
{
    mQ = 1.0;
    mk0 = 1.0;
    mx0 = 0.0;
    my0 = 0.0;
    mz0 = 0.0;
}

Charge::Charge(double Q)
{
    mQ = Q;
    mk0 = 1;
    mx0 = 0.0;
    my0 = 0.0;
    mz0 = 0.0;
}

Charge::Charge(double Q, double x, double y, double z)
{
    mQ = Q;
    mk0 = 1.0;
    mx0 = x;
    my0 = y;
    mz0 = z;
}

Charge::Charge(double Q, double k, double x, double y, double z)
{
    mQ = Q;
    mk0 = k;
    mx0 = x;
    my0 = y;
    mz0 = z;
}

double Charge::GetQ()
{
    return mQ;
}

double* Charge::returnQ()
{
    return &mQ;
}

void Charge::GetUnit_R(double r[3], double x, double y, double z)
{
    double r_x = x - mx0;
    double r_y = y - my0;
    double r_z = z - mz0;
    double mag_r = GetR(x, y, z);
    if (mag_r == 0)
    {
        r[0] = 0.0;
        r[1] = 0.0;
        r[2] = 0.0;
    }
    else
    {
        double u_x = r_x / mag_r;
        double u_y = r_y / mag_r;
        double u_z = r_z / mag_r;
        r[0] = u_x;
        r[1] = u_y;
        r[2] = u_z;
    }
}

double Charge::GetR(double x, double y, double z)
{
    double mag_r = CalculateR(mx0, my0, mz0, x, y, z);
    return mag_r;
}

void Charge::Grad_Greens_r(std::complex<double> output[3], double x, double y, double z)
{
    double r = GetR(x, y, z);
    double u_r[3];
    GetUnit_R(u_r, x, y, z);
    std::complex<double> alpha = 1i*mk0*r;
    std::complex<double> Grad_Greens_r_numerator = (alpha+1.0)*std::exp(-alpha);
    std::complex<double> Grad_Greens_r_denominator = 4*PI*r*r;
    std::complex<double> Grad_G_r = Grad_Greens_r_numerator / Grad_Greens_r_denominator;
    output[0] = Grad_G_r * u_r[0];
    output[1] = Grad_G_r * u_r[1];
    output[2] = Grad_G_r * u_r[2];
}

void Charge::GetE(std::complex<double> output[3], double x, double y, double z)
{
    double RightSide = -mQ / EPS0; // the right hand side of PDE
    std::complex<double> Grad_G_ur[3];
    Grad_Greens_r(Grad_G_ur, x, y, z); // complex vector
    output[0] = -RightSide * Grad_G_ur[0];
    output[1] = -RightSide * Grad_G_ur[1];
    output[2] = -RightSide * Grad_G_ur[2];
}
