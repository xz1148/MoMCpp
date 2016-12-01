#include "Greens.hpp"
#include <cmath>

double CalculateR(double x0, double y0, double z0, double x, double y, double z)
{
    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;
    double mag_r_sqr = dx*dx + dy*dy + dz*dz;
    double mag_r = std::sqrt(mag_r_sqr);
    return mag_r;
}

void Grad_Greens_r(std::complex<double> output[3],
                   double x0, double y0, double z0,
                   double x, double y, double z);

//double CalculateUnitR(double x0, double y0, double z0, double x, double y, double z )

//void Grad_Greens_r(std::complex<double> output[3], double x, double y, double z)
//{
//    double r = GetR(x, y, z);
//    double u_r[3];
//    GetUnit_R(u_r, x, y, z);
//    std::complex<double> alpha = 1i*mk0*r;
//    std::complex<double> Grad_Greens_r_numerator = (alpha+1.0)*std::exp(-alpha);
//    std::complex<double> Grad_Greens_r_denominator = 4*PI*r*r;
//    std::complex<double> Grad_G_r = Grad_Greens_r_numerator / Grad_Greens_r_denominator;
//    output[0] = Grad_G_r * u_r[0];
//    output[1] = Grad_G_r * u_r[1];
//    output[2] = Grad_G_r * u_r[2];
//}

