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

double Current::Getk0() const
{
    return mk0;
}

double Current::GetI() const
{
    return mI;
}


double Current::Getux() const {return mux;}
double Current::Getuy() const {return muy;}
double Current::Getuz() const {return muz;}



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
    double omega = mk0 / std::sqrt(EPS0 * MU0);  //calculate angular frequency
    std::complex<double> pure_imag(0.0,1.0);
    output[0] = -pure_imag*omega*A_vec[0];
    output[1] = -pure_imag*omega*A_vec[1];
    output[2] = -pure_imag*omega*A_vec[2];
}


CubeCurrent::CubeCurrent(double I, double k0,
            double x0, double y0, double z0,
            double vx, double vy, double vz, double dl):Current(I,k0,x0,y0,z0,vx,vy,vz)
{
    mdl = dl;
}

void CubeCurrent::GetRemoteE(std::complex<double> output[3], double x, double y, double z)
{
    double volume = mdl*mdl*mdl;
    GetE(output, x, y, z);
    for (int i=0; i<3; i++)
    {
        output[i] = output[i]*volume;
    }
}

void CubeCurrent::GetLocalE(std::complex<double> output[3])
{
    double volume=mdl*mdl*mdl;
    double k0 = Getk0();
    double r_eq = std::pow(3.0/(4.0*PI)*volume, 1.0/3.0); // the equvalent radius of each element;
    double omega = Getk0() / std::sqrt(EPS0 * MU0);  //calculate angular frequency
    std::complex<double> C0 = 1.0i*omega*MU0*GetI(); //now mI is the current density
    std::complex<double> alpha = 1.0i*k0;            //1j*k0
    std::complex<double> term1 = r_eq * std::exp(-alpha*r_eq) / alpha;
    std::complex<double> term2 = std::exp(-alpha*r_eq) / (k0*k0);
    std::complex<double> term3 = 1 / (k0*k0);
    std::complex<double> E = C0*(term1 - term2 + term3);
    output[0] = E*Getux(); output[1] = E*Getuy(); output[2] = E*Getuz();
}
