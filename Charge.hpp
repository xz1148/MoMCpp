#ifndef CHARGE_HPP_INCLUDED
#define CHARGE_HPP_INCLUDED
#include <cmath>
#include "EMConst.hpp"
#include "Array2D.hpp"

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
    void SetXYZ(double x0, double y0, double z0);
    void SetQ(double Q);
    void Setk0(double k0);
    double Getx0() const;
    double Gety0() const;
    double Getz0() const;
    void GetUnit_R(double r[3], double x, double y, double z);
    // this function calculates the unit vector pointing from (mx0, my0, mz0) to (x, y, z) and it returns an array
    double GetR(double x, double y, double z);
    //the function calculates the distance form charge to position (x,y,z)
    void Grad_Greens_r(std::complex<double> output[3], double x, double y, double z);
    // this calculates the gradient of Greens function as a complex vector
    void GetE(std::complex<double> output[3], double x, double y, double z);
    // this function calculates the E field as a complex vector at ( x, y, z )
    double GetQ() const;
    double* returnQ();
};

class SquareCharge
{
private:
    Charge* mCharge; // the charges in a square
    Array2D m_xyz_c;
    int mnt;  // the resolution
    int mdir;  // direction of the square face, reference to AbsAndWtsSquare()
    double mdl;
    double mk0;
    double mQ; // the charge density, so the total amount of charges are mQ*ds, ds is area
public:
    int GetChargeNumber() const;
    SquareCharge();
    SquareCharge(int nt);
    SquareCharge(int nt, int dir, double dl, double k0, double Q);
    SquareCharge(int nt, int dir, double dl, double k0, double Q, double x0, double y0, double z0);
    void GetE(std::complex<double> output[3], const Array2D& xyz);
    ~SquareCharge();
    friend std::ostream& operator<<(std::ostream& output, const SquareCharge& SqrCharge);
};

#endif // CHARGE_HPP_INCLUDED
