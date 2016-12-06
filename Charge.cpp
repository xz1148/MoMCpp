#include "Charge.hpp"
#include "Greens.hpp"
#include "GLQuad.hpp"
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
    mk0 = 1.0;
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

double Charge::Getx0() const
{
    return mx0;
}

double Charge::Gety0() const
{
    return my0;
}

double Charge::Getz0() const
{
    return mz0;
}

double Charge::GetQ() const
{
    return mQ;
}

double* Charge::returnQ()
{
    return &mQ;
}

void Charge::SetXYZ(double x0, double y0, double z0)
{
    mx0 = x0;
    my0 = y0;
    mz0 = z0;
}

void Charge::SetQ(double Q)
{
    mQ = Q;
}

void Charge::Setk0(double k0)
{
    mk0 = k0;
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
//    std::cout << mQ <<  std::endl << std::flush;
//    std::cout << GetQ() << std::endl << std::flush;
    double RightSide = -mQ / EPS0; // the right hand side of PDE
    std::complex<double> Grad_G_ur[3];
    Grad_Greens_r(Grad_G_ur, x, y, z); // complex vector
    output[0] = -RightSide * Grad_G_ur[0];
    output[1] = -RightSide * Grad_G_ur[1];
    output[2] = -RightSide * Grad_G_ur[2];
}


SquareCharge::SquareCharge():m_xyz_c(1,3), mnt(2)
{
    int N_charges = mnt*mnt;
    mdl = 1.0;
    mdir = 1;
    mk0 = 1.0;
    mQ = 1.0;
    m_xyz_c[0][0] = 0.0;m_xyz_c[0][1] = 0.0;m_xyz_c[0][1] = 0.0;
    Array2D Charge_xyz_Temp(mnt*mnt, 3);  //stores the coordinate of charges
    Array2D Charge_Wts(1, mnt*mnt);  // stores the weights of charges
    AbsAndWtsSquare(mnt, mdl, mdir, m_xyz_c, Charge_xyz_Temp, Charge_Wts); //calculates all the xyzs
    mCharge = new Charge [N_charges];
    for (int i=0; i<mnt*mnt; i++)
    {
        mCharge[i].SetXYZ(Charge_xyz_Temp[i][0], Charge_xyz_Temp[i][1], Charge_xyz_Temp[i][2]);
        mCharge[i].SetQ(mQ*Charge_Wts[0][i]);
        mCharge[i].Setk0(mk0);
    }
}

SquareCharge::SquareCharge(int nt):m_xyz_c(1,3)
{
    mnt = nt;
//    m_xyz_c(1,3);
    mdl = 1.0;
    mdir = 1;
    mk0 = 1.0;
    mQ = 1.0;

    m_xyz_c[0][0] = 0.0;m_xyz_c[0][1] = 0.0;m_xyz_c[0][1] = 0.0;
    Array2D Charge_xyz_Temp(mnt*mnt, 3);  //stores the coordinate of charges
    Array2D Charge_Wts(1, mnt*mnt);  // stores the weights of charges
    AbsAndWtsSquare(mnt, mdl, mdir, m_xyz_c, Charge_xyz_Temp, Charge_Wts); //calculates all the xyzs
    mCharge = new Charge [mnt*mnt];
    for (int i=0; i<mnt*mnt; i++)
    {
        mCharge[i].SetXYZ(Charge_xyz_Temp[i][0], Charge_xyz_Temp[i][1], Charge_xyz_Temp[i][2]);
        mCharge[i].SetQ(mQ*Charge_Wts[0][i]);
        mCharge[i].Setk0(mk0);
    }
}

SquareCharge::SquareCharge(int nt, int dir, double dl, double k0, double Q):m_xyz_c(1,3)
{
    mnt = nt;
    mdir = dir;
    mdl = dl;
    mk0 = k0;
    mQ = Q;
    m_xyz_c[0][0] = 0.0;m_xyz_c[0][1] = 0.0;m_xyz_c[0][1] = 0.0;
    Array2D Charge_xyz_Temp(mnt*mnt, 3);  //stores the coordinate of charges
    Array2D Charge_Wts(1, mnt*mnt);  // stores the weights of charges
    AbsAndWtsSquare(mnt, mdl, mdir, m_xyz_c, Charge_xyz_Temp, Charge_Wts); //calculates all the xyzs
    mCharge = new Charge [mnt*mnt];
    for (int i=0; i<mnt*mnt; i++)
    {
        mCharge[i].SetXYZ(Charge_xyz_Temp[i][0], Charge_xyz_Temp[i][1], Charge_xyz_Temp[i][2]);
        mCharge[i].SetQ(mQ*Charge_Wts[0][i]);
        mCharge[i].Setk0(mk0);
    }

}

SquareCharge::SquareCharge(int nt,
                           int dir,
                           double dl,
                           double k0,
                           double Q,
                           double x0,
                           double y0,
                           double z0):m_xyz_c(1,3)
{
    mnt = nt;
    mdir = dir;
    mdl = dl;
    mk0 = k0;
    mQ = Q;
    m_xyz_c[0][0] = x0;m_xyz_c[0][1] = y0;m_xyz_c[0][2] = z0;
    Array2D Charge_xyz_Temp(mnt*mnt, 3);  //stores the coordinate of charges
    Array2D Charge_Wts(1, mnt*mnt);  // stores the weights of charges
    AbsAndWtsSquare(mnt, mdl, mdir, m_xyz_c, Charge_xyz_Temp, Charge_Wts); //calculates all the xyzs
    mCharge = new Charge [mnt*mnt];
    for (int i=0; i<mnt*mnt; i++)
    {
        mCharge[i].SetXYZ(Charge_xyz_Temp[i][0], Charge_xyz_Temp[i][1], Charge_xyz_Temp[i][2]);
        mCharge[i].SetQ(mQ*Charge_Wts[0][i]);
        mCharge[i].Setk0(mk0);
    }
}



int SquareCharge::GetChargeNumber() const
{
    return mnt*mnt;
}

void SquareCharge::GetE(std::complex<double> output[3], const Array2D& xyz)
{
    int N_charge = mnt*mnt;
    double x = xyz.GetData(0,0);
    double y = xyz.GetData(0,1);
    double z = xyz.GetData(0,2);
    std::complex<double> E_temp[3];
    std::complex<double> E_sum[3];
    for (int i=0; i<3; i++)
    {
        E_temp[i] = std::complex<double>(0.0, 0.0);
        E_sum[i] = std::complex<double>(0.0, 0.0);
    }

    for (int i=0; i<N_charge; i++)
    {
        mCharge[i].GetE(E_temp, x, y, z);
        E_sum[0] += E_temp[0];
        E_sum[1] += E_temp[1];
        E_sum[2] += E_temp[2];
    }
    output[0] = E_sum[0]; output[1] =  E_sum[1]; output[2] = E_sum[2];
}

std::ostream& operator<<(std::ostream& output, const SquareCharge& SqrCharge)
{
    double Q_sum = 0.0;
    int ChargeNumber = SqrCharge.GetChargeNumber();
    output << "The square contains " << SqrCharge.GetChargeNumber() << " Charges:" << std::endl;
    output << "The coordinate of the charges are: " << std::endl;
    for (int i=0; i<ChargeNumber; i++)
    {
        output << SqrCharge.mCharge[i].Getx0() << ", ";
        output << SqrCharge.mCharge[i].Gety0() << ", ";
        output << SqrCharge.mCharge[i].Getz0() << std::endl;
    }
    output << "And the charges are:" << std::endl;
    for (int i=0; i<ChargeNumber; i++)
    {
        output << SqrCharge.mCharge[i].GetQ() << std::endl;
    }
    for (int i=0; i<ChargeNumber; i++)
    {
        Q_sum = Q_sum + SqrCharge.mCharge[i].GetQ();
    }
    output << "And the total charge is: " << Q_sum << std::endl;

    return output;
}

SquareCharge::~SquareCharge()
{
    delete[] mCharge;
}
