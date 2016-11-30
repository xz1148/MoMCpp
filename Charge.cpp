#include "Charge.hpp"

Charge::Charge()
{
    mQ = 1.0;
    mx0 = 0.0;
    my0 = 0.0;
    mz0 = 0.0;
}

Charge::Charge(double Q)
{
    mQ = Q;
    mx0 = 0.0;
    my0 = 0.0;
    mz0 = 0.0;
}

Charge::Charge(double Q, double x, double y, double z)
{
    mQ = Q;
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

void Charge::unit_r(double r[3], double x, double y, double z)
{
    double r_x = x - mx0;
    double r_y = y - my0;
    double r_z = z - mz0;
    double mag_r_sqr = r_x*r_x + r_y*r_y + r_z*r_z;
    if (mag_r_sqr == 0)
    {
        r[0] = 0.0;
        r[1] = 0.0;
        r[2] = 0.0;
    }
    else
    {
        double mag_r = std::sqrt(mag_r_sqr);
        double u_x = r_x / mag_r;
        double u_y = r_y / mag_r;
        double u_z = r_z / mag_r;
        r[0] = u_x;
        r[1] = u_y;
        r[2] = u_z;
    }
}

double Charge::r(double x, double y, double z)
{
    double r_x = x - mx0;
    double r_y = y - my0;
    double r_z = z - mz0;
    double mag_r_sqr = r_x*r_x + r_y*r_y + r_z*r_z;
    double mag_r = std::sqrt(mag_r_sqr);
    return mag_r;
}
