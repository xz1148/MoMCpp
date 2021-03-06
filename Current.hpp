#ifndef CURRENT_HPP_INCLUDED
#define CURRENT_HPP_INCLUDED
#include <complex>

class Current
{
private:
    double mI; // this is the value of current
    double mk0; // the wave number in free space
    double mx0;
    double my0;
    double mz0; //the position of current
    double mux;
    double muy;
    double muz; //current point to (ux, uy, uz)

public:
    Current(double I, double k0,
            double x0, double y0, double z0,
            double vx, double vy, double vz);
    double GetR(double x, double y, double z);
    double Getk0() const;
    double GetI() const;
    double Getux() const;
    double Getuy() const;
    double Getuz() const;
    std::complex<double> Greens(double k0, double r);
    void A(std::complex<double> output[3], double x, double y, double z);
    void GetE(std::complex<double> output[3], double x, double y, double z);
};

class CubeCurrent:public Current
{
private:
    double mdl; // the edge length of cube
public:
    CubeCurrent(double I, double k0,
                double x0, double y0, double z0,
                double vx, double vy, double vz, double dl);
    void GetRemoteE(std::complex<double> output[3], double x, double y, double z);
    void GetLocalE(std::complex<double> output[3]);
};


#endif // CURRENT_HPP_INCLUDED
