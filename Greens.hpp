/*this header contains the basic greens functions and their
useful transformations*/
#ifndef GREENS_HPP_INCLUDED
#define GREENS_HPP_INCLUDED
#include <cmath>
#include <complex>
double CalculateR(double x0,
                  double y0,
                  double z0,
                  double x,
                  double y,
                  double z);

void Grad_Greens_r(std::complex<double> output[3],
                   double x0, double y0, double z0,
                   double x, double y, double z);


#endif // GREENS_HPP_INCLUDED
