#include <iostream>
#include "Charge.hpp"
#include "PrintArray.hpp"
#include "Current.hpp"
#include "GLQuad.hpp"
#include "Array2D.hpp"

using namespace std;



int main()
{
//    double freq = 300e6;
//    double omega = 2*PI*freq;
//    double k0 = omega * std::sqrt(MU0 * EPS0);
//    double x0 = 1.0;
//    double y0 = 1.0;
//    double z0 = 1.0;
//    double Q = 1.0;
//    Charge q1(Q, k0, x0, y0, z0);
//    std::cout << q1.GetQ() << std::endl;
//    double x = 6.0;
//    double y = 2.0;
//    double z = 4.0;
//
//
//    std::complex<double> output[3];
//    q1.GetE(output, x, y, z);
//    PrintArray< std::complex<double> >(output, 3);
//    std::cout << q1.GetR(x, y, z) << std::endl;
//
//    double ux = 1.0;
//    double uy = 1.0;
//    double uz = 3.0;
//    double I = 1.0;
//    Current J(I, k0, x0, y0, z0, ux, uy, uz);  //give all the properties to current density
//    std::complex<double> E[3];
//    J.GetE(E, x, y, z);
//    std::cout << J.GetR(x, y, z) << std::endl;
//    PrintArray< std::complex<double> > (E, 3);
//    std::cout << sizeof(J) << std::endl;
//    std::cout << sizeof(q1) << std::endl;
//    std::cout << sizeof(E) << std::endl;

//    int nt = 4;
//    double t[4];
//    double wts[4];
//    double A = 2.0;
//    double B = 5.0;
//    double sum;
//    AbsAndWtsAB(nt, A, B, t, wts);
//    for (int i =0; i<nt; i++)
//    {
//        cout << t[i] << endl;
//        cout << wts[i] << endl;
//        sum+= wts[i];
//    }
//    std::cout << output[0] << std::endl;
//    std::cout << output[1] << std::endl;
//    std::cout << output[2] << std::endl;
//    cout << sum << endl;
//    double a = {{1,2},{2,3}};

    int nt = 2;
    double dl = 0.5;
    int direction = 1;
    Array2D xyz_c(1,3);
    Array2D xyz_out(4,3);
    Array2D wts(1,3);
    AbsAndWtsSquare(nt, dl, direction, xyz_c, xyz_out, wts);

    Array2D x(5,3);


    x[4][1] = 2.0;
    cout << x << endl;
//    cout << x << endl;

//    double** x;
//    x = new double* [10];;
//    for (int i=0; i<10; i++)
//    {
//        x[i] = new double [2];
//        for (int j=0; j<2; j++)
//            x[i][j] = i+j;
//    }
//    for (int i=0; i<10; i++)
//    {
//        delete[] x[i];
//    }
//    delete[] x;

    return 0;
}
