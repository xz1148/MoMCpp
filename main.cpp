#include <iostream>
#include "Charge.hpp"
#include "PrintArray.hpp"
#include "Current.hpp"
#include "GLQuad.hpp"
#include "Array2D.hpp"

using namespace std;



int main()
{

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

//    int nt = 5;
//    double dl = 0.25;
//    int direction = 2;
//    Array2D xyz_c(1,3);
//    xyz_c[0][0] = 0;
//    xyz_c[0][1] = 1;
//    xyz_c[0][2] = 3;
//    Array2D xyz_out(nt*nt,3);
//    Array2D xyz_out2(nt*nt, 3);
//
//
//
//
//    Array2D wts(1, nt*nt);
//    for (int i=0; i<1; i++)
//    {
//        AbsAndWtsSquare(nt, dl, direction, xyz_c, xyz_out, wts);
//    }
//    xyz_out2 = xyz_out;
//    double sum;
//    for (int i=0; i<nt*nt; i++)
//    {
//        sum+=wts[0][i];
//    }
//    xyz_out[1][1] = 3.0;
//    cout << sum << endl;
//    cout << xyz_out2[1][1] << endl;
//    cout << xyz_out[1][1] << endl;
//    cout << endl;
//    cout << xyz_c << endl;
////    cout << x+3.0 << endl << flush;
   // cout << x[2][2] << endl << flush;
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
    double freq = 100e6;
    double omega = 2*PI*freq;
    double k0 = omega * std::sqrt(MU0 * EPS0);
    int nt = 10;
    int dir = 1;
    double dl = 0.02;
    double Q = 1.0;
    double x = 1.0;
    double y = 2.0;
    double z = 3.0;
//    SquareCharge sqc1(nt, dir, dl, k0, Q, x, y, z);
    Charge sc1(dl*dl, k0, x, y, z);
    Array2D xyz_E(1,3);
    xyz_E[0][0] = 1.0; xyz_E[0][1] = 2.0; xyz_E[0][2] = 3+dl/2;
    complex<double> cplx[3];
    complex<double> cplx1[3];
    complex<double> cplx2[3];
    complex<double> cplx3[3];
    complex<double> cplx4[3];


    sc1.GetE(cplx1, xyz_E[0][0], xyz_E[0][1], xyz_E[0][2]);
//    for (int i = 0; i<500000; i++)
//    {
    CubeCurrent cbc1(1.0,k0,0,0,0,0,0,1,dl);
    SquareCharge sqc1(nt, dir, dl, k0, Q, x, y, z);
    sqc1.GetE(cplx, xyz_E);
    cbc1.GetRemoteE(cplx2, xyz_E[0][0], xyz_E[0][1], xyz_E[0][2]);
//    }

//    CubeCurrent cbc1(1.0,k0,0,0,0,3,0,0,dl);
    Current c1(1.0, k0, 0,0,0,0,0,1);
    c1.GetE(cplx3, xyz_E[0][0], xyz_E[0][1], xyz_E[0][2]);
    cbc1.GetLocalE(cplx4);
    for (int i=0; i<3; i++)
    {
        cout << setprecision(15) << cplx1[i] << endl;
//        cout << setprecision(15) << cplx1[i] << endl;
//        cout << setprecision(15) << cplx2[i] << endl;
    }
    for (int i=0; i<3; i++)
    {
        cout << setprecision(15) << cplx2[i] << endl;
    }
    for (int i=0; i<3; i++)
    {
        cout << setprecision(15) << cplx3[i]*dl*dl*dl << endl;
    }
    for (int i=0; i<3; i++)
    {
        cout << setprecision(15) << cplx4[i] << endl;
    }

//    cout << sqc1 << endl;
//    cout << sqc1 << endl;

//    cout << cplx[0] << endl;
//    cout << sqc1 << endl;
    return 0;
}
