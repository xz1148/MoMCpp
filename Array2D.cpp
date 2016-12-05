#include "Array2D.hpp"
#include <cassert>

Array2D::Array2D()
{
    mDim1 = 2;
    mDim2 = 2;
    mp_Data = new double* [mDim1];  // give  a
//    for (int i=0; i<mDim1; i++)
//    {
//        mp_Data[i] = new double [mDim2];
//    }
    for (int i=0; i<mDim1; i++)
    {
        mp_Data[i] = new double [mDim2];
        for (int j=0; j<mDim2; j++)
        {
            mp_Data[i][j] = 0.0;
        }
    }
}

Array2D::Array2D(const Array2D& otherArray)
{
    mDim1 = otherArray.GetDim1();
    mDim2 = otherArray.GetDim2();
    mp_Data = new double* [mDim1];
    for (int i=0; i<mDim1; i++)
    {
        mp_Data[i] = new double [mDim2];
    }
    for (int i=0; i<mDim1; i++)
    {
        for (int j=0; j<mDim2; j++)
        {
            mp_Data[i][j] = otherArray.mp_Data[i][j];
        }
    }
}

Array2D::Array2D(int dim1, int dim2)
{
    mDim1 = dim1;
    mDim2 = dim2;
    mp_Data = new double* [mDim1];
    for (int i=0; i<mDim1; i++)
    {
        mp_Data[i] = new double [mDim2];
    }
    for (int i=0; i<mDim1; i++)
    {
        for (int j=0; j<mDim2; j++)
        {
            mp_Data[i][j] = 0.0;
        }
    }
}

Array2D::~Array2D()
{
    std::cout << "memory freeed" << std::endl << std::flush;
//    delete[] mp_Data[0];
    for (int i=0; i<mDim1; i++)
    {
        delete [] mp_Data[i];
    }
    delete [] mp_Data;
}

int Array2D::GetDim1() const
{
    return mDim1;
}

int Array2D::GetDim2() const
{
    return mDim2;
}

double Array2D::GetData(int dim1, int dim2) const
{
    assert(dim1>=0);
    assert(dim1<mDim1);

    assert(dim2>=0);
    assert(dim2<mDim2);

    return mp_Data[dim1][dim2];
}

double* Array2D::operator[](int i)
{
//    std::cout << i << std::endl;
    assert(i>-1);
    assert(i<mDim1);
    return mp_Data[i];
}

Array2D Array2D::operator+(double db)
{
    Array2D temp(mDim1, mDim2);  // create a temporary array
    for (int i=0; i<mDim1; i++)
    {
        for (int j=0; j<mDim2; j++)
        {
            temp[i][j] = mp_Data[i][j] + db;
        }
    }
    return temp;
}

Array2D Array2D::operator=(const Array2D& otherArray)
{
    assert(mDim1 = otherArray.mDim1);
    assert(mDim2 = otherArray.mDim2); // two arrays should have the same dimension

    for (int i=0; i<mDim1; i++)
    {
        for (int j=0; j<mDim2; j++)
        {
            mp_Data[i][j] = otherArray.GetData(i, j);
        }
    }
}

std::ostream& operator<<(std::ostream& output, const Array2D& a2D)
{
    int dim1 = a2D.mDim1;
    int dim2 = a2D.mDim2;
    for (int i=0; i<dim1; i++)
    {
        for (int j=0; j<dim2; j++)
        {
            output << a2D.mp_Data[i][j] << " ";
        }
        output << std::endl;
    }
    return output;
}
