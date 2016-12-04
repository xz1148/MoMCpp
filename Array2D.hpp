#ifndef ARRAY2D_HPP_INCLUDED
#define ARRAY2D_HPP_INCLUDED
#include <iostream>
class Array2D
{
private:
    double** mp_Data;  // the data
    int mDim1;
    int mDim2;
public:
//    double** mp_Data;  // the data
    Array2D();
//    Array2D(const Array2D& otherArray);
    Array2D(int dim1, int dim2);

    ~Array2D();  //destructor
    int GetDim1() const;
    int GetDim2() const;
    double* operator[](int i);  // the quick indexing
    friend std::ostream& operator<<(std::ostream& output, const Array2D& a2D);
};


#endif // ARRAY2D_HPP_INCLUDED
