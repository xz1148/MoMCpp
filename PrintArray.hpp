#ifndef PRINTARRAY_HPP_INCLUDED
#define PRINTARRAY_HPP_INCLUDED
template<class T> void PrintArray(T* p_array, int length)
{
    for (int i=0; i<length; i++)
    {
        std::cout << p_array[i] << std::endl << std::flush;
    }
}

#endif // PRINTARRAY_HPP_INCLUDED
