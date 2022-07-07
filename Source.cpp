#define _MATRIX_MANUAL_OPTIMISATION
#define _MATRIX_USE_FORCED_INLINES
#define _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N 3
#include "Matrix_Combined.h"

using namespace Matrix;

int main()
{
    // Assigning a constant value will set all values of 
    // the matrix to the provided value:
    Matrix2D<double, 3, 8> allPi = 3.141592658979323;

    // You can also use an initializer list to specifically  
    // set each value of the
    Matrix2D<unsigned char, 2, 6> allPrime =
    {
        { 2,  3,  5,  7, 11, 13},
        {17, 19, 23, 29, 31, 37},
    };

    // You can also assign a matrix to another matrix. They
    // need to have the same dimensions:
    Matrix2D<unsigned char, 2, 6> alsoAllPrime = allPrime;

    // But they don't need to be of the same type
    Matrix2D<double, 2, 6> alsoAllPrimeButDoubles = allPrime;

    return 0;
}
