# 2D-Matrix-Template-Class-C++
A 2D matrix template class with many fast matrix operations.

## Preface
This project has arisen on a boring day that I had to spend sick at home. The goal was to create a simple matrix class that could do basic matrix operations.

However, since I felt quite an enjoyment when I started overengineering this class by introducing increasingly more templates, my plan changed: now I wanted to make the most overengineering project I had done for a while. If you feel like there is any functionality missing for the matrix class for accompanying functions, you may add a suggestion in the [suggestions file](SUGGESTIONS.md).

If you discover a bug, please check if your headers are up-to-date. In case they are, you can submit a bug in the [known bugs file](BUGS.md). Please first check if the bug is already known.

## Class-Operators Overview
You can use the following operators with the class itself:
|Operator        | Description |
|----------------|:----------- |
| = Matrix2D     | You can assign one matrix with another one. This will copy the matrix values to the destination matrix. |
| = const \_Elem | You can set the entire matrix to one specific value. typename \_Elem is the type of the matrix. |
| + Matrix2D     | You can add two matrices together. _**WARNING:** If you use this operator, make sure you are also assigning the result to a third matrix, otherwise there will be a memory leak!_ |
| - Matrix2D     | You can subtract two matrices. _**WARNING:** If you use this operator, make sure you are also assigning the result to a third matrix, otherwise there will be a memory leak!_ |
| += Matrix2D    | Use this to add two matrices and store the result in the first matrix. This is recommended over the single + operation, since no memory leak can occur. |
| -= Matrix2D    | Use this to subtract two matrices and store the result in the first matrix. This is recommended over the single - operation, since no memory leak can occur. |
| \*= Matrix2D   | Use this to multiply two matrices and store the result in the first matrix. |
| == Matrix2D    | Check if the values of two matrices with the same sized are identical. |
| != Matrix2D    | Check if the values of two matrices with the same sized are not identical. |
| []             | Index the array as if it would be 1-dimensional. Alternative functions at(index) or at(row, column) |
| cout()         | Not really an operator, but also quite useful. This function will output the matrix to the console. |

## Accompanying Functions
Please note that most accompanying functions in the Matrix-namespace are also declared and defined in the Matrix2D class. However, it is much simpler to use the functions outside the class. They will do the exact same calculations and there will be no performance to gain nor lose.

| Function Name | Arguments | Description |
| :------------ | :-------- | :---------- |
| Mul           | Matrix A, Matrix B, Matrix C | A and B will be multiplied and the result will be stored in C. Make sure you obey the matrix multiplication rule for the supplied matrices dimensions |
| Add           | Matrix A, Matrix B, Matrix C | A and B will be added and stored in C. All 3 matrices must have the same dimensions! |
| Add           | Matrix A, typename T scalar, Matrix B | The scalar value will be added to A and the result will be stored in B. Both matrices must have the same dimensions! |
| Sub           | Matrix A, Matrix B, Matrix C | A and B will be subtracted and stored in C. All 3 matrices must have the same dimensions! |
| Sub           | Matrix A, typename T scalar, Matrix B | The scalar value will be subtracted from A and the result will be stored in B. Both matrices must have the same dimensions! |
| Sub           | typename T scalar, Matrix A, Matrix B | A will be subtracted from the scalar and the result will be stored in B. Both matrices must have the same dimensions! |
| MulEBE        | Matrix A, Matrix B, Matrix C | Multiplying matrices A and B element-by-element, the result will be stored in C. All matrices must have the same dimensions. |
| MulEBE        | Matrix A, typename T scalar, Matrix B | Multiplying matrix A with the scalar value element-by-element, the result will be stored in B. All matrices must have the same dimensions. |
| DivEBE        | Matrix A, Matrix B, Matrix C | Dividing matrices A and B element-by-element, the result will be stored in C. All matrices must have the same dimensions. |
| DivEBE        | Matrix A, typename T scalar, Matrix B | Dividing matrix A by the scalar value element-by-element, the result will be stored in B. All matrices must have the same dimensions. |
| DivEBE        | typename T scalar, Matrix A, Matrix B | Dividing the scalar value by matrix A element-by-element, the result will be stored in B. All matrices must have the same dimensions. |
| PowEBE        | Matrix A, Matrix B, Matrix C | Calculating pow(A, B) element-by-element and storing result in C. |
| PowEBE        | Matrix A, typename T scalar, Matrix B | Calculating pow(A, scalar) element-by-element and storing result in B. |
| PowEBE        | typename T scalar, Matrix A, Matrix B | Calculating pow(scalar, A) element-by-element and storing result in B. |
| ExpEBE        | Matrix A, Matrix B | Calculating exp(A) element-by-element and storing the result in B. |
| SqrtEBE       | Matrix A, Matrix B | Calculating sqrt(A) element-by-element and storing the result in B. |
| Transpose90deg | Matrix A, Matrix B | Rotating A 90 degrees clockwise and storing the result in B. |
| Transpose180deg | Matrix A, Matrix B | Rotating A 180 degrees clockwise and storing the result in B. |
| Transpose270deg | Matrix A, Matrix B | Rotating A 270 degrees clockwise and storing the result in B. |

## Definitions
You can use the following definitions to adjust code compilation to your need:
| Name      | Description    |
| :-------- | :------------- |
| \_MATRIX_MANUAL_OPTIMISATION | You can define this label to use manual loop-unrolling for the multiplication operation. In my testings, I achieved a 40% speedup. |
| \_MATRIX_USE_FORCED_INLINES | You can define this label to force the compiler to inline addressing methods such as the .at(row, column) function. This may increase the program's speed. |
| \_MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N | This will limit the dimensions of the manual unrolled loops to max nxn * nxn. <br>The default value is 4. <br>Maximum: 5, minimum: 3. <br>It is advantages to keep this as small as possible. The max value is 5, however, some compilers may fail at such a high number to optimize away dead code, and it may slow down the program. If you can, keep it at 3 or 4. Everything below 3 will disable [_\_MATRIX_MANUAL_OPTIMISATION_] |
 
## Code Examples
### Getting Started
```cpp
#define _MATRIX_MANUAL_OPTIMISATION
#define _MATRIX_USE_FORCED_INLINES
#define _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N 3
#include "Matrix_Combined.h"

using namespace Matrix;

int main()
{
    // You can use any type and dimensions on a Matrix. You can 
    // learn more about the creation of a matrix at: "Creating a Matrix"
    Matrix2D<int, 3, 5> myIntMatrix =
    {
        { 1,  2,  3,  4,  5},
        { 6,  7,  8,  9, 10},
        {11, 12, 13, 14, 15},
    };

    // If you want to use an operation on two matrices, they
    // don't need to have the same type. However, make sure 
    // the dimensions are suitable for the desired operation.
    Matrix2D<float, 3, 5> myFloatMatrix =
    {
        { 1.1f,  2.2f,  3.3f,  4.4f,  5.5f},
        { 6.1f,  7.2f,  8.3f,  9.4f, 10.5f},
        {11.1f, 12.2f, 13.3f, 14.4f, 15.5f},
    };

    Matrix2D<float, 3, 5> myMatricAddionResult = 0.0f;

    // Add both matrices together and store them into [myMatricAddionResult].
    // You don't need to worry about the type of the matrices, the result will
    // automatically be casted to the correct type.
    Add(myIntMatrix, myFloatMatrix, myMatricAddionResult);

    return 0;
}
```
### Creating a Matrix
```cpp
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
```
### Multiplication
```cpp
#define _MATRIX_MANUAL_OPTIMISATION
#define _MATRIX_USE_FORCED_INLINES
#define _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N 3
#include "Matrix_Combined.h"
#include <iostream>

using namespace Matrix;

int main()
{
    Matrix2D<int, 2, 3> A =
    { 
        {1,2,3},
        {4,5,6} 
    };

    Matrix2D<int, 3, 2> B =
    {
        {9,8,7},
        {6,5,4},
    };

    Matrix2D<int, 2, 2> Result = 0.f;
    
    std::cout << "Performing matrix multiplication. A * B = Result" << std::endl;

    // Multiplying A * B and storing result in Result.
    Mul(A, B, Result);

    std::cout << std::endl << "Matrix A:" << std::endl;
    A.cout();

    std::cout << std::endl << "Matrix B:" << std::endl;
    B.cout();

    std::cout << std::endl << "Result:" << std::endl;
    Result.cout();

    return 0;
}
```
### Transposition
```cpp
#define _MATRIX_MANUAL_OPTIMISATION
#define _MATRIX_USE_FORCED_INLINES
#define _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N 3
#include "Matrix_Combined.h"
#include <iostream>

using namespace Matrix;

int main()
{
    Matrix2D<int, 2, 3> A =
    { 
        {1,2,3},
        {4,5,6} 
    };

    Matrix2D<int, 3, 2> Result90deg = 0;
    Matrix2D<int, 2, 3> Result180deg = 0;
    Matrix2D<int, 3, 2> Result270deg = 0;


    // Transposing Matrix 90 deg
    std::cout << "Performing matrix transposition. A 90 deg, strong result in Result90deg" << std::endl;

    Transpose90deg(A, Result90deg);

    std::cout << std::endl << "Matrix A:" << std::endl;
    A.cout();

    std::cout << std::endl << "Result 90 deg:" << std::endl;
    Result90deg.cout();

    // Transposing Matrix 180 deg
    std::cout << std::endl << std::endl << 
        "Performing matrix transposition. A 180 deg, strong result in Result180deg" << std::endl;

    Transpose180deg(A, Result180deg);

    std::cout << std::endl << "Matrix A:" << std::endl;
    A.cout();

    std::cout << std::endl << "Result 180 deg:" << std::endl;
    Result180deg.cout();


    // Transposing Matrix 180 deg
    std::cout << std::endl << std::endl << 
        "Performing matrix transposition. A 270 deg, strong result in Result270deg" << std::endl;

    Transpose270deg(A, Result270deg);

    std::cout << std::endl << "Matrix A:" << std::endl;
    A.cout();

    std::cout << std::endl << "Result 270 deg:" << std::endl;
    Result270deg.cout();


    return 0;
}
```
