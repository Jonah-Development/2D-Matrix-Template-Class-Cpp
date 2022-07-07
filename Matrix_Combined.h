/*
    Matrix2D, a 2D matrix template class with many fast matrix operations
    Copyright (C) 2022 Jonah Kable

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
    IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <vector>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <string>
#include <sstream>

// using manually unrolled loop for common size matrices will increase speed by ~40%,
// so we use it if we are in release mode
#if !defined(_MATRIX_MANUAL_OPTIMISATION) && defined(NDEBUG)
#define _MATRIX_MANUAL_OPTIMISATION
#endif // !!defined(_MATRIX_MANUAL_OPTIMISATION) && defined(NDEBUG)

#ifdef _INLINE
#define _MATRIX_INLINE_TMP
#endif

#ifdef _MATRIX_USE_FORCED_INLINES
#define _INLINE __forceinline
#elif defined(_MATRIX_DISABLE_INLINE)
#define _INLINE
#else 
#define _INLINE inline
#endif // defined()

#ifndef _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N
#define _MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N 4
#endif // !_MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N


namespace Matrix
{
    template<typename _Elem, size_t _Rows, size_t _Cols>
    class Matrix2D
    {
    public:
        Matrix2D(void) = delete;
        template <typename _ElemB> Matrix2D(const _ElemB);
        Matrix2D(const Matrix2D&);
        template <typename _ElemB> Matrix2D(Matrix2D<_ElemB, _Rows, _Cols>&);
        Matrix2D(const std::initializer_list<std::initializer_list<_Elem>>);
        Matrix2D(_Elem*);
        ~Matrix2D(void);

        Matrix2D& operator = (const Matrix2D&);
        Matrix2D& operator = (_Elem*);
        Matrix2D& operator = (const _Elem);
        _Elem* operator + (Matrix2D&);
        _Elem* operator - (Matrix2D&);
        _Elem* operator * (Matrix2D&) = delete; // unable to ensure that both matrices have the correct size. Please use Matrix::Mul()
        _Elem* operator / (Matrix2D&) = delete;
        _Elem* operator % (Matrix2D&) = delete;
        _Elem* operator << (Matrix2D&) = delete;
        _Elem* operator >> (Matrix2D&) = delete;

        void operator += (Matrix2D&);
        void operator -= (Matrix2D&);
        template <typename _Elm2, size_t _Rows2, size_t _Cols2>
        void operator *= (Matrix2D<_Elm2, _Rows2, _Cols2>&);
        bool operator == (Matrix2D&);
        bool operator != (Matrix2D&);
        void operator /= (Matrix2D&) = delete;
        void operator %= (Matrix2D&) = delete;
        void operator <<= (Matrix2D&) = delete;
        void operator >>= (Matrix2D&) = delete;
        Matrix2D& operator ++ (void) = delete;
        Matrix2D operator ++ (int) = delete;
        Matrix2D& operator -- (void) = delete;
        Matrix2D operator -- (int) = delete;
        bool operator | (Matrix2D&) = delete;
        bool operator & (Matrix2D&) = delete;
        bool operator < (Matrix2D&) = delete;
        bool operator <= (Matrix2D&) = delete;
        bool operator > (Matrix2D&) = delete;
        bool operator >= (Matrix2D&) = delete;
        _Elem* operator |= (Matrix2D&) = delete;
        _Elem* operator &= (Matrix2D&) = delete;
        _Elem* operator ^= (Matrix2D&) = delete;
        void operator ^ (Matrix2D&) = delete;

        _INLINE _Elem& operator [] (const size_t index) { return p_mat[index]; }
        _INLINE _Elem& at(const size_t index) { return p_mat[index]; }
        _INLINE _Elem& at(const size_t row, const size_t column) { return p_mat[POS_XY(row, column)]; }

        template<typename _Cast = _Elem>
        void cout(void);

        template <typename element, size_t rows, size_t columns>
        static void Print(Matrix2D<element, rows, columns>& mat);

        //
        // All functions needed to perform all kinds of matrix math stuff
        //

        // A * B = C
        template<typename _ElemA, size_t _RowsA, size_t _ColsA, typename _ElemB, size_t _ColsB, typename _ElemC>
        static void Mul(Matrix2D<_ElemA, _RowsA, _ColsA>& A, Matrix2D<_ElemB, _ColsA, _ColsB>& B, Matrix2D<_ElemC, _RowsA, _ColsB>& C);

        // A + B = C
        template<typename _ElemA, typename _ElemB>
        static void Add(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C);

        // A + scalar = B
        template<typename _ElemA, typename _ElemB>
        static void Add(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A - B = C
        template<typename _ElemA, typename _ElemB>
        static void Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C);

        // A - scalar = B
        template<typename _ElemA, typename _ElemB>
        static void Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B);

        // scalar - A = B
        template<typename _ElemA, typename _ElemB>
        static void Sub(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A .* B = C
        template<typename _ElemA, typename _ElemB>
        static void MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C);

        // A .* scalar = B
        template<typename _ElemA, typename _ElemB>
        static void MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A ./ B = C
        template<typename _ElemA, typename _ElemB>
        static void DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C);

        // A ./ scalar = B
        template<typename _ElemA, typename _ElemB>
        static void DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B);

        // scalar ./ A = B
        template<typename _ElemA, typename _ElemB>
        static void DivEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A .^ B = C
        template<typename _ElemA, typename _ElemB>
        static void PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C);

        // A .^ scalar = B
        template<typename _ElemA, typename _ElemB>
        static void PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B);

        // scalar .* A = B
        template<typename _ElemA, typename _ElemB>
        static void PowEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A .exp = B
        template<typename _ElemA>
        static void ExpEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B);

        // A .sqrt = B
        template<typename _ElemA>
        static void SqrtEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B);

    private:
        _Elem* p_mat;
        const static size_t c_elements = _Rows * _Cols;
        const static size_t c_rows = _Rows;
        const static size_t c_cols = _Cols;

    private:

        constexpr static _INLINE size_t POS_XY(const size_t row, const size_t column)
        {

            //std::cout << "Checking dimensions, column: " << column << ",\trows: " << row << std::endl;
            //std::cout << "            expected: _Cols: " << _Cols << ",\t_Rows: " << _Rows << std::endl;

            if (column > _Cols) throw std::runtime_error("X out of range!");
            if (row > _Rows) throw std::runtime_error("Y out of range!");

            return (row * _Cols) + column;
        }

        bool SameDimensions(const Matrix2D&);
    };

    template <class T>
    constexpr static auto SIZE_OF_MATRIX2D = sizeof(Matrix2D<T, 1, 1>);

    template <size_t _Width>
    constexpr static _INLINE size_t POS_XY(const size_t row, const size_t column)
    {
        return (row * _Width) + column;
    }

    // A * B = C
    template<typename _ElemA, size_t _RowsA, size_t _ColsA, typename _ElemB, size_t _ColsB, typename _ElemC>
    static void Mul(Matrix2D<_ElemA, _RowsA, _ColsA>& A, Matrix2D<_ElemB, _ColsA, _ColsB>& B, Matrix2D<_ElemC, _RowsA, _ColsB>& C);

    // A + B = C
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void Add(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_ElemC, _Rows, _Cols>& C);

    // A + scalar = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void Add(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // A - B = C
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_ElemC, _Rows, _Cols>& C);

    // A - scalar = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // scalar - A = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void Sub(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // A .* B = C
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_ElemC, _Rows, _Cols>& C);

    // A .* scalar = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // A ./ B = C
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_ElemC, _Rows, _Cols>& C);

    // A ./ scalar = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // scalar ./ A = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void DivEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_ElemC, _Rows, _Cols>& B);

    //  A .pow B = C
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_ElemC, _Rows, _Cols>& C);

    // A .pow scalar = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB = void, typename _ElemC>
    static void PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // scalar .pow A = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
    static void PowEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_ElemC, _Rows, _Cols>& B);

    // A .exp = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
    static void ExpEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B);

    // A .sqrt = B
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
    static void SqrtEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B);

    // Transpose matrix 90 degrees to the right
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
    static void Transpose90deg(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Cols, _Rows>& B);

    // Transpose Matrix 180 degrees
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
    static void Transpose180deg(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B);

    // Transpose matrix 270 degrees to the right, or 90 degrees to the left
    template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
    static void Transpose270deg(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Cols, _Rows>& B);
}

#undef _INLINE

#ifdef _MATRIX_INLINE_TMP
#define _INLINE _MATRIX_INLINE_TMP
#undef _MATRIX_INLINE_TMP
#endif // _MATRIX_INLINE_TMP


//
///
//// Function Definitions:
///
//


template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemB>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::Matrix2D(const _ElemB init) :
    p_mat(new _Elem[_Rows * _Cols])
{
    for (size_t i = 0; i < c_elements; i++)
        p_mat[i] = _Elem(init);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::Matrix2D(const Matrix2D& _initilizer) :
    p_mat(new _Elem[_Rows * _Cols])
{
    memcpy(p_mat, _initilizer.p_mat, _Cols * _Rows * sizeof(_Elem));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemB>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::Matrix2D(Matrix2D<_ElemB, _Rows, _Cols>& _initilizer) :
    p_mat(new _Elem[_Rows * _Cols])
{
    for (size_t i = 0; i < _Cols * _Rows; i++)
        p_mat[i] = _Elem(_initilizer[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::Matrix2D(const std::initializer_list<std::initializer_list<_Elem>> lists) :
    p_mat(new _Elem[_Rows * _Cols])
{
    size_t i = 0;
    for (const auto list : lists)
        for (const auto element : list)
        {
            if (i < c_elements) p_mat[i++] = element;
            else i++;
        }

    if (i != c_elements) std::cout << "ERROR! Initializing matrix with " << i << " elements, " << c_elements << " expected!" << std::endl;
    if (i < c_elements) memset(p_mat + i, 0, (c_elements - i) * sizeof(_Elem));
}

#undef _MATRIX_P_MAT_ASSIGNMENT

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::Matrix2D(_Elem* ptr) :
    p_mat(ptr)
{}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>::~Matrix2D(void)
{
    delete[] p_mat;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _Cast>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::cout(void)
{
    // calculate the longest element per [_Cols] (column)
    // generate a string with all elements, separated by a space [' ']
    std::stringstream ss;
    for (size_t i = 0; i < (_Cols * _Rows); i++)
        ss << _Cast(p_mat[i]) << ' ';
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector <std::string> mat_vals(begin, end);

    // find out the longest string per [_Cols] 
    size_t max_length[_Cols] = { 0 };
    for (size_t row = 0; row < _Rows; row++)
        for (size_t column = 0; column < _Cols; column++)
            max_length[column] = std::max(max_length[column], mat_vals[POS_XY(row, column)].size());


    // increase all by one to add 1 space of padding
    for (auto& element : max_length)
        element += 1;

    // print each element with the needed padding corresponding to the longest
    // element in the same [_Cols]
    for (size_t i = 0; i < (_Cols * _Rows); i++)
    {
        if ((i % _Cols) == 0) std::cout << '|';
        std::cout << std::setw(max_length[i % _Cols]) << _Cast(p_mat[i]);
        if ((i % _Cols) == _Cols - 1) std::cout << ' ' << '|' << std::endl;
        else std::cout << ',';
    }
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _Elm2, size_t _Rows2, size_t _Cols2>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator *= (Matrix2D<_Elm2, _Rows2, _Cols2>& _mat2)
{
    static_assert(_Cols == _Rows2, "MATRIX ASSERT: Cannot multiply matrices A and B where A.Colums != B.Rows");
    static_assert(_Cols == _Cols2, "MATRIX ASSERT: The destination matrix must have the specific dimension if A.Rows x B.Columns!");

    // temporary container needed, because the Mul()-function will malfunction if two inputs are the same object
    Matrix2D tmp = *this;
    Matrix::Mul<_Elem, _Rows, _Cols, _Elm2, _Cols2, _Elem>(tmp, _mat2, *this);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename T, size_t rows, size_t columns>
static void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Print(Matrix2D<T, rows, columns>& mat)
{
    mat.cout();
}

template<typename _Elem, size_t _Rows, size_t _Cols>
bool Matrix::Matrix2D<_Elem, _Rows, _Cols>::SameDimensions(const Matrix2D& _mat)
{
    return ((_Cols == _mat.c_cols) and (_Rows == _mat.c_rows));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>& Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator = (const Matrix2D& _mat2)
{
    if (p_mat == NULL) p_mat = new _Elem[_Cols * _Rows];

    memcpy(p_mat, _mat2.p_mat, _Cols * _Rows * sizeof(_Elem));
    return *this;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>& Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator = (_Elem* ptr)
{
    if (p_mat) delete[] ptr;
    p_mat = ptr;
    return *this;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
Matrix::Matrix2D<_Elem, _Rows, _Cols>& Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator = (const _Elem val)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        p_mat[i] = val;
    return *this;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
_Elem* Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator + (Matrix2D& _mat2)
{
    _Elem* tmp = new _Elem[_Cols * _Rows];

    for (size_t i = 0; i < (_Cols * _Rows); i++)
        tmp[i] = this->at(i) + _mat2[i];

    return tmp;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
_Elem* Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator - (Matrix2D& _mat2)
{
    _Elem* tmp = new _Elem[_Cols * _Rows];

    for (size_t i = 0; i < (_Cols * _Rows); i++)
        tmp[i] = this->at(i) - _mat2[i];

    return tmp;
}

template<typename _Elem, size_t _Rows, size_t _Cols>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator += (Matrix2D& _mat2)
{
    for (size_t i = 0; i < (_Cols * _Rows); i++)
        p_mat[i] += _mat2[i];
}

template<typename _Elem, size_t _Rows, size_t _Cols>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator -= (Matrix2D& _mat2)
{
    for (size_t i = 0; i < (_Cols * _Rows); i++)
        p_mat[i] -= _mat2[i];
}

template<typename _Elem, size_t _Rows, size_t _Cols>
bool Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator == (Matrix2D& _mat2)
{
    if (!SameDimensions(_mat2)) return false;

    return (memcmp(p_mat, _mat2.p_mat, _Cols * _Rows * sizeof(_Elem)) == 0);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
bool Matrix::Matrix2D<_Elem, _Rows, _Cols>::operator != (Matrix2D& _mat2)
{
    if (!SameDimensions(_mat2)) return true;

    return (bool)(memcmp(p_mat, _mat2.p_mat, _Cols * _Rows * sizeof(_Elem)));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, size_t _RowsA, size_t _ColsA, typename _ElemB, size_t _ColsB, typename _ElemC>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Mul(Matrix2D<_ElemA, _RowsA, _ColsA>& A, Matrix2D<_ElemB, _ColsA, _ColsB>& B, Matrix2D<_ElemC, _RowsA, _ColsB>& C)
{
    static_assert(_Cols == _ColsA, "MATRIX ASSERT: Cannot multiply matrices A and B where A.Colums != B.Rows");
    static_assert(_Cols == _ColsB, "MATRIX ASSERT: The destination matrix must have the specific dimension if A.Rows x B.Columns!");

    Matrix::Mul<_ElemA, _RowsA, _ColsA, _ElemB, _ColsB, _ElemC>(A, B, C);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Add(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _Elem(A[i] + B[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Add(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(A[i] + scalar);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _Elem(A[i] - B[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Sub(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(A[i] - scalar);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::Sub(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(scalar - A[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _Elem(A[i] * B[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::MulEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(A[i] * scalar);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _Elem(A[i] / B[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::DivEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(A[i] / scalar);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::DivEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(scalar / A[i]);
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_ElemB, _Rows, _Cols>& B, Matrix2D<_Elem, _Rows, _Cols>& C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _Elem(std::pow(A[i], B[i]));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::PowEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, const _ElemB scalar, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(std::pow(A[i], scalar));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA, typename _ElemB>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::PowEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(std::pow(scalar, A[i]));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::ExpEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(std::exp(A[i]));
}

template<typename _Elem, size_t _Rows, size_t _Cols>
template<typename _ElemA>
void Matrix::Matrix2D<_Elem, _Rows, _Cols>::SqrtEBE(Matrix2D<_ElemA, _Rows, _Cols>& A, Matrix2D<_Elem, _Rows, _Cols>& B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _Elem(std::sqrt(A[i]));
}

template<typename _ElemA, size_t _RowsA, size_t _ColsA, typename _ElemB, size_t _ColsB, typename _ElemC>
void Matrix::Mul(Matrix2D<_ElemA, _RowsA, _ColsA>& A, Matrix2D<_ElemB, _ColsA, _ColsB>& B, Matrix2D<_ElemC, _RowsA, _ColsB>& C)
{
    constexpr auto _RowsB = _ColsA;
    constexpr auto _RowsC = _RowsA;
    constexpr auto _ColsC = _ColsB;

    // test if an input-Matrix is the same as the output-Matrix, because this would be bad
    if (_RowsA == _RowsC and _ColsA == _ColsC)
    {
        if ((void*)(&A) == (void*)(&C)) std::cerr << "MATRIX WARNING: You are trying to multiply two Matrices "
            "and input A is the same object as output C. If and input is the same as "
            "an output, the multiplication may fail and give an invalid result!" << std::endl;

        if ((void*)(&B) == (void*)(&C)) std::cerr << "MATRIX WARNING: You are trying to multiply two Matrices "
            "and input B is the same object as output C. If and input is the same as "
            "an output, the multiplication may fail and give an invalid result!" << std::endl;
    }

    /*
    A quick reminder on how to multiply matrices

    1. A.Cols == B.Rows (must equal)
    2. The final matrix will have the dimensions A.Rows x B.Cols

    Matrix A: (2 x 3)
    -         -
    | a, b, c |
    | d, e, f |
    -         -

    Matrix B: (3 x 3)
    -         -
    | g, h, i |
    | j, k, l |
    | m, n, o |
    -         -

    Matrix C: ( 2 * 3)
    -                                    -
    | (ag+bi+ck), (ah+bk+cn), (ai+bl+co) |
    | (dg+ei+fk), (dh+ek+fn), (di+el+fo) |
    -                                    -

    */

    // Full amount of functions for Release mode
#if defined(_MATRIX_MANUAL_OPTIMISATION) && defined(NDEBUG) && (_MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N > 4)
#pragma message  ( "Compiling Matrix2D unrolled multiplication function for up to 5x5" )
    // using multiple if-statements to detect common Matrix sizes. The compiler will optimize all
    // necessary if-statements for each matrix type away, since the arguments are known at compile time.

    // 1x1 * 1xn
    if (true)
    {
        // 1x1 * 1x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            return;
        }

        // 1x1 * 1x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            return;
        }

        // 1x1 * 1x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            return;
        }

        // 1x1 * 1x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            return;
        }

        // 1x1 * 1x5 -> 1x5
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)));

            return;
        }

    }

    // 2x1 * 1xn
    if (true)
    {
        // 2x1 * 1x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            return;
        }

        // 2x1 * 1x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            return;
        }

        // 2x1 * 1x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            return;
        }

        // 2x1 * 1x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            return;
        }

        // 2x1 * 1x5 -> 2x5
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)));

            return;
        }

    }

    // 3x1 * 1xn
    if (true)
    {
        // 3x1 * 1x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            return;
        }

        // 3x1 * 1x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            return;
        }

        // 3x1 * 1x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            return;
        }

        // 3x1 * 1x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));

            return;
        }

        // 3x1 * 1x5 -> 3x5
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)));

            return;
        }

    }

    // 4x1 * 1xn
    if (true)
    {
        // 4x1 * 1x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));

            return;
        }

        // 4x1 * 1x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));

            return;
        }

        // 4x1 * 1x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));

            return;
        }

        // 4x1 * 1x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)));

            return;
        }

        // 4x1 * 1x5 -> 4x5
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)));

            return;
        }

    }

    // 5x1 * 1xn
    if (true)
    {
        // 5x1 * 1x1 -> 5x1
        if (_RowsA == 5 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)));

            return;
        }

        // 5x1 * 1x2 -> 5x2
        if (_RowsA == 5 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)));

            return;
        }

        // 5x1 * 1x3 -> 5x3
        if (_RowsA == 5 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)));

            return;
        }

        // 5x1 * 1x4 -> 5x4
        if (_RowsA == 5 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)));

            return;
        }

        // 5x1 * 1x5 -> 5x5
        if (_RowsA == 5 and _ColsA == 1 and _RowsB == 1 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)));
            C.at(4, 4) = _ElemC((A.at(4, 0) * B.at(0, 4)));

            return;
        }

    }

    // Matrices n1x1 * 1xn2
    if (_ColsA == 1 and _RowsB == 1)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)));
        return;
    }

    // 1x2 * 2xn
    if (true)
    {
        // 1x2 * 2x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            return;
        }

        // 1x2 * 2x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            return;
        }

        // 1x2 * 2x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            return;
        }

        // 1x2 * 2x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            return;
        }

        // 1x2 * 2x5 -> 1x5
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)));

            return;
        }

    }

    // 2x2 * 2xn
    if (true)
    {
        // 2x2 * 2x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            return;
        }

        // 2x2 * 2x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            return;
        }

        // 2x2 * 2x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            return;
        }

        // 2x2 * 2x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            return;
        }

        // 2x2 * 2x5 -> 2x5
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)));

            return;
        }

    }

    // 3x2 * 2xn
    if (true)
    {
        // 3x2 * 2x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            return;
        }

        // 3x2 * 2x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            return;
        }

        // 3x2 * 2x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            return;
        }

        // 3x2 * 2x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));

            return;
        }

        // 3x2 * 2x5 -> 3x5
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)));

            return;
        }

    }

    // 4x2 * 2xn
    if (true)
    {
        // 4x2 * 2x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));

            return;
        }

        // 4x2 * 2x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));

            return;
        }

        // 4x2 * 2x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));

            return;
        }

        // 4x2 * 2x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)));

            return;
        }

        // 4x2 * 2x5 -> 4x5
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)));

            return;
        }

    }

    // 5x2 * 2xn
    if (true)
    {
        // 5x2 * 2x1 -> 5x1
        if (_RowsA == 5 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)));

            return;
        }

        // 5x2 * 2x2 -> 5x2
        if (_RowsA == 5 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)));

            return;
        }

        // 5x2 * 2x3 -> 5x3
        if (_RowsA == 5 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)));

            return;
        }

        // 5x2 * 2x4 -> 5x4
        if (_RowsA == 5 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)));

            return;
        }

        // 5x2 * 2x5 -> 5x5
        if (_RowsA == 5 and _ColsA == 2 and _RowsB == 2 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)));
            C.at(4, 4) = _ElemC((A.at(4, 0) * B.at(0, 4)) + (A.at(4, 1) * B.at(1, 4)));

            return;
        }

    }

    // Matrices n1x2 * 2xn2
    if (_ColsA == 2 and _RowsB == 2)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)));
        return;
    }

    // 1x3 * 3xn
    if (true)
    {
        // 1x3 * 3x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            return;
        }

        // 1x3 * 3x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            return;
        }

        // 1x3 * 3x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            return;
        }

        // 1x3 * 3x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            return;
        }

        // 1x3 * 3x5 -> 1x5
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)));

            return;
        }

    }

    // 2x3 * 3xn
    if (true)
    {
        // 2x3 * 3x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            return;
        }

        // 2x3 * 3x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            return;
        }

        // 2x3 * 3x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            return;
        }

        // 2x3 * 3x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            return;
        }

        // 2x3 * 3x5 -> 2x5
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)));

            return;
        }

    }

    // 3x3 * 3xn
    if (true)
    {
        // 3x3 * 3x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            return;
        }

        // 3x3 * 3x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            return;
        }

        // 3x3 * 3x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            return;
        }

        // 3x3 * 3x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));

            return;
        }

        // 3x3 * 3x5 -> 3x5
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)));

            return;
        }

    }

    // 4x3 * 3xn
    if (true)
    {
        // 4x3 * 3x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));

            return;
        }

        // 4x3 * 3x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));

            return;
        }

        // 4x3 * 3x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));

            return;
        }

        // 4x3 * 3x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)));

            return;
        }

        // 4x3 * 3x5 -> 4x5
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)));

            return;
        }

    }

    // 5x3 * 3xn
    if (true)
    {
        // 5x3 * 3x1 -> 5x1
        if (_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)));

            return;
        }

        // 5x3 * 3x2 -> 5x2
        if (_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)));

            return;
        }

        // 5x3 * 3x3 -> 5x3
        if (_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)));

            return;
        }

        // 5x3 * 3x4 -> 5x4
        if (_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)));

            return;
        }

        // 5x3 * 3x5 -> 5x5
        if (_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)));
            C.at(4, 4) = _ElemC((A.at(4, 0) * B.at(0, 4)) + (A.at(4, 1) * B.at(1, 4)) + (A.at(4, 2) * B.at(2, 4)));

            return;
        }

    }

    // Matrices n1x3 * 3xn2
    if (_ColsA == 3 and _RowsB == 3)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)));
        return;
    }

    // 1x4 * 4xn
    if (true)
    {
        // 1x4 * 4x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            return;
        }

        // 1x4 * 4x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            return;
        }

        // 1x4 * 4x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            return;
        }

        // 1x4 * 4x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            return;
        }

        // 1x4 * 4x5 -> 1x5
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)));

            return;
        }

    }

    // 2x4 * 4xn
    if (true)
    {
        // 2x4 * 4x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            return;
        }

        // 2x4 * 4x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            return;
        }

        // 2x4 * 4x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            return;
        }

        // 2x4 * 4x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            return;
        }

        // 2x4 * 4x5 -> 2x5
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)));

            return;
        }

    }

    // 3x4 * 4xn
    if (true)
    {
        // 3x4 * 4x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));

            return;
        }

        // 3x4 * 4x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));

            return;
        }

        // 3x4 * 4x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));

            return;
        }

        // 3x4 * 4x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));

            return;
        }

        // 3x4 * 4x5 -> 3x5
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)));

            return;
        }

    }

    // 4x4 * 4xn
    if (true)
    {
        // 4x4 * 4x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));

            return;
        }

        // 4x4 * 4x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));

            return;
        }

        // 4x4 * 4x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));

            return;
        }

        // 4x4 * 4x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)));

            return;
        }

        // 4x4 * 4x5 -> 4x5
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)) + (A.at(3, 3) * B.at(3, 4)));

            return;
        }

    }

    // 5x4 * 4xn
    if (true)
    {
        // 5x4 * 4x1 -> 5x1
        if (_RowsA == 5 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)));

            return;
        }

        // 5x4 * 4x2 -> 5x2
        if (_RowsA == 5 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)));

            return;
        }

        // 5x4 * 4x3 -> 5x3
        if (_RowsA == 5 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)));

            return;
        }

        // 5x4 * 4x4 -> 5x4
        if (_RowsA == 5 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)) + (A.at(4, 3) * B.at(3, 3)));

            return;
        }

        // 5x4 * 4x5 -> 5x5
        if (_RowsA == 5 and _ColsA == 4 and _RowsB == 4 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)) + (A.at(3, 3) * B.at(3, 4)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)) + (A.at(4, 3) * B.at(3, 3)));
            C.at(4, 4) = _ElemC((A.at(4, 0) * B.at(0, 4)) + (A.at(4, 1) * B.at(1, 4)) + (A.at(4, 2) * B.at(2, 4)) + (A.at(4, 3) * B.at(3, 4)));

            return;
        }

    }

    // Matrices n1x4 * 4xn2
    if (_ColsA == 4 and _RowsB == 4)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)) + (A.at(rows, 3) * B.at(3, cols)));
        return;
    }

    // 1x5 * 5xn
    if (true)
    {
        // 1x5 * 5x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 5 and _RowsB == 5 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));

            return;
        }

        // 1x5 * 5x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 5 and _RowsB == 5 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));

            return;
        }

        // 1x5 * 5x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 5 and _RowsB == 5 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));

            return;
        }

        // 1x5 * 5x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 5 and _RowsB == 5 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));

            return;
        }

        // 1x5 * 5x5 -> 1x5
        if (_RowsA == 1 and _ColsA == 5 and _RowsB == 5 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)) + (A.at(0, 4) * B.at(4, 4)));

            return;
        }

    }

    // 2x5 * 5xn
    if (true)
    {
        // 2x5 * 5x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 5 and _RowsB == 5 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));

            return;
        }

        // 2x5 * 5x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 5 and _RowsB == 5 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));

            return;
        }

        // 2x5 * 5x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 5 and _RowsB == 5 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));

            return;
        }

        // 2x5 * 5x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 5 and _RowsB == 5 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));

            return;
        }

        // 2x5 * 5x5 -> 2x5
        if (_RowsA == 2 and _ColsA == 5 and _RowsB == 5 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)) + (A.at(0, 4) * B.at(4, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)) + (A.at(1, 4) * B.at(4, 4)));

            return;
        }

    }

    // 3x5 * 5xn
    if (true)
    {
        // 3x5 * 5x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 5 and _RowsB == 5 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));

            return;
        }

        // 3x5 * 5x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 5 and _RowsB == 5 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));

            return;
        }

        // 3x5 * 5x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 5 and _RowsB == 5 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));

            return;
        }

        // 3x5 * 5x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 5 and _RowsB == 5 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));

            return;
        }

        // 3x5 * 5x5 -> 3x5
        if (_RowsA == 3 and _ColsA == 5 and _RowsB == 5 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)) + (A.at(0, 4) * B.at(4, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)) + (A.at(1, 4) * B.at(4, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)) + (A.at(2, 4) * B.at(4, 4)));

            return;
        }

    }

    // 4x5 * 5xn
    if (true)
    {
        // 4x5 * 5x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 5 and _RowsB == 5 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));

            return;
        }

        // 4x5 * 5x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 5 and _RowsB == 5 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));

            return;
        }

        // 4x5 * 5x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 5 and _RowsB == 5 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));

            return;
        }

        // 4x5 * 5x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 5 and _RowsB == 5 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)) + (A.at(3, 4) * B.at(4, 3)));

            return;
        }

        // 4x5 * 5x5 -> 4x5
        if (_RowsA == 4 and _ColsA == 5 and _RowsB == 5 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)) + (A.at(0, 4) * B.at(4, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)) + (A.at(1, 4) * B.at(4, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)) + (A.at(2, 4) * B.at(4, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)) + (A.at(3, 4) * B.at(4, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)) + (A.at(3, 3) * B.at(3, 4)) + (A.at(3, 4) * B.at(4, 4)));

            return;
        }

    }

    // 5x5 * 5xn
    if (true)
    {
        // 5x5 * 5x1 -> 5x1
        if (_RowsA == 5 and _ColsA == 5 and _RowsB == 5 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)) + (A.at(4, 4) * B.at(4, 0)));

            return;
        }

        // 5x5 * 5x2 -> 5x2
        if (_RowsA == 5 and _ColsA == 5 and _RowsB == 5 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)) + (A.at(4, 4) * B.at(4, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)) + (A.at(4, 4) * B.at(4, 1)));

            return;
        }

        // 5x5 * 5x3 -> 5x3
        if (_RowsA == 5 and _ColsA == 5 and _RowsB == 5 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)) + (A.at(4, 4) * B.at(4, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)) + (A.at(4, 4) * B.at(4, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)) + (A.at(4, 4) * B.at(4, 2)));

            return;
        }

        // 5x5 * 5x4 -> 5x4
        if (_RowsA == 5 and _ColsA == 5 and _RowsB == 5 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)) + (A.at(3, 4) * B.at(4, 3)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)) + (A.at(4, 4) * B.at(4, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)) + (A.at(4, 4) * B.at(4, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)) + (A.at(4, 4) * B.at(4, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)) + (A.at(4, 3) * B.at(3, 3)) + (A.at(4, 4) * B.at(4, 3)));

            return;
        }

        // 5x5 * 5x5 -> 5x5
        if (_RowsA == 5 and _ColsA == 5 and _RowsB == 5 and _ColsB == 5)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)) + (A.at(0, 4) * B.at(4, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)) + (A.at(0, 4) * B.at(4, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)) + (A.at(0, 4) * B.at(4, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)) + (A.at(0, 4) * B.at(4, 3)));
            C.at(0, 4) = _ElemC((A.at(0, 0) * B.at(0, 4)) + (A.at(0, 1) * B.at(1, 4)) + (A.at(0, 2) * B.at(2, 4)) + (A.at(0, 3) * B.at(3, 4)) + (A.at(0, 4) * B.at(4, 4)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)) + (A.at(1, 4) * B.at(4, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)) + (A.at(1, 4) * B.at(4, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)) + (A.at(1, 4) * B.at(4, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)) + (A.at(1, 4) * B.at(4, 3)));
            C.at(1, 4) = _ElemC((A.at(1, 0) * B.at(0, 4)) + (A.at(1, 1) * B.at(1, 4)) + (A.at(1, 2) * B.at(2, 4)) + (A.at(1, 3) * B.at(3, 4)) + (A.at(1, 4) * B.at(4, 4)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)) + (A.at(2, 4) * B.at(4, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)) + (A.at(2, 4) * B.at(4, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)) + (A.at(2, 4) * B.at(4, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)) + (A.at(2, 4) * B.at(4, 3)));
            C.at(2, 4) = _ElemC((A.at(2, 0) * B.at(0, 4)) + (A.at(2, 1) * B.at(1, 4)) + (A.at(2, 2) * B.at(2, 4)) + (A.at(2, 3) * B.at(3, 4)) + (A.at(2, 4) * B.at(4, 4)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)) + (A.at(3, 4) * B.at(4, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)) + (A.at(3, 4) * B.at(4, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)) + (A.at(3, 4) * B.at(4, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)) + (A.at(3, 4) * B.at(4, 3)));
            C.at(3, 4) = _ElemC((A.at(3, 0) * B.at(0, 4)) + (A.at(3, 1) * B.at(1, 4)) + (A.at(3, 2) * B.at(2, 4)) + (A.at(3, 3) * B.at(3, 4)) + (A.at(3, 4) * B.at(4, 4)));

            C.at(4, 0) = _ElemC((A.at(4, 0) * B.at(0, 0)) + (A.at(4, 1) * B.at(1, 0)) + (A.at(4, 2) * B.at(2, 0)) + (A.at(4, 3) * B.at(3, 0)) + (A.at(4, 4) * B.at(4, 0)));
            C.at(4, 1) = _ElemC((A.at(4, 0) * B.at(0, 1)) + (A.at(4, 1) * B.at(1, 1)) + (A.at(4, 2) * B.at(2, 1)) + (A.at(4, 3) * B.at(3, 1)) + (A.at(4, 4) * B.at(4, 1)));
            C.at(4, 2) = _ElemC((A.at(4, 0) * B.at(0, 2)) + (A.at(4, 1) * B.at(1, 2)) + (A.at(4, 2) * B.at(2, 2)) + (A.at(4, 3) * B.at(3, 2)) + (A.at(4, 4) * B.at(4, 2)));
            C.at(4, 3) = _ElemC((A.at(4, 0) * B.at(0, 3)) + (A.at(4, 1) * B.at(1, 3)) + (A.at(4, 2) * B.at(2, 3)) + (A.at(4, 3) * B.at(3, 3)) + (A.at(4, 4) * B.at(4, 3)));
            C.at(4, 4) = _ElemC((A.at(4, 0) * B.at(0, 4)) + (A.at(4, 1) * B.at(1, 4)) + (A.at(4, 2) * B.at(2, 4)) + (A.at(4, 3) * B.at(3, 4)) + (A.at(4, 4) * B.at(4, 4)));

            return;
        }

    }

    // Matrices n1x5 * 5xn2
    if (_ColsA == 5 and _RowsB == 5)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)) + (A.at(rows, 3) * B.at(3, cols)) + (A.at(rows, 4) * B.at(4, cols)));
        return;
    }
#elif defined(_MATRIX_MANUAL_OPTIMISATION) && defined(NDEBUG) && (_MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N > 3)
#pragma message ("Compiling Matrix2D unrolled multiplication function for up to 4x4")
    if (true)
    {
        // 1x1 * 1x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            return;
        }

        // 1x1 * 1x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            return;
        }

        // 1x1 * 1x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            return;
        }

        // 1x1 * 1x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            return;
        }

    }

    // 2x1 * 1xn
    if (true)
    {
        // 2x1 * 1x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            return;
        }

        // 2x1 * 1x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            return;
        }

        // 2x1 * 1x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            return;
        }

        // 2x1 * 1x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            return;
        }

    }

    // 3x1 * 1xn
    if (true)
    {
        // 3x1 * 1x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            return;
        }

        // 3x1 * 1x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            return;
        }

        // 3x1 * 1x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            return;
        }

        // 3x1 * 1x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));

            return;
        }

    }

    // 4x1 * 1xn
    if (true)
    {
        // 4x1 * 1x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));

            return;
        }

        // 4x1 * 1x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));

            return;
        }

        // 4x1 * 1x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));

            return;
        }

        // 4x1 * 1x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 1 and _RowsB == 1 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)));

            return;
        }

    }

    // Matrices n1x1 * 1xn2
    if (_ColsA == 1 and _RowsB == 1)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)));
        return;
    }

    // 1x2 * 2xn
    if (true)
    {
        // 1x2 * 2x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            return;
        }

        // 1x2 * 2x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            return;
        }

        // 1x2 * 2x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            return;
        }

        // 1x2 * 2x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            return;
        }

    }

    // 2x2 * 2xn
    if (true)
    {
        // 2x2 * 2x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            return;
        }

        // 2x2 * 2x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            return;
        }

        // 2x2 * 2x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            return;
        }

        // 2x2 * 2x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            return;
        }

    }

    // 3x2 * 2xn
    if (true)
    {
        // 3x2 * 2x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            return;
        }

        // 3x2 * 2x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            return;
        }

        // 3x2 * 2x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            return;
        }

        // 3x2 * 2x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));

            return;
        }

    }

    // 4x2 * 2xn
    if (true)
    {
        // 4x2 * 2x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));

            return;
        }

        // 4x2 * 2x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));

            return;
        }

        // 4x2 * 2x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));

            return;
        }

        // 4x2 * 2x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 2 and _RowsB == 2 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)));

            return;
        }

    }

    // Matrices n1x2 * 2xn2
    if (_ColsA == 2 and _RowsB == 2)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)));
        return;
    }

    // 1x3 * 3xn
    if (true)
    {
        // 1x3 * 3x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            return;
        }

        // 1x3 * 3x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            return;
        }

        // 1x3 * 3x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            return;
        }

        // 1x3 * 3x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            return;
        }

    }

    // 2x3 * 3xn
    if (true)
    {
        // 2x3 * 3x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            return;
        }

        // 2x3 * 3x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            return;
        }

        // 2x3 * 3x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            return;
        }

        // 2x3 * 3x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            return;
        }

    }

    // 3x3 * 3xn
    if (true)
    {
        // 3x3 * 3x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            return;
        }

        // 3x3 * 3x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            return;
        }

        // 3x3 * 3x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            return;
        }

        // 3x3 * 3x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));

            return;
        }

    }

    // 4x3 * 3xn
    if (true)
    {
        // 4x3 * 3x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));

            return;
        }

        // 4x3 * 3x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));

            return;
        }

        // 4x3 * 3x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));

            return;
        }

        // 4x3 * 3x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 3 and _RowsB == 3 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)));

            return;
        }

    }

    // Matrices n1x3 * 3xn2
    if (_ColsA == 3 and _RowsB == 3)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)));
        return;
    }

    // 1x4 * 4xn
    if (true)
    {
        // 1x4 * 4x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            return;
        }

        // 1x4 * 4x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            return;
        }

        // 1x4 * 4x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            return;
        }

        // 1x4 * 4x4 -> 1x4
        if (_RowsA == 1 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            return;
        }

    }

    // 2x4 * 4xn
    if (true)
    {
        // 2x4 * 4x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            return;
        }

        // 2x4 * 4x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            return;
        }

        // 2x4 * 4x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            return;
        }

        // 2x4 * 4x4 -> 2x4
        if (_RowsA == 2 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            return;
        }

    }

    // 3x4 * 4xn
    if (true)
    {
        // 3x4 * 4x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));

            return;
        }

        // 3x4 * 4x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));

            return;
        }

        // 3x4 * 4x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));

            return;
        }

        // 3x4 * 4x4 -> 3x4
        if (_RowsA == 3 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));

            return;
        }

    }

    // 4x4 * 4xn
    if (true)
    {
        // 4x4 * 4x1 -> 4x1
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));

            return;
        }

        // 4x4 * 4x2 -> 4x2
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));

            return;
        }

        // 4x4 * 4x3 -> 4x3
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));

            return;
        }

        // 4x4 * 4x4 -> 4x4
        if (_RowsA == 4 and _ColsA == 4 and _RowsB == 4 and _ColsB == 4)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)) + (A.at(0, 3) * B.at(3, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)) + (A.at(0, 3) * B.at(3, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)) + (A.at(0, 3) * B.at(3, 2)));
            C.at(0, 3) = _ElemC((A.at(0, 0) * B.at(0, 3)) + (A.at(0, 1) * B.at(1, 3)) + (A.at(0, 2) * B.at(2, 3)) + (A.at(0, 3) * B.at(3, 3)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)) + (A.at(1, 3) * B.at(3, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)) + (A.at(1, 3) * B.at(3, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)) + (A.at(1, 3) * B.at(3, 2)));
            C.at(1, 3) = _ElemC((A.at(1, 0) * B.at(0, 3)) + (A.at(1, 1) * B.at(1, 3)) + (A.at(1, 2) * B.at(2, 3)) + (A.at(1, 3) * B.at(3, 3)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)) + (A.at(2, 3) * B.at(3, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)) + (A.at(2, 3) * B.at(3, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)) + (A.at(2, 3) * B.at(3, 2)));
            C.at(2, 3) = _ElemC((A.at(2, 0) * B.at(0, 3)) + (A.at(2, 1) * B.at(1, 3)) + (A.at(2, 2) * B.at(2, 3)) + (A.at(2, 3) * B.at(3, 3)));

            C.at(3, 0) = _ElemC((A.at(3, 0) * B.at(0, 0)) + (A.at(3, 1) * B.at(1, 0)) + (A.at(3, 2) * B.at(2, 0)) + (A.at(3, 3) * B.at(3, 0)));
            C.at(3, 1) = _ElemC((A.at(3, 0) * B.at(0, 1)) + (A.at(3, 1) * B.at(1, 1)) + (A.at(3, 2) * B.at(2, 1)) + (A.at(3, 3) * B.at(3, 1)));
            C.at(3, 2) = _ElemC((A.at(3, 0) * B.at(0, 2)) + (A.at(3, 1) * B.at(1, 2)) + (A.at(3, 2) * B.at(2, 2)) + (A.at(3, 3) * B.at(3, 2)));
            C.at(3, 3) = _ElemC((A.at(3, 0) * B.at(0, 3)) + (A.at(3, 1) * B.at(1, 3)) + (A.at(3, 2) * B.at(2, 3)) + (A.at(3, 3) * B.at(3, 3)));

            return;
        }

    }

    // Matrices n1x4 * 4xn2
    if (_ColsA == 4 and _RowsB == 4)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)) + (A.at(rows, 3) * B.at(3, cols)));
        return;
    }

#elif (defined(_MATRIX_MANUAL_OPTIMISATION) && !defined(NDEBUG)) or (defined(_MATRIX_MANUAL_OPTIMISATION) && defined(NDEBUG) && (_MATRIX_LIMIT_MANUAL_OPTIMISATION_TO_N > 2)) 
#pragma message ("Compiling Matrix2D unrolled multiplication function for up to 3x3")
    // Smaller amount of functions for Debug Mode or explicit release mode

    // 1x1 * 1xn
    if (true)
    {
        // 1x1 * 1x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            return;
        }

        // 1x1 * 1x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            return;
        }

        // 1x1 * 1x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            return;
        }

    }

    // 2x1 * 1xn
    if (true)
    {
        // 2x1 * 1x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            return;
        }

        // 2x1 * 1x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            return;
        }

        // 2x1 * 1x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            return;
        }

    }

    // 3x1 * 1xn
    if (true)
    {
        // 3x1 * 1x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));

            return;
        }

        // 3x1 * 1x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));

            return;
        }

        // 3x1 * 1x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 1 and _RowsB == 1 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)));

            return;
        }

    }

    // Matrices n1x1 * 1xn2
    if (_ColsA == 1 and _RowsB == 1)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)));
        return;
    }

    // 1x2 * 2xn
    if (true)
    {
        // 1x2 * 2x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            return;
        }

        // 1x2 * 2x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            return;
        }

        // 1x2 * 2x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            return;
        }

    }

    // 2x2 * 2xn
    if (true)
    {
        // 2x2 * 2x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            return;
        }

        // 2x2 * 2x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            return;
        }

        // 2x2 * 2x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            return;
        }

    }

    // 3x2 * 2xn
    if (true)
    {
        // 3x2 * 2x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));

            return;
        }

        // 3x2 * 2x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));

            return;
        }

        // 3x2 * 2x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 2 and _RowsB == 2 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)));

            return;
        }

    }

    // Matrices n1x2 * 2xn2
    if (_ColsA == 2 and _RowsB == 2)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)));
        return;
    }

    // 1x3 * 3xn
    if (true)
    {
        // 1x3 * 3x1 -> 1x1
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            return;
        }

        // 1x3 * 3x2 -> 1x2
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            return;
        }

        // 1x3 * 3x3 -> 1x3
        if (_RowsA == 1 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            return;
        }

    }

    // 2x3 * 3xn
    if (true)
    {
        // 2x3 * 3x1 -> 2x1
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            return;
        }

        // 2x3 * 3x2 -> 2x2
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            return;
        }

        // 2x3 * 3x3 -> 2x3
        if (_RowsA == 2 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            return;
        }

    }

    // 3x3 * 3xn
    if (true)
    {
        // 3x3 * 3x1 -> 3x1
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 1)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));

            return;
        }

        // 3x3 * 3x2 -> 3x2
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));

            return;
        }

        // 3x3 * 3x3 -> 3x3
        if (_RowsA == 3 and _ColsA == 3 and _RowsB == 3 and _ColsB == 3)
        {
            C.at(0, 0) = _ElemC((A.at(0, 0) * B.at(0, 0)) + (A.at(0, 1) * B.at(1, 0)) + (A.at(0, 2) * B.at(2, 0)));
            C.at(0, 1) = _ElemC((A.at(0, 0) * B.at(0, 1)) + (A.at(0, 1) * B.at(1, 1)) + (A.at(0, 2) * B.at(2, 1)));
            C.at(0, 2) = _ElemC((A.at(0, 0) * B.at(0, 2)) + (A.at(0, 1) * B.at(1, 2)) + (A.at(0, 2) * B.at(2, 2)));

            C.at(1, 0) = _ElemC((A.at(1, 0) * B.at(0, 0)) + (A.at(1, 1) * B.at(1, 0)) + (A.at(1, 2) * B.at(2, 0)));
            C.at(1, 1) = _ElemC((A.at(1, 0) * B.at(0, 1)) + (A.at(1, 1) * B.at(1, 1)) + (A.at(1, 2) * B.at(2, 1)));
            C.at(1, 2) = _ElemC((A.at(1, 0) * B.at(0, 2)) + (A.at(1, 1) * B.at(1, 2)) + (A.at(1, 2) * B.at(2, 2)));

            C.at(2, 0) = _ElemC((A.at(2, 0) * B.at(0, 0)) + (A.at(2, 1) * B.at(1, 0)) + (A.at(2, 2) * B.at(2, 0)));
            C.at(2, 1) = _ElemC((A.at(2, 0) * B.at(0, 1)) + (A.at(2, 1) * B.at(1, 1)) + (A.at(2, 2) * B.at(2, 1)));
            C.at(2, 2) = _ElemC((A.at(2, 0) * B.at(0, 2)) + (A.at(2, 1) * B.at(1, 2)) + (A.at(2, 2) * B.at(2, 2)));

            return;
        }

    }

    // Matrices n1x3 * 3xn2
    if (_ColsA == 3 and _RowsB == 3)
    {
        for (size_t rows = 0; rows < _RowsC; rows++)
            for (size_t cols = 0; cols < _ColsC; cols++)
                C.at(rows, cols) = _ElemC((A.at(rows, 0) * B.at(0, cols)) + (A.at(rows, 1) * B.at(1, cols)) + (A.at(rows, 2) * B.at(2, cols)));
        return;
    }



#endif // _MATRIX_MANUAL_OPTIMISATION

    // If we don't use [_MATRIX_MANUAL_OPTIMISATION] or for matrices that don't fit
    // into the cases describes for [_MATRIX_MANUAL_OPTIMISATION]:

    // clear destination for preparation
    memset(&C[0], 0, _ColsC * _RowsC * sizeof(_ElemC));

    // Now multiply all and add to corresponding index
    for (size_t rowsA = 0, index = 0; rowsA < _RowsA; rowsA++)
        for (size_t colsB = 0; colsB < _ColsB; colsB++, index++)
            for (size_t rowCol = 0; rowCol < _ColsA; rowCol++)
                C[index] += _ElemC(A.at(rowsA, rowCol) * B.at(rowCol, colsB));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::Add(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B, Matrix2D<_ElemC, _Rows, _Cols>&C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _ElemC(A[i] + B[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::Add(Matrix2D<_ElemA, _Rows, _Cols>&A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(A[i] + scalar);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::Sub(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B, Matrix2D<_ElemC, _Rows, _Cols>&C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _ElemC(A[i] - B[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::Sub(Matrix2D<_ElemA, _Rows, _Cols>&A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(A[i] - scalar);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::Sub(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>&A, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(scalar - A[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::MulEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B, Matrix2D<_ElemC, _Rows, _Cols>&C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _ElemC(A[i] * B[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::MulEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(A[i] * scalar);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::DivEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B, Matrix2D<_ElemC, _Rows, _Cols>&C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _ElemC(A[i] / B[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::DivEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(A[i] / scalar);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::DivEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>&A, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(scalar / A[i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::PowEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B, Matrix2D<_ElemC, _Rows, _Cols>&C)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        C[i] = _ElemC(std::pow(A[i], B[i]));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::PowEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, const _ElemB scalar, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(std::pow(A[i], scalar));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB, typename _ElemC>
void Matrix::PowEBE(const _ElemA scalar, Matrix2D<_ElemB, _Rows, _Cols>&A, Matrix2D<_ElemC, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemC(std::pow(scalar, A[i]));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
void Matrix::ExpEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemB(std::exp(A[i]));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
void Matrix::SqrtEBE(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B)
{
    for (size_t i = 0; i < _Rows * _Cols; i++)
        B[i] = _ElemB(std::sqrt(A[i]));
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
void Matrix::Transpose90deg(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Cols, _Rows>&B)
{
    constexpr auto c_width = (_Rows * _Cols);

    constexpr auto _RowsA = _Rows;
    constexpr auto _ColsA = _Cols;
    constexpr auto _RowsB = _Cols;
    constexpr auto _ColsB = _Rows;
    constexpr auto a = _ColsA;
    constexpr auto b = _ColsB;

    for (size_t i = 0; i < c_width; i++)
        B[i] = _ElemB(A[(((a * b) - a) - (a * i)) + ((i / b) * (b * a)) + (i / b)]);

    //B = p_tmp;
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
void Matrix::Transpose180deg(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Rows, _Cols>&B)
{
    constexpr auto c_width = (_Rows * _Cols);

    for (size_t i = 0; i < c_width; i++)
        B[i] = _ElemB(A[(c_width - 1) - i]);
}

template<typename _ElemA, size_t _Rows, size_t _Cols, typename _ElemB>
void Matrix::Transpose270deg(Matrix2D<_ElemA, _Rows, _Cols>&A, Matrix2D<_ElemB, _Cols, _Rows>&B)
{
    constexpr auto c_width = (_Rows * _Cols);

    constexpr auto _RowsA = _Rows;
    constexpr auto _ColsA = _Cols;
    constexpr auto _RowsB = _Cols;
    constexpr auto _ColsB = _Rows;
    constexpr auto a = _ColsA;
    constexpr auto b = _ColsB;

    for (size_t i = 0; i < c_width; i++)
        B[i] = _ElemB(A[(i * a) - ((i / b) * (b * a + 1)) + (a - 1)]);
}
