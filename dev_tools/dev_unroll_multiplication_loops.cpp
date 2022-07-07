#include <iostream>

// generating unrolled loop from (1x1 * 1x1) to (nxn * nxn)
void GenerateInlineCodeMAT_MUL(const size_t max)
{
    constexpr auto tab = "    ";

    for (size_t i = 0; i < max; i++)
    {
        const auto _ColsA = (i + 1);
        const auto _RowsB = _ColsA;

        for (size_t k = 0; k < max; k++)
        {
            const auto _RowsA = (k + 1);

            std::cout << "// " << _RowsA << 'x' << _ColsA << " * " << _RowsB << "xn" << std::endl;
            std::cout << "if (true)" << std::endl;
            std::cout << '{' << std::endl;

            for (size_t c = 0; c < max; c++)
            {
                const auto _ColsB = (c + 1);

                const auto _RowsC = _RowsA;
                const auto _ColsC = _ColsB;

                std::cout << tab << "// " << _RowsA << 'x' << _ColsA << " * " << _RowsB << 'x' << _ColsB << " -> " << _RowsC << 'x' << _ColsC << std::endl;
                std::cout << tab << "if (_RowsA == " << _RowsA << " and _ColsA == " << _ColsA << " and _RowsB == " << _RowsB << " and _ColsB == " << _ColsB << ")" << std::endl;
                std::cout << tab << '{' << std::endl;
                //_RowsA == 5 and _ColsA == 3 and _RowsB == 3 and _ColsB == 2
                for (size_t rowsA = 0, index = 0; rowsA < _RowsA; rowsA++)
                {
                    for (size_t colsB = 0; colsB < _ColsB; colsB++, index++)
                    {
                        std::cout << tab << tab << "C.at(" << rowsA << ", " << colsB << ") = _ElemC(";
                        for (size_t rowCol = 0; rowCol < _ColsA; rowCol++)
                        {
                            //C[index] += _ElemC(A.at(rowsA, rowCol) * B.at(rowCol, colsB));
                            std::cout << "(A.at(" << rowsA << ", " << rowCol << ") * B.at(" << rowCol << ", " << colsB << "))";
                            if (rowCol < (_ColsA - 1)) std::cout << " + ";
                        }
                        std::cout << ");" << std::endl;
                    }
                    std::cout << std::endl;
                }
                std::cout << tab << tab << "return;" << std::endl;
                std::cout << tab << '}';

                std::cout << std::endl << std::endl;
            }

            std::cout << '}' << std::endl << std::endl;
        }

        std::cout << "// Matrices n1x" << _ColsA << " * " << _RowsB << "xn2" << std::endl;
        std::cout << "if (_ColsA == " << _ColsA << " and _RowsB == " << _RowsB << ')' << std::endl;
        std::cout << '{' << std::endl;

        std::cout << tab << "for (size_t rows = 0; rows < _RowsC; rows++)" << std::endl;
        std::cout << tab << tab << "for (size_t cols = 0; cols < _ColsC; cols++)" << std::endl;
        std::cout << tab << tab << tab << "C.at(rows, cols) = _ElemC(";


        for (size_t c = 0; c < _ColsA; c++)
        {
            std::cout << "(A.at(rows, " << c << ") * B.at(" << c << ", cols))";
            if (c < _ColsA - 1) std::cout << " + ";
        }

        std::cout << ");" << std::endl;
        std::cout << tab << "return;" << std::endl;

        std::cout << '}' << std::endl << std::endl;
    }
}
