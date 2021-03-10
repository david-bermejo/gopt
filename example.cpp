#include "math/gopt.hpp"
using namespace gopt;

int main()
{
    Matrix<double> m(4, 3);
    m.fill(4, 3, 2, 1, 6, 7, 8, 9, 0, 5, 4, 7);
    m(1) = Vector<double>(3).fill(3);
    m(3) = Vector<double>(3).fill(2);

    std::cout << m(3) << std::endl;
    std::cout << m << std::endl;

    return 0;
}