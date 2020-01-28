#include <iostream>

#include "math/matrix.hpp"
#include "math/algorithms.hpp"

int main()
{
	using namespace math;

	Matrix_t<float, 3, 3> M(4, 12, -16, 12, 37, -43, -16, -43, 98);

	std::cout << M.toString() << std::endl;
}