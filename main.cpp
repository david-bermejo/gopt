#include <iostream>

#include "math/vector.hpp"
#include "math/matrix.hpp"
#include "math/algorithms.hpp"

int main()
{
	using namespace math;

	Vector<3> v(3, 3, 3);
	Matrix<3, 3> M(4, 12, -16, 12, 37, -43, -16, -43, 98);

	std::cout << v.toString() << std::endl;
}