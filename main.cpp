#include <iostream>

#include "math/vector.hpp"
#include "math/matrix.hpp"
#include "math/algorithms.hpp"
#include "math/solvers.hpp"
#include "math/constants.hpp"

#include <chrono>

using namespace gopt;

Vector<2> f(const Vector<2>& x, const Vector<2>& dx, double t)
{
	return
	{
		-4*t*t*x[0] + 2*dx[0]/(dx[0]*dx[0] + dx[1]*dx[1]),
		-4*t*t*x[1] + 2*dx[1]/(dx[0]*dx[0] + dx[1]*dx[1])
	};
}

int main()
{
	Vector<2> x0 { 0, 1 };
	Vector<2> dx0 { -std::sqrt(2*pi<>), 0 };
	double t0 = std::sqrt(pi<>/2);

	auto start = std::chrono::system_clock::now();
	// TEST CODE
	for (int i = 0; i < 300; i++)
		RKN5(f, x0, dx0, t0, 10.0, 0.1e-16);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<float, std::milli> duration = end - start;

    std::cout << duration.count()/1000.0f << " ms" << std::endl;
}