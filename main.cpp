#include <iostream>

#include "math/vector.hpp"
#include "math/matrix.hpp"
#include "math/algorithms.hpp"
#include "math/solvers.hpp"
#include "math/constants.hpp"

#include <chrono>
#include <iomanip>

using namespace gopt;

Vector<4> f(const Vector<4>& x, double t)
{
	return
	{
		x[2],
		x[3],
		-4*t*t*x[0] + 2*x[2] / (std::sqrt(x[0]*x[0] + x[1]*x[1]) * std::sqrt(x[2]*x[2] + x[3]*x[3])),
		-4*t*t*x[1] + 2*x[3] / (std::sqrt(x[0]*x[0] + x[1]*x[1]) * std::sqrt(x[2]*x[2] + x[3]*x[3])),
	};
}

Vector<2> f2(const Vector<2>& x, const Vector<2>& dx, double t)
{
	return
	{
		-4*t*t*x[0] + 2*dx[0] / (std::sqrt(x[0]*x[0] + x[1]*x[1]) * std::sqrt(dx[0]*dx[0] + dx[1]*dx[1])),
		-4*t*t*x[1] + 2*dx[1] / (std::sqrt(x[0]*x[0] + x[1]*x[1]) * std::sqrt(dx[0]*dx[0] + dx[1]*dx[1])),
	};
}

int main()
{
	Vector<4> x0 { 0, 1, -std::sqrt(2*pi<>), 0 };
	Vector<2> x00 { 0, 1 };
	Vector<2> dx00 {-std::sqrt(2 * pi<>), 0 };
	double t0 = std::sqrt(pi<>/2);

	auto start = std::chrono::system_clock::now();

	// TEST CODE
#if 1
	Vector<4> res;
	for (int i = 0; i < 500; i++)
		res = DP45(f, x0, t0, 10.0, 1.0e-9, 1.0e-9);
#else
	Vector<4> res;
	for (int i = 0; i < 500; i++)
		res = RKN5(f2, x00, dx00, t0, 10.0, 0.1e-16);
#endif

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<float, std::milli> duration = end - start;

	std::cout << "Analytical: {" << std::cos(100.0) << ", " << std::sin(100.0) << "}\n";
	std::cout << "Results: " << res << std::endl;
    std::cout << duration.count()/1000.0f << " ms" << std::endl;

	system("pause");
}