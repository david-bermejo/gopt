#include "math/gopt.hpp"
using namespace gopt;

#include <random>
#include <vector>

int main()
{
    Matrix<double> m(3, 4);

    std::default_random_engine gen;
    std::normal_distribution<double> dist(0.0, 1.0);
    m.fill([&]() { return dist(gen); });

    std::cout << m << std::endl;

    return 0;
}