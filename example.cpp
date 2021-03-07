#include "math/gopt.hpp"
using namespace gopt;

#include <random>
#include <vector>

int main()
{
    Vector<double> v(5);
    v.fill([]() { return 122; });
    Vector<double> u(v);

    u = v;

    std::cout << v << ", " << u << std::endl;

    return 0;
}