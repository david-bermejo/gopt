#include "math/gopt.hpp"
using namespace gopt;

#include <random>
#include <vector>

#if 0
double fcn(const Vector<double>& x)
{
    return x.magnitude();
}
#else
double fcn(const Vec2& x)
{
    return x.magnitude();
}
#endif

int main()
{
#if 0
    Vector<double> ub(2);
    Vector<double> lb(2);
    ub.fill(1);
    lb.fill(-1);

    const auto res = particleswarm(fcn, lb, ub, 2, 200, 100, {0.5, 0.5, 0.9});
    std::cout << res << std::endl;
#else
    Vec2 ub(1);
    Vec2 lb(-1);

    const auto res = particleswarm(fcn, lb, ub, 200, 100, {0.5, 0.5, 0.9, false});
    std::cout << res << std::endl;
#endif

    return 0;
}