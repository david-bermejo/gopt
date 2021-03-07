#include "math/gopt.hpp"
using namespace gopt;

#include <vector>

int main()
{
    Vector<double> v(10);
    std::vector<double> src(10);
    for (int i = 0; i < 10; i++)
        src[i] = i,
    v.fill(src);

    std::cout << v << std::endl;

    return 0;
}