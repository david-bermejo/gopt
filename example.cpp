#include "math/gopt.hpp"
using namespace gopt;

double fcn(const Vector<double>& x)
{
    return -(4*x[0] -2*x[1] + 7*x[2] + 5*x[3] + 11*x[4] + x[5]);
}

int main()
{
    unsigned int N = 1000;
    Vector<double> lb(6);
    lb.fill(-10);
    Vector<double> ub(6);
    ub.fill(10);

    std::default_random_engine gen(time(0));
    std::uniform_real_distribution<double> dist(lb[0], ub[0]);
    std::vector<Vector<double>> x0(N);

    for (int i = 0; i < N; i++)
        x0[i] = Vector<double>(6).fill(dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen));
    
    auto [best, score] = ga(fcn, x0, lb, ub, 0.2, 0.8, 100.0, 100.0, 1000);
    std::cout << best << ", score: " << score << std::endl;

    return 0;
}