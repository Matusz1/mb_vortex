#include <cassert>
#include "fock/Numeric.h"

namespace Mbs {

std::vector<double> linspace(double start, double end, uint n, bool endpoint)
{
    assert(n > 1 && "n < 2 in linspace algorithm\n");
    std::vector<double> res(n);
    const double h = (end - start) / (endpoint ? n - 1 : n);
    for (uint i = 0; i < n; ++i)
        res[i] = start + i * h;
    return res;
}

std::vector<double> logspace(double start, double end, uint n)
{
    assert(n > 1 && "n < 2 in logspace algorithm\n");
    assert(start > 0 && end > 0 && "start or end is not positive in logspace algorithm\n");
    start = std::log10(start);
    end = std::log10(end);
    std::vector<double> res = linspace(start, end, n);
    for (uint i = 0; i < n; ++i)
        res[i] = std::pow(10, res[i]);
    return res;
}

double simpson(double* f, uint n, double h)
{
    assert(n > 2 && "n < 3 in simpson integration algorithm\n");
    double s0 = 0;
    double s1 = 0;
    double s2 = 0;
    for (uint i = 1; i < n-1; i += 2) {
        s0 += f[i-1];
        s1 += f[i];
        s2 += f[i+1];
    }
    const double sum = (s0 + 4*s1 + s2) / 3;

    // If n is even, the interval is divided into odd numer of subintervals
    // and so the last interval is not included in the sum and we have to correct it
    if (IS_EVEN(n))
        return h * (sum + (5*f[n-1] + 8*f[n-2] - f[n-3]) / 12);
    else
        return h * sum;
}

} // namespace Mbs
