#ifndef _PARTICLE_DRAWING_H_
#define _PARTICLE_DRAWING_H_

#include "FockState.h"
#include <limits>

namespace Mbs {

struct RPhi {
    double r;
    double phi;
};

struct XY {
    double x;
    double y;
};

inline XY polar_to_cartesian(double r, double phi)
{
    return {r * std::cos(phi), r * std::sin(phi)};
}

inline RPhi cartesian_to_polar(double x, double y)
{
    return {std::sqrt(x*x + y*y), std::atan2(y, x)};
}

std::vector<RPhi> draw_positions(const FockLinear<double>& state, double c_cutoff = std::numeric_limits<double>::min());
std::vector<Complex> fill_wave_function(const std::vector<RPhi>& positions, const FockLinear<double>& state, uint size = 101, double c_cutoff = std::numeric_limits<double>::min());

}

#endif // _PARTICLE_DRAWING_H_
