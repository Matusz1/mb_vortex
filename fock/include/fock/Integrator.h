#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "Core.h"
#include <array>

namespace Mbs {

double potential_integrate(const std::array<u16,4>& idxs, uint eta);
double potential_integrate_npts(const std::array<u16,4>& idxs, uint eta, uint num_pts);

/* Cached integration functionality */

// This can be redefined to change the cache directory.
#define CONFIG_INT_CACHE_DIR "./Cache"

void integrator_setup(uint eta);
void integrator_cleanup(void);

/* The cached integrator, don't use it directly,
 * unless you know what you are doing */

double potential_integrate_cached_unchecked(std::array<u16,4> idxs);
double potential_integrate_cached(const std::array<u16,4>& idxs);

} // namespace Mbs

#endif // _INTEGRATOR_H_
