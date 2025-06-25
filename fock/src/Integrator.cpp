#include "fock/Integrator.h"
#include "fock/DiskState.h"
#include "fock/Numeric.h"

#include <gsl/gsl_sf_bessel.h>

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <map>
#include <unordered_map>

#define FOCK_INTEGRATOR_INT_PTS 255

/* Functions to directly integrate the gauss/delta function,
 * in practice those will be covered by caching the integrals,
 * such functionality is also provided */

namespace Mbs {

/* The caching functionality,
 * the data is saved to files described by 'eta' parameter */

using _CachePair = std::pair<std::array<u16,4>, double>;

// Thank you StackOverflow: https://stackoverflow.com/questions/42701688/using-an-unordered-map-with-arrays-as-keys
struct _ArrHash {
    std::size_t operator()(const std::array<u16,4>& arr) const
    {
        std::size_t h = 0;
        for (auto i : arr)
            h ^= std::hash<u16>{}(i) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct _IntCacheDat {
    int eta{0};
    std::string filename;
    std::unordered_map<std::array<u16,4>, double, _ArrHash> cache;
};

static struct _IntCacheDat _int_cache{};

// Some more caching but only for the wave functions, not saved to file
using _WfCache = std::array<double, FOCK_INTEGRATOR_INT_PTS>;
static std::map<u16, _WfCache> _wave_cache;

// Caching for the bessel functions, this is eta dependent so clear whenerver eta changes
// The cache is not saved to file
using _IBessCache = std::array<double, FOCK_INTEGRATOR_INT_PTS*FOCK_INTEGRATOR_INT_PTS>;
static std::map<uint, _IBessCache> _ibess_cache;

static const _WfCache& _get_wave_cache(u16 idx)
{
    auto it = _wave_cache.find(idx);
    if (it != _wave_cache.end())
        return it->second;

    _WfCache& cache = _wave_cache[idx];
    constexpr double h = 1.0 / (FOCK_INTEGRATOR_INT_PTS - 1.0);
    for (uint i = 0; i < FOCK_INTEGRATOR_INT_PTS; ++i)
        cache[i] = wave_function(idx, i*h);

    return cache;
}

static const _IBessCache& _get_ibess_cache(uint delta_m)
{
    auto it = _ibess_cache.find(delta_m);
    if (it != _ibess_cache.end())
        return it->second;

    constexpr double h = 1.0 / (FOCK_INTEGRATOR_INT_PTS - 1.0);
    const uint eta = _int_cache.eta;

    _IBessCache& cache = _ibess_cache[delta_m];
    for (uint i = 0; i < cache.size(); ++i)
        cache[i] = gsl_sf_bessel_In(delta_m, 2*eta*eta*i*h*h);
    return cache;
}

static double _int_delta(const std::array<u16, 4>& idxs)
{
    constexpr uint npts = FOCK_INTEGRATOR_INT_PTS;
    constexpr double h = 1.0 / (npts - 1);

    const _WfCache& wf0 = _get_wave_cache(idxs[0]);
    const _WfCache& wf1 = _get_wave_cache(idxs[1]);
    const _WfCache& wf2 = _get_wave_cache(idxs[2]);
    const _WfCache& wf3 = _get_wave_cache(idxs[3]);

    double f[npts];
    for (uint i = 0; i < npts; ++i)
        f[i] = wf0[i] * wf1[i] * wf2[i] * wf3[i] * i*h;

    return 2 * M_PI * simpson(f, npts, h);
}

static double _int_gauss(const std::array<u16,4>& idxs, uint eta)
{
    constexpr uint npts = FOCK_INTEGRATOR_INT_PTS;
    constexpr double h = 1.0 / (npts - 1);

    double f_inner[npts];
    double f_outer[npts];

    const _WfCache& wf0 = _get_wave_cache(idxs[0]);
    const _WfCache& wf1 = _get_wave_cache(idxs[1]);
    const _WfCache& wf2 = _get_wave_cache(idxs[2]);
    const _WfCache& wf3 = _get_wave_cache(idxs[3]);
    const _IBessCache& ibess = _get_ibess_cache(std::abs(angular_momentum(idxs[3]) - angular_momentum(idxs[1])));

    // Perform the double integral
    for (uint i = 0; i < npts; ++i) {
        f_outer[i] = i*h * wf0[i] * wf2[i];
        f_outer[i] *= std::exp(-(i*h) * (i*h) * eta * eta);
        for (uint j = 0; j < npts; ++j) {
            const double val = ibess[i*j] * std::exp(-(j*h) * (j*h) * eta * eta);
            f_inner[j] = j*h * val * wf1[j] * wf3[j];
        }
        f_outer[i] *= simpson(f_inner, npts, h);
    }

    return 4 * M_PI * eta * eta * simpson(f_outer, npts, h);
}

double potential_integrate(const std::array<u16,4>& idxs, uint eta)
{
    if (eta == 0)
        return _int_delta(idxs);
    return _int_gauss(idxs, eta);
}


static void _int_cache_read(void)
{
    FILE *f = fopen(_int_cache.filename.c_str(), "r");
    if (f == NULL)
        return;

    uint len;
    if (!fread(&len, sizeof(len), 1, f)) {
        fclose(f);
        return;
    }

    printf("Found cache file: %s (size = %u)\n", _int_cache.filename.c_str(), len);

    _int_cache.cache.reserve(len);
    for (uint i = 0; i < len; ++i) {
        _CachePair elem;
        if (!fread(&elem, sizeof(elem), 1, f)) {
            fclose(f);
            throw std::runtime_error("Failed to read cache file");
            return;
        }
        _int_cache.cache.insert(elem);
    }
    fclose(f);
}

static void _int_cache_write(void)
{
    if (_int_cache.cache.empty())
        return;
    FILE *f = fopen(_int_cache.filename.c_str(), "w");
    const uint size = _int_cache.cache.size();
    printf("Writing cache file: %s (size = %u)\n", _int_cache.filename.c_str(), size);

    fwrite(&size, sizeof(size), 1, f);
    for (const auto& elem : _int_cache.cache)
        fwrite(&elem, sizeof(elem), 1, f);
    fclose(f);
}

void integrator_setup(uint eta)
{
    integrator_cleanup();

    std::ostringstream oss;
    oss << CONFIG_INT_CACHE_DIR << "/eta_" << eta << ".dat";
    _int_cache.filename = oss.str();
    _int_cache.eta = eta;
    _int_cache_read();
}

void integrator_cleanup(void)
{
    if (!_int_cache.cache.empty()) {
        _int_cache_write();
        _int_cache.cache.clear();
    }
    _ibess_cache.clear(); // This cache is eta dependent, so clear it
    _int_cache.eta = 0;
    _int_cache.filename.clear();
}

static void _normalize_idx(std::array<u16,4>& idxs)
{
    if (std::max(idxs[0], idxs[1]) > std::max(idxs[2], idxs[3])) {
        std::swap(idxs[0], idxs[2]);
        std::swap(idxs[1], idxs[3]);
    }
    if (idxs[2] > idxs[3]) {
        std::swap(idxs[0], idxs[1]);
        std::swap(idxs[2], idxs[3]);
    }
}

double potential_integrate_cached_unchecked(std::array<u16,4> idxs)
{
    _normalize_idx(idxs);

    // Check if the integral is cached
    auto it = _int_cache.cache.find(idxs);
    if (it != _int_cache.cache.end())
        return it->second;

    // Not cached, evaluate the integral
    const double val = potential_integrate(idxs, _int_cache.eta);
    _int_cache.cache[idxs] = val;

    return val;
}

double potential_integrate_cached(const std::array<u16,4>& idxs)
{
    std::array<int,4> m = {
        angular_momentum(idxs[0]),
        angular_momentum(idxs[1]),
        angular_momentum(idxs[2]),
        angular_momentum(idxs[3]),
    };
    if (m[0] + m[1] != m[2] + m[3])
        return 0;

    return potential_integrate_cached_unchecked(idxs);
}

} // namespace Mbs
