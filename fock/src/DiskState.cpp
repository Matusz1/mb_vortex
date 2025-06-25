#include "fock/DiskState.h"
#include "gsl/gsl_sf_bessel.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <vector>

namespace Mbs {

struct _StateInfo {
    int ang_mom;
    uint zero_num;
    double zero_val;
    double norm_const;

    _StateInfo(int m, uint k);
    [[nodiscard]] bool operator<(const _StateInfo& other) const;
    [[nodiscard]] bool operator>(const _StateInfo& other) const { return other < *this; }
};

inline _StateInfo::_StateInfo(int m, uint k) :
    ang_mom(m),
    zero_num(k),
    zero_val(gsl_sf_bessel_zero_Jnu(m, k))
{
    const double v = std::sqrt(M_PI) * gsl_sf_bessel_Jn(m+1, zero_val);
    norm_const = 1.0 / std::abs(v);
}

inline bool _StateInfo::operator<(const _StateInfo& other) const {
    if (zero_num == other.zero_num && ang_mom == -other.ang_mom)
        return ang_mom > 0;
    return zero_val < other.zero_val;
}

static std::vector<_StateInfo> _states_info;

static void _allocate_state_info(u16 max)
{
    std::priority_queue<_StateInfo, std::vector<_StateInfo>, std::greater<_StateInfo>> queue;
    queue.push(_StateInfo(0, 1));

    _states_info.clear();
    while (_states_info.size() < max) {
        _StateInfo elem = queue.top();
        queue.pop();

        queue.push(_StateInfo(elem.ang_mom+1, elem.zero_num));
        if (elem.ang_mom == 0)
            queue.push(_StateInfo(0, elem.zero_num+1));

        _states_info.push_back(elem);
        if (elem.ang_mom != 0) {
            elem.ang_mom = -elem.ang_mom;
            _states_info.push_back(elem);
        }
    }
}

static inline void _ensure_enough_state_info(u16 n)
{
    if (_states_info.size() <= n)
        _allocate_state_info(std::max(1, 2*n));
}

double energy(u16 n)
{
    _ensure_enough_state_info(n);
    const double k = _states_info[n].zero_val;
    return 0.5 * k * k;
}

Complex wave_function(u16 n, double r, double phi)
{
    _ensure_enough_state_info(n);
    const _StateInfo& info = _states_info[n];
    const double k = info.zero_val;
    const double rad = info.norm_const * gsl_sf_bessel_Jn(info.ang_mom, k*r);
    const double angle = info.ang_mom*phi;
    return Complex(rad * std::cos(angle), rad * std::sin(angle));
}

double wave_function(u16 n, double r)
{
    _ensure_enough_state_info(n);
    const _StateInfo& info = _states_info[n];
    const double k = info.zero_val;
    return info.norm_const * gsl_sf_bessel_Jn(info.ang_mom, k*r);
}

int angular_momentum(u16 n)
{
    _ensure_enough_state_info(n);
    return _states_info[n].ang_mom;
}

u16 stateIndex(int ang_mom, uint zero_num)
{
    assert(zero_num != 0 && "Cannot have 0-th 'zero of bessel funtion'.");
    u16 i = 0;
    while (true) {
        _ensure_enough_state_info(i);
        const _StateInfo& info = _states_info[i];
        if (info.ang_mom == ang_mom && info.zero_num == zero_num)
            return i;
        ++i;
    }
}

LinearCombination::LinearCombination(uint n, int mom) :
    m_coeff(n, 0.0),
    m_basis(n)
{
    u16 i = 0;
    for (uint j = 0; j < n; ++j) {
        while (angular_momentum(i) != mom)
            ++i;
        m_basis[j] = i;
        ++i;
    }
    m_coeff[0] = 1.0;
}

LinearCombination::LinearCombination(const double* coeff, const u16* basis, uint size) :
    m_coeff(coeff, coeff + size),
    m_basis(basis, basis + size)
{ }

void LinearCombination::normalize()
{
    double norm = 0.0;
    for (double c : m_coeff)
        norm += c*c;
    norm = std::sqrt(norm);
    for (double& c : m_coeff)
        c /= norm;
}

double LinearCombination::wave_function(double r) const
{
    double res = 0.0;
    for (uint i = 0; i < m_coeff.size(); ++i)
        res += m_coeff[i] * ::Mbs::wave_function(m_basis[i], r);
    return res;
}

Complex LinearCombination::wave_function(double r, double phi) const
{
    Complex res = 0.0;
    for (uint i = 0; i < m_coeff.size(); ++i)
        res += m_coeff[i] * ::Mbs::wave_function(m_basis[i], r, phi);
    return res;
}

double LinearCombination::density(double r) const
{
    const double res = wave_function(r);
    return res*res;
}

double LinearCombination::density(double r, double phi) const
{
    const Complex res = wave_function(r, phi);
    return res.real()*res.real() + res.imag()*res.imag();
}

} // namespace Mbs
