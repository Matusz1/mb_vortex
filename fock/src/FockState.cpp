#include <cassert>
#include <vector>
#include "fock/FockState.h"
#include "fock/DiskState.h"
#include "fock/Numeric.h"
#include "fock/AnnichilationAdvancer.h"
#include "fock/Integrator.h"

namespace Mbs {

bool Fock::operator==(const Fock& other) const
{
    const u16* left = m_state.data();
    const u16* right = other.m_state.data();
    return std::equal(left, left + m_size, right);
}

bool Fock::operator<(const Fock& other) const
{
    const u16* left = m_state.data();
    const u16* right = other.m_state.data();
    return std::lexicographical_compare(std::reverse_iterator(left + m_size),
        std::reverse_iterator(left), std::reverse_iterator(right + m_size),
        std::reverse_iterator(right));
}

Complex wave_function(const u16* state, const double* r, const double* phi, uint size)
{
    std::vector<Complex> arr(size*size);
    for (uint i = 0; i < size; ++i)
        for (uint j = 0; j < size; ++j)
            arr[i*size + j] = wave_function(state[i], r[j], phi[j]);

    uint prev = 0;
    uint prod = 1;
    for (uint i = 1; i < size; ++i) {
        prod *= i+1;
        if (state[i] == state[i - 1])
            prod *= i - prev;
        else
            prev = i;
    }
    return permanent(arr.data(), size) / std::sqrt(prod);
}

double Fock::annihilate(u16 i)
{
    // TODO: Check if simple find is faster (most likely for our small sizes)
    const auto end = m_state.begin() + m_size;
    const auto lb = std::lower_bound(m_state.begin(), end, i);
    const auto ub = std::upper_bound(lb, end, i);
    if (lb == ub) {
        m_size = 0;
        return 0.0;
    }
    std::move(ub, end, ub - 1);
    --m_size;
    return std::sqrt(ub - lb);
}

double Fock::create(u16 i)
{
    // TODO: As in annihilate, check if simple find is faster
    const auto end = m_state.begin() + m_size;
    const auto lb = std::lower_bound(m_state.begin(), end, i);
    const auto ub = std::upper_bound(lb, end, i);
    std::move_backward(ub, end, end + 1);
    *ub = i;
    ++m_size;
    return std::sqrt(ub - lb + 1);
}

double potential_element_unchecked(const Fock& left, const Fock& right)
{
    uint size = left.size();
    Fock r = right;
    r.m_size = size - 2;

    // <ij|V|kl>
    std::array<u16,4> ijkl = {0, 0, 0, 0};

    double ret = 0.0;
    for (AnnichilationAdvancer<false> a1(r.m_state.data(), size); a1.valid(); a1.advance()) {
        const uint v1 = a1.count();
        ijkl[3] = a1.index();
        for (AnnichilationAdvancer<true> a2(r.m_state.data(), size-1); a2.valid(); a2.advance()) {
            // Check if contains
            if (!fock_reachable(left.m_state.data(), r.m_state.data(), size-2))
                continue;

            const uint v2 = a2.count();
            ijkl[2] = a2.index();
            find_creation_operators(left.m_state.data(), r.m_state.data(), size-2, ijkl.data());
            if (ijkl[0] != ijkl[1]) {
                const uint v3 = r.count(ijkl[0]) + 1;
                const uint v4 = r.count(ijkl[1]) + 1;
                double current = potential_integrate_cached(ijkl);
                std::swap(ijkl[0], ijkl[1]);
                current += potential_integrate_cached(ijkl);
                ret += std::sqrt(v1*v2*v3*v4) * current;
            } else {
                const uint v3 = r.count(ijkl[0]) + 1;
                double current = potential_integrate_cached(ijkl);
                ret += std::sqrt(v1*v2*v3*(v3+1)) * current;
            }
        }
    }

    return 0.5 * ret;
}

double potential_element(const Fock& left, const Fock& right)
{
    if (fock_mismatch(left.m_state.data(), right.m_state.data(), left.size()) > 2)
        return 0.0;
    return potential_element_unchecked(left, right);
}

double Fock::permanent_norm() const
{
    uint prev = 0;
    uint prod = 1;
    for (uint i = 1; i < m_size; ++i) {
        prod *= i+1;
        if (m_state[i] == m_state[i - 1])
            prod *= i - prev;
        else
            prev = i;
    }
    return 1.0 / std::sqrt(prod);
}

} // namespace Mbs
