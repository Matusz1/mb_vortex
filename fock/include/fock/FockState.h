#ifndef _FOCK_STATE_H_
#define _FOCK_STATE_H_

#include <array>
#include <map>
#include <cassert>
#include <Eigen/Dense>

#include "DiskState.h"


namespace Mbs {

constexpr uint CONFING_MAX_FOCK_SIZE = 10;

template<typename Scalar>
class FockLinear;

class Fock {
public:
    std::array<u16, CONFING_MAX_FOCK_SIZE+1> m_state;
    uint m_size;

public:
    explicit Fock(uint size) : m_size(size) { m_state.fill(0); }
    Fock(uint size, const u16* state) : m_size(size) { std::copy(state, state + size, m_state.begin()); }

    [[nodiscard]] double energy() const;
    [[nodiscard]] int angular_momentum() const;
    [[nodiscard]] uint count(u16 i) const;

    [[nodiscard]] uint size() const { return m_size; }
    [[nodiscard]] u16 operator[](uint i) const { return m_state[i]; }
    [[nodiscard]] u16& operator[](uint i) { return m_state[i]; }

    double permanent_norm() const;

    double annihilate(u16 i);
    double create(u16 i);

    [[nodiscard]] bool operator==(const Fock& other) const;
    [[nodiscard]] bool operator!=(const Fock& other) const { return !(*this == other); };
    [[nodiscard]] bool operator<(const Fock& other) const;
    [[nodiscard]] bool operator>(const Fock& other) const { return other < *this; };
    [[nodiscard]] bool operator<=(const Fock& other) const { return !(*this > other); };
    [[nodiscard]] bool operator>=(const Fock& other) const { return !(*this < other); };

    class Hash {
    public:
        std::size_t operator()(const Fock& fock) const
        {
            std::size_t hash = 0;
            for (uint i = 0; i < fock.size(); ++i)
                hash ^= std::hash<u16>{}(fock[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            return hash;
        }
    };

    template<typename Scalar>
    friend class FockLinear;
};

template<typename Scalar>
FockLinear<Complex> apply_psi(double r, double phi, const FockLinear<Scalar>& psi);

template<typename Scalar>
class FockLinear {
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Complex>::value);

private:
    std::map<Fock, Scalar> m_states;

    FockLinear() = default;
public:
    using DensityMatrix = Eigen::MatrixX<Scalar>;

    std::size_t size() const { return m_states.size(); }
    uint npart() const { return m_states.begin()->first.size(); }

    DensityMatrix density_matrix() const;

    [[nodiscard]] FockLinear annihilate(u16 i) const;
    [[nodiscard]] FockLinear create(u16 i) const;

    void reduce(double cutoff);
    void normalize();
    Scalar product(const FockLinear& other) const;
    FockLinear& operator+=(const FockLinear& other);
    FockLinear& operator*=(Scalar factor);

    FockLinear<Complex> complexify() const;

    auto begin() const { return m_states.begin(); }
    auto end() const { return m_states.end(); }

    
    template<typename Ty>
    friend FockLinear<Complex> apply_psi(double r, double phi, const FockLinear<Ty>& psi);

    template<typename Ty>
    friend class FockLinear;

    friend class Workspace;
    friend class Fock;
};

template<typename Scalar>
double density_matrix_compute(double r, double phi, const Eigen::MatrixX<Scalar>& rho)
{
    const uint size = rho.rows();
    Eigen::VectorXcd vec(size);
    for (u16 i = 0; i < size; ++i)
        vec(i) = wave_function(i, r, phi);
    Complex tmp = vec.adjoint() * rho * vec;
    return tmp.real();
}

template<typename Scalar>
inline double density_matrix_compute(double r, const Eigen::MatrixX<Scalar>& rho)
{
    return density_matrix_compute(r, 0.0, rho);
}

[[nodiscard]] double potential_element(const Fock& left, const Fock& right);
[[nodiscard]] double potential_element_unchecked(const Fock& left, const Fock& right);
[[nodiscard]] Complex wave_function(const u16* state, const double* r, const double* phi, uint size);

inline double Fock::energy() const
{
    double ener = 0.0;
    for (uint i = 0; i < size(); ++i)
        ener += ::Mbs::energy(m_state[i]);
    return ener;
}

inline int Fock::angular_momentum() const
{
    int ang = 0;
    for (uint i = 0; i < size(); ++i)
        ang += ::Mbs::angular_momentum(m_state[i]);
    return ang;
}

inline uint Fock::count(u16 i) const
{
    uint num = 0;
    for (uint j = 0; j < size(); ++j)
        if (m_state[j] == i)
            ++num;
    return num;
}

inline double density(const u16* state, const double* r, const double* phi, uint size)
{
    const double v = std::abs(wave_function(state, r, phi, size));
    return v*v;
}

inline double phase(const u16* state, const double* r, const double* phi, uint size)
{
    return std::arg(wave_function(state, r, phi, size));
}

inline void find_creation_operators(const u16* large, const u16* small, uint ssize, u16* out)
{
    uint i = 0;
    while (i < ssize && large[i] == small[i])
        ++i;
    out[0] = large[i];
    while (i < ssize && large[i+1] == small[i])
        ++i;
    out[1] = large[i+1];
}

inline uint fock_mismatch(const u16* left, const u16* right, uint size)
{
    uint num_pairs = 0;
    const u16* lend = left + size;
    const u16* rend = right + size;
    while (left != lend && right != rend) {
        if (*left == *right) {
            ++num_pairs;
            ++left;
            ++right;
        } else if (*left < *right) {
            ++left;
        } else {
            ++right;
        }
    }
    return size - num_pairs;
}

inline bool fock_reachable(const u16* large, const u16* small, uint size)
{
    uint i = 0;
    uint j = 0;
    while (i < size && j < size+2) {
        if (small[i] == large[j])
            ++i;
        ++j;
    }
    return i == size;
}

template<typename Scalar>
FockLinear<Scalar> FockLinear<Scalar>::annihilate(u16 i) const
{
    FockLinear<Scalar> ret;
    for (auto it = m_states.begin(); it != m_states.end(); ++it) {
        Fock fock = it->first;
        const double factor = fock.annihilate(i);
        if (factor != 0.0)
            ret.m_states[fock] += it->second * factor;
    }
    return ret;
}

template<typename Scalar>
FockLinear<Scalar> FockLinear<Scalar>::create(u16 i) const
{
    FockLinear<Scalar> ret;
    for (auto it = m_states.begin(); it != m_states.end(); ++it) {
        Fock fock = it->first;
        const double factor = fock.create(i);
        ret.m_states[fock] += it->second * factor;
    }
    return ret;
}

template<typename Scalar>
void FockLinear<Scalar>::reduce(double cutoff)
{
    for (auto it = m_states.begin(); it != m_states.end();) {
        if (std::abs(it->second) < cutoff)
            it = m_states.erase(it);
        else
            ++it;
    }
}

template<typename Scalar>
void FockLinear<Scalar>::normalize()
{
    double norm = 0.0;
    for (const auto& [_, factor] : m_states)
        norm += std::abs(factor * factor);
    norm = 1.0 / std::sqrt(norm);
    for (auto& [_, factor] : m_states)
        factor *= norm;
}

template<typename Scalar>
Scalar FockLinear<Scalar>::product(const FockLinear& other) const
{
    Scalar ret = 0.0;
    for (const auto& [right, right_factor] : m_states) {
        const auto it = other.m_states.find(right);
        if (it != other.m_states.end()) {
            if constexpr (std::is_same<Scalar, double>::value)
                ret += right_factor * it->second;
            else
                ret += std::conj(right_factor) * it->second;
        }
    }
    return ret;
}

template<typename Scalar>
FockLinear<Scalar>& FockLinear<Scalar>::operator+=(const FockLinear& other)
{
    for (const auto& [right, right_factor] : other.m_states)
        m_states[right] += right_factor;
    return *this;
}

template<typename Scalar>
FockLinear<Scalar>& FockLinear<Scalar>::operator*=(Scalar factor)
{
    for (auto& [right, right_factor] : m_states)
        right_factor *= factor;
    return *this;
}

template<typename Scalar>
typename FockLinear<Scalar>::DensityMatrix FockLinear<Scalar>::density_matrix() const
{
    u16 size = 0;
    for (const auto& [fock, _] : m_states)
        if (fock[fock.size() - 1] > size)
            size = fock[fock.size() - 1];

    ++size;
    DensityMatrix ret = DensityMatrix::Zero(size, size);

    for (u16 i = 0; i < size; ++i) {
        FockLinear fock = this->annihilate(i);
        for (u16 j = i; j < size; ++j) {
            const Scalar value = fock.create(j).product(*this);
            ret(i, j) = value;
            ret(j, i) = value;
        }
    }
    return ret;
}

template<typename Scalar>
FockLinear<Complex> apply_psi(double r, double phi, const FockLinear<Scalar>& psi)
{
    u16 size = 0;
    for (const auto& [fock, _] : psi.m_states)
        if (fock[fock.size() - 1] > size)
            size = fock[fock.size() - 1];
    ++size;

    FockLinear<Complex> ret;
    for (u16 i = 0; i < size; ++i) {
        FockLinear<Scalar> fock = psi.annihilate(i);
        fock *= wave_function(i, r, phi);
        ret += fock;
    }
    return ret;
}

template<typename Scalar>
FockLinear<Complex> FockLinear<Scalar>::complexify() const
{
    FockLinear<Complex> ret;
    for (const auto& [fock, factor] : m_states)
        ret.m_states[fock] = factor;
    return ret;
}

} // namespace Mbs

#endif // _FOCK_STATE_H_
