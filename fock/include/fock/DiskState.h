#ifndef _DISK_STATE_H_
#define _DISK_STATE_H_

#include "Core.h"
#include <vector>

namespace Mbs {

/*
 * Single particle bessel states, they are indexed by
 * u16's, and for more detail about a state you can call
 * appropriate function (state_info(indx))
 */

[[nodiscard]] Complex wave_function(u16 n, double r, double phi);
[[nodiscard]] double wave_function(u16 n, double r);

[[nodiscard]] double energy(u16 n);
[[nodiscard]] int angular_momentum(u16 n);

[[nodiscard]] u16 disk_state_index(int ang_mom, uint zero_num);

class LinearCombination {
private:
    std::vector<double> m_coeff;
    std::vector<u16> m_basis;

public:
    LinearCombination(uint n, int mom);
    LinearCombination(const double* coeff, const u16* basis, uint size);

    uint size() const { return m_coeff.size(); }
    double& operator[](uint i) { return m_coeff[i]; }
    double operator[](uint i) const { return m_coeff[i]; }

    void normalize();

    double wave_function(double r) const;
    Complex wave_function(double r, double phi) const;

    double density(double r) const;
    double density(double r, double phi) const;
};

} // namespace Mbs

#endif // _DISK_STATE_H_
