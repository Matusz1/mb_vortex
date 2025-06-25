#ifndef _BOGOLIUBOV_H_
#define _BOGOLIUBOV_H_
#include <Eigen/Dense>

namespace Mbs {

class SolutionBdG {
public:
    enum class Family {
        Minus = -1,
        Zero = 0,
        Plus = 1
    };

private:
    uint m_num;
    std::vector<uint> m_u_basis;
    std::vector<uint> m_v_basis;

    Eigen::MatrixXcd m_uv_matrix;
    Eigen::VectorXcd m_energy;
    std::vector<Family> m_family;
public:
    SolutionBdG(uint n, int upper_mom, int lower_mom);

    std::size_t size() const { return 2*m_num; }

    Family family(uint i) const { return m_family[i]; }
    auto u(uint i) const { return m_uv_matrix.block(0, i, m_num, 1); }
    auto v(uint i) const { return m_uv_matrix.block(m_num, i, m_num, 1); }
    double compute_u(uint i, double r) const;
    double compute_v(uint i, double r) const;
    auto energy(uint i) const { return m_energy(i); }

    friend SolutionBdG solve_BdG(uint N, double g, int mom, int upper_mom, double omega);
};

SolutionBdG solve_BdG(uint N, double g, int mom, int upper_mom, double omega = 0.0);

struct BdGcomputed {
    double dN;
    double dE;
};

BdGcomputed bdg_correction(uint N, double g, int mom);

}

#endif // _BOGOLIUBOV_H_
