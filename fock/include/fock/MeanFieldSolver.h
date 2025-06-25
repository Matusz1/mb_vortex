#ifndef _MEAN_FIELD_SOLVER_H_
#define _MEAN_FIELD_SOLVER_H_

#include "Eigen/Dense"
#include "DiskState.h"

namespace Mbs {

class MeanFieldSolver {
private:
    uint m_lattice_size;
    uint m_sp_states_size;
    double m_dr;
    int m_sp_ang_mom;

    std::vector<double> m_r;
    std::vector<u16> m_sp_states;
    std::vector<double> m_bessel;

    std::vector<double> m_rho;
    std::vector<double> m_veff;
    std::vector<double> m_work;
    Eigen::VectorXd m_coeff;

    double m_mu;
    double m_E;
public:
    explicit MeanFieldSolver(int sp_ang_mom, uint lattice_size = 101, uint sp_states = 20);

    void solve(double g, uint N, uint eta);

    [[nodiscard]] double mu() const { return m_mu; }
    [[nodiscard]] double energy() const { return m_E; }

    LinearCombination get_state() const;
    std::vector<double> veff(uint lattice_size, double g, uint N, uint eta) const;

private:
    void compute_rho();
    void compute_veff(double g, uint N, uint eta);
};

}

#endif // _MEAN_FIELD_SOLVER_H_
