#include "fock/MeanFieldSolver.h"
#include "fock/DiskState.h"
#include "fock/Numeric.h"
#include "Eigen/Dense"
#include "Spectra/SymEigsShiftSolver.h"
#include "gsl/gsl_sf_bessel.h"

namespace Mbs {

MeanFieldSolver::MeanFieldSolver(int ang_mom, uint lattice_size, uint sp_states) :
    m_lattice_size(lattice_size),
    m_sp_states_size(sp_states),
    m_dr(1.0 / (lattice_size - 1)),
    m_sp_ang_mom(ang_mom),
    m_r(lattice_size),
    m_sp_states(sp_states),
    m_bessel(lattice_size*lattice_size),
    m_rho(lattice_size),
    m_veff(lattice_size),
    m_work(lattice_size),
    m_coeff(Eigen::VectorXd::Zero(sp_states))
{
    for (uint i = 0; i < lattice_size; ++i)
        m_r[i] = i * m_dr;

    u16 j = 0;
    for (uint i = 0; i < sp_states; ++i) {
        while (angular_momentum(j) != ang_mom)
            ++j;
        m_sp_states[i] = j;
        ++j;
    }

    m_coeff(0) = 1.0;
}

void MeanFieldSolver::solve(double g, uint N, uint eta)
{
    Eigen::MatrixXd H(m_sp_states_size, m_sp_states_size);

    for (uint i = 0; i < m_lattice_size*m_lattice_size; ++i)
        m_bessel[i] = gsl_sf_bessel_In(0, 2 * i * m_dr * m_dr * eta * eta);

    constexpr uint max_iter = 10'000;
    uint iter = 0;

    double mu_old = 0.0;
    double mu_diff = 1.0;
    while (mu_diff > 1e-10 && iter < max_iter) {
        compute_rho();
        compute_veff(g, N, eta);

        // Fill the 'Hamiltonian' matrix
        for (uint i = 0; i < m_sp_states_size; ++i) {
            for (uint j = i; j < m_sp_states_size; ++j) {
                for (uint k = 0; k < m_lattice_size; ++k) {
                    const double v = wave_function(m_sp_states[i], m_r[k]) * wave_function(m_sp_states[j], m_r[k]);
                    m_work[k] = v * m_r[k] * m_veff[k];
                }
                const double v = 2 * M_PI * simpson(m_work.data(), m_lattice_size, m_dr);
                H(i, j) = v;
                H(j, i) = H(i, j);
            }
            H(i, i) += ::Mbs::energy(m_sp_states[i]);
        }

        // Diagonalize the 'Hamiltonian' matrix
        
        Spectra::DenseSymShiftSolve<double> op(H);
        Spectra::SymEigsShiftSolver<Spectra::DenseSymShiftSolve<double>> solver(op, 3, 6, 0.0);
        solver.init();
        solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestAlge);
        if (solver.info() != Spectra::CompInfo::Successful)
            throw std::runtime_error("Failed to solve the eigenvalue problem");

        // See the progress
        constexpr double mixing = 0.8;
        m_coeff = mixing * solver.eigenvectors().col(0) + (1 - mixing) * m_coeff;
        m_coeff.normalize();
        double mu = solver.eigenvalues()[0];
        mu_diff = std::abs((mu - mu_old)/mu);
        mu_old = mu;
        ++iter;
    }

    if (iter == max_iter)
        throw std::runtime_error("Failed to converge the mean field solution");

    // COMPUTING ENERGY

    compute_rho();
    compute_veff(g, N, eta);

    // Compute energy
    for (uint j = 0; j < m_lattice_size; ++j)
        m_work[j] = m_veff[j] * m_rho[j] * m_r[j];
    double Epp = mu_old - M_PI * simpson(m_work.data(), m_lattice_size, m_dr);

    m_E = Epp * N;
    m_mu = mu_old;
}

LinearCombination MeanFieldSolver::get_state() const {
    LinearCombination lc(m_sp_states_size, m_sp_ang_mom);
    for (uint i = 0; i < m_sp_states_size; ++i)
        lc[i] = m_coeff[i];
    lc.normalize();
    return lc;
}

std::vector<double> MeanFieldSolver::veff(uint lattice_size, double g, uint N, uint eta) const
{
    const double dr = 1.0 / (lattice_size - 1);
    std::vector<double> veff(lattice_size);
    if (eta == 0) {
        for (uint i = 0; i < lattice_size; ++i)
            veff[i] = g * (N-1) * m_rho[i];
        return veff;
    }

    std::vector<double> work(lattice_size);
    for (uint i = 0; i < lattice_size; ++i) {
        for (uint j = 0; j < lattice_size; ++j) {
            double v = (j*dr) * std::exp(-(j*dr) * (j*dr) * eta * eta) * m_rho[j];
            v *= gsl_sf_bessel_In(0, 2*eta*eta*(i*dr)*(j*dr));
            work[j] = v;
        }

        veff[i] = simpson(work.data(), lattice_size, dr);
        veff[i] *= 2.0*g*eta*eta*(N - 1)*std::exp(-(i*dr)*(i*dr)*eta*eta);
    }
    return veff;
}

void MeanFieldSolver::compute_rho()
{
    for (uint i = 0; i < m_lattice_size; ++i) {
        double v = 0.0;
        for (uint j = 0; j < m_sp_states_size; ++j)
            v += m_coeff[j] * wave_function(m_sp_states[j], m_r[i]);
        m_rho[i] = v * v;
    }
}

void MeanFieldSolver::compute_veff(double g, uint N, uint eta)
{
    if (eta == 0) {
        for (uint i = 0; i < m_lattice_size; ++i)
            m_veff[i] = g * (N-1) * m_rho[i];
        return;
    }

    for (uint i = 0; i < m_lattice_size; ++i) {
        for (uint j = 0; j < m_lattice_size; ++j) {
            double v = m_r[j] * std::exp(-m_r[j] * m_r[j] * eta * eta) * m_rho[j];
            m_work[j] = v * m_bessel[i*j];
        }

        m_veff[i] = simpson(m_work.data(), m_lattice_size, m_dr);
        m_veff[i] *= 2.0*g*eta*eta*(N - 1)*std::exp(-m_r[i]*m_r[i]*eta*eta);
    }
}

}
