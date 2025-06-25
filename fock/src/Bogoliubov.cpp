#include "fock/Bogoliubov.h"
#include "fock/MeanFieldSolver.h"
#include "fock/DiskState.h"
#include "fock/Numeric.h"
#include "gsl/gsl_sf_bessel.h"
#include <map>

namespace Mbs {

SolutionBdG::SolutionBdG(uint n, int upper_mom, int lower_mom)
    : m_num(n)
    , m_u_basis(n)
    , m_v_basis(n)
    , m_uv_matrix(2 * n, 2 * n)
    , m_energy(2 * n)
    , m_family(2 * n)
{
    u16 i = 0;
    for (uint j = 0; j < n; ++j) {
        while (angular_momentum(i) != upper_mom)
            ++i;
        m_u_basis[j] = i;
        ++i;
    }

    i = 0;
    for (uint j = 0; j < n; ++j) {
        while (angular_momentum(i) != lower_mom)
            ++i;
        m_v_basis[j] = i;
        ++i;
    }
}

double SolutionBdG::compute_u(uint i, double r) const
{
    double val = 0.0;
    for (uint j = 0; j < m_num; ++j)
        val += m_uv_matrix(j, i).real() * wave_function(m_u_basis[j], r);
    return val;
}

double SolutionBdG::compute_v(uint i, double r) const
{
    double val = 0.0;
    for (uint j = 0; j < m_num; ++j)
        val += m_uv_matrix(j + m_num, i).real() * wave_function(m_v_basis[j], r);
    return val;
}

SolutionBdG solve_BdG(uint N, double g, int mom, int upper_mom, double omega)
{
    assert((mom == 0 || mom == 1) && "Only mom = 0 or mom = 1 is supported in solve_BdG");
    constexpr uint lattice_size = 401;
    constexpr double dr = 1.0 / (lattice_size - 1);
    constexpr double eta = 10;
    constexpr uint n_sp = 50;

    MeanFieldSolver mfs(mom, lattice_size, n_sp);
    mfs.solve(g, N, eta);
    const LinearCombination state = mfs.get_state();

    std::vector<double> work_in(lattice_size);
    std::vector<double> work_out(lattice_size);
    std::vector<double> psi0(lattice_size);
    std::vector<double> veff = mfs.veff(lattice_size, g, N, eta);

    for (uint i = 0; i < lattice_size; ++i)
        psi0[i] = state.wave_function(i * dr);

    const int lower_mom = upper_mom - 2*mom;
    SolutionBdG sol(n_sp, upper_mom, lower_mom);

    // Create and fille the L_GP matrix
    Eigen::MatrixXd lgp = Eigen::MatrixXd::Zero(2*n_sp, 2*n_sp);

    std::vector<double> bessel_up(lattice_size*lattice_size);
    std::vector<double> bessel_down(lattice_size*lattice_size);
    for (uint i = 0; i < lattice_size*lattice_size; ++i) {
        bessel_up[i] = gsl_sf_bessel_In(std::abs(upper_mom-mom), 2*eta*eta*i*dr*dr);
        bessel_down[i] = gsl_sf_bessel_In(std::abs(lower_mom+mom), 2*eta*eta*i*dr*dr);
    }

    std::map<u16, std::vector<double>> wfs_up;
    std::map<u16, std::vector<double>> wfs_down;
    for (uint i = 0; i < n_sp; ++i) {
        const u16 u = sol.m_u_basis[i];
        const u16 v = sol.m_v_basis[i];
        std::vector<double> wf_up(lattice_size);
        std::vector<double> wf_down(lattice_size);
        for (uint j = 0; j < lattice_size; ++j) {
            wf_up[j] = wave_function(u, j*dr);
            wf_down[j] = wave_function(v, j*dr);
        }
        wfs_up[u] = std::move(wf_up);
        wfs_down[v] = std::move(wf_down);
    }


    for (uint i = 0; i < n_sp; ++i) {
        const std::vector<double>& lu = wfs_up[sol.m_u_basis[i]];
        const std::vector<double>& lv = wfs_down[sol.m_v_basis[i]];
        for (uint j = 0; j < n_sp; ++j) {
            const std::vector<double>& ru = wfs_up[sol.m_u_basis[j]];
            const std::vector<double>& rv = wfs_down[sol.m_v_basis[j]];

            // u*u integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                for (uint jj = 0; jj < lattice_size; ++jj) {
                    work_in[jj] = jj*dr * psi0[jj] * ru[jj];
                    work_in[jj] *= bessel_up[ii*jj] * std::exp(-eta*eta*(jj*dr)*(jj*dr));
                }
                work_out[ii] = simpson(work_in.data(), lattice_size, dr) * std::exp(-eta*eta*(ii*dr)*(ii*dr));
                work_out[ii] *= ii*dr * psi0[ii] * lu[ii];
            }
            double val = 4 * M_PI * eta * eta * simpson(work_out.data(), lattice_size, dr);
            lgp(i, j) = g * (N-1) * val;

            // v*v integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                for (uint jj = 0; jj < lattice_size; ++jj) {
                    work_in[jj] = jj*dr * psi0[jj] * rv[jj];
                    work_in[jj] *= bessel_down[ii*jj] * std::exp(-eta*eta*(jj*dr)*(jj*dr));
                }
                work_out[ii] = simpson(work_in.data(), lattice_size, dr) * std::exp(-eta*eta*(ii*dr)*(ii*dr));
                work_out[ii] *= ii*dr * psi0[ii] * lv[ii];
            }
            val = 4 * M_PI * eta * eta * simpson(work_out.data(), lattice_size, dr);
            lgp(n_sp+i, n_sp+j) = - g * (N-1) * val;

            // u*v integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                for (uint jj = 0; jj < lattice_size; ++jj) {
                    work_in[jj] = jj*dr * psi0[jj] * rv[jj];
                    work_in[jj] *= bessel_up[ii*jj] * std::exp(-eta*eta*(jj*dr)*(jj*dr));
                }
                work_out[ii] = simpson(work_in.data(), lattice_size, dr) * std::exp(-eta*eta*(ii*dr)*(ii*dr));
                work_out[ii] *= ii*dr * psi0[ii] * lu[ii];
            }
            val = 4 * M_PI * eta * eta * simpson(work_out.data(), lattice_size, dr);
            lgp(i, n_sp+j) = g * (N-1) * val;

            // v*u integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                for (uint jj = 0; jj < lattice_size; ++jj) {
                    work_in[jj] = jj*dr * psi0[jj] * ru[jj];
                    work_in[jj] *= bessel_down[ii*jj] * std::exp(-eta*eta*(jj*dr)*(jj*dr));
                }
                work_out[ii] = simpson(work_in.data(), lattice_size, dr) * std::exp(-eta*eta*(ii*dr)*(ii*dr));
                work_out[ii] *= ii*dr * psi0[ii] * lv[ii];
            }
            val = 4 * M_PI * eta * eta * simpson(work_out.data(), lattice_size, dr);
            lgp(n_sp+i, j) = - g * (N-1) * val;

            // u*u Veff integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                work_out[ii] = ii*dr * veff[ii];
                work_out[ii] *= ru[ii] * lu[ii];
            }
            val = 2 * M_PI * simpson(work_out.data(), lattice_size, dr);
            lgp(i,j) += val;

            // v*v Veff integral
            for (uint ii = 0; ii < lattice_size; ++ii) {
                work_out[ii] = ii*dr * veff[ii];
                work_out[ii] *= rv[ii] * lv[ii];
            }
            val = 2 * M_PI * simpson(work_out.data(), lattice_size, dr);
            lgp(n_sp+i,n_sp+j) -= val;
        }

        // Kinetic energy and chemical potential
        lgp(i,i) += energy(sol.m_u_basis[i]) - omega*upper_mom - (mfs.mu() - omega*mom);
        lgp(n_sp+i,n_sp+i) += -energy(sol.m_v_basis[i]) - omega*lower_mom + (mfs.mu() - omega*mom);
    }

    // The projector
    if (upper_mom == mom) {
        Eigen::MatrixXd projQ = Eigen::MatrixXd::Identity(2*n_sp, 2*n_sp);
        for (uint i = 0; i < n_sp; ++i)
            for (uint j = 0; j < n_sp; ++j)
                projQ(i,j) -= state[i] * state[j];
        projQ.bottomRightCorner(n_sp, n_sp) = projQ.topLeftCorner(n_sp, n_sp);
        lgp = lgp * projQ;
    }

    // Solve the eigenvalue problem
    Eigen::EigenSolver<Eigen::MatrixXd> ces(lgp);
    sol.m_uv_matrix = ces.eigenvectors();
    sol.m_energy = ces.eigenvalues();

    // Determine the family, normalize to +-1
    for (uint i = 0; i < 2*n_sp; ++i) {
        double norm = sol.u(i).squaredNorm() - sol.v(i).squaredNorm();
        // TODO: Check this
        if (std::abs(sol.energy(i)) < 1e-6 || std::abs(sol.energy(i).imag()) > 1e-6) {
            sol.m_family[i] = SolutionBdG::Family::Zero;
            continue;
        } else if (norm > 0) {
            sol.m_family[i] = SolutionBdG::Family::Plus;
        } else {
            sol.m_family[i] = SolutionBdG::Family::Minus;
        }
        if (norm != 0.0) {
            norm = 1.0 / std::sqrt(norm);
            sol.m_uv_matrix.block(0, i, 2*n_sp, 1) *= norm;
        }
    }

    return sol;
}

BdGcomputed bdg_correction(uint N, double g, int mom)
{
    BdGcomputed res;
    res.dN = 0.0;
    res.dE = 0.0;

    for (int umom = -5; umom < 5; ++umom) {
        SolutionBdG sol = solve_BdG(N, g, mom, umom, mom*5.5);
        double dN = 0.0;
        double dE = 0.0;
        for (uint i = 0; i < sol.size(); ++i) {
            if (sol.family(i) != SolutionBdG::Family::Plus)
                continue;
            const double norm = sol.v(i).squaredNorm();
            dN -= norm;
            dE -= norm * sol.energy(i).real();
        }
        res.dN += dN;
        res.dE += dE;
    }
    return res;
}

}
