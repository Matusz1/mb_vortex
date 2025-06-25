#include <iostream>
#include <fstream>
#include "fock/DiskState.h"
#include "fock/MeanFieldSolver.h"
#include "fock/Integrator.h"
#include "fock/ImportanceTruncation.h"

// Retarded compiler producing warnings ...
static Mbs::LinearCombination calc_lowest(const Eigen::Ref<const Eigen::MatrixXd>& rho)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(rho);
    const auto& last_evec = es.eigenvectors().col(rho.rows()-1);
    const double* lowest = last_evec.data();
    std::vector<u16> sp_states(rho.rows());
    std::iota(sp_states.begin(), sp_states.end(), 0);
    Mbs::LinearCombination ret(lowest, sp_states.data(), rho.rows());
    return ret;
}

int main()
{
    std::vector<int> gns = {0, 6, 12, 18};

    Mbs::TruncationParameters params;
    params.ener_cutoff = 800.0;
    params.C_min = 1e-4;

    Mbs::integrator_setup(10);
    for (double gn : gns) {
        std::cout << "### Running with gn " << gn << " ###" << std::endl;

        // The L = 0
        Mbs::Workspace ws = Mbs::importance_truncation_scheme(gn/6, 6, 0, params);
        Mbs::FockLinear evec = ws.eigenvector_combination(0);
        Eigen::MatrixXd rho = evec.density_matrix();
        Mbs::LinearCombination lowest = calc_lowest(rho);
        Mbs::MeanFieldSolver mfs(0);
        mfs.solve(gn/6, 6, 10);
        Mbs::LinearCombination mfsol = mfs.get_state();
        std::fstream file("density_" + std::to_string(gn) + "_0.txt", std::ios::out);
        file << "# r rho(r) rho0(r) mf" << std::endl;
        for (double r = 0; r < 1.005; r += 0.01) {
            double diag = Mbs::density_matrix_compute(r, rho);
            double phi0 = lowest.density(r);
            double mf = mfsol.density(r);
            file << r << " " << diag << " " << phi0 << " " << mf << std::endl;
        }
        file.close();

        // The L = N
        ws = Mbs::importance_truncation_scheme(gn/6, 6, 6, params);
        evec = ws.eigenvector_combination(0);
        rho = evec.density_matrix();
        lowest = calc_lowest(rho);
        mfs = Mbs::MeanFieldSolver(1);
        mfs.solve(gn/6, 6, 10);
        mfsol = mfs.get_state();
        file.open("density_" + std::to_string(gn) + "_1.txt", std::ios::out);
        file << "# r rho(r) rho0(r) mf" << std::endl;
        for (double r = 0; r < 1.005; r += 0.01) {
            double diag = Mbs::density_matrix_compute(r, rho);
            double phi0 = lowest.density(r);
            double mf = mfsol.density(r);
            file << r << " " << diag << " " << phi0 << " " << mf << std::endl;
        }
        file.close();

        // The L = N/2
        ws = Mbs::importance_truncation_scheme(gn/6, 6, 3, params);
        evec = ws.eigenvector_combination(0);
        rho = evec.density_matrix();
        lowest = calc_lowest(rho);
        file.open("density_" + std::to_string(gn) + "_Nh.txt", std::ios::out);
        file << "# r rho(r) rho0(r)" << std::endl;
        for (double r = 0; r < 1.005; r += 0.01) {
            double diag = Mbs::density_matrix_compute(r, rho);
            double phi0 = lowest.density(r);
            file << r << " " << diag << " " << phi0 << std::endl;
        }
    }

    return 0;
}
