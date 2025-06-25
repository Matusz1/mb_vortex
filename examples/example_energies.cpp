#include <iostream>
#include <fstream>
#include "fock/ImportanceTruncation.h"
#include "fock/Integrator.h"
#include "fock/Numeric.h"

// Retarded compiler producing warnings ...
static double calc_n0(const Eigen::Ref<const Eigen::MatrixXd>& rho)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(rho);
    return es.eigenvalues().maxCoeff();
}

int main()
{
    std::vector<double> gns = Mbs::linspace(0.0, 18.0, 30);

    std::ofstream file0("energies_mb_0.txt");
    std::ofstream file1("energies_mb_1.txt");

    file0 << "# gn mb lambda0" << std::endl;
    file1 << "# gn mb lambda0" << std::endl;

    Mbs::TruncationParameters params;
    params.ener_cutoff = 800.0;
    params.C_min = 1e-4;

    Mbs::integrator_setup(10);
    for (double gn : gns) {
        std::cout << "### Running with gn " << gn << " ###" << std::endl;
        Mbs::Workspace ws = Mbs::importance_truncation_scheme(gn/6, 6, 0, params);
        double ener = ws.eigenvalues()[0];
        Mbs::FockLinear evec = ws.eigenvector_combination(0);
        Eigen::MatrixXd rho = evec.density_matrix();
        double lambda0 = calc_n0(rho) / 6;
        file0 << gn << " " << ener << " " << lambda0 << std::endl;

        ws = Mbs::importance_truncation_scheme(gn/6, 6, 6, params);
        ener = ws.eigenvalues()[0];
        evec = ws.eigenvector_combination(0);
        rho = evec.density_matrix();
        lambda0 = calc_n0(rho) / 6;
        file1 << gn << " " << ener << " " << lambda0 << std::endl;
    }
    file0.close();
    file1.close();
}
