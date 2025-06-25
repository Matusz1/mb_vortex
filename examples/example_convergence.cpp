#include <iostream>
#include <iomanip>
#include <fstream>
#include "fock/ImportanceTruncation.h"
#include "fock/Integrator.h"
#include "fock/Numeric.h"
#include <vector>

int main()
{
    std::ofstream file("convergence_cutoffs.txt");
    std::vector<double> cutoffs = Mbs::linspace(100.0, 1200.0, 12);

    Mbs::TruncationParameters params;
    params.kappa_min = 1e-5;
    params.C_min = 1e-4;

    Mbs::integrator_setup(10);
    for (double cut : cutoffs) {
        std::cout << "### Running with cutoff " << cut << " ###" << std::endl;
        params.ener_cutoff = cut;
        auto ws = Mbs::importance_truncation_scheme(18.0/6, 6, 0, params);
        const double energy = ws.eigenvalues()[0];
        file << cut << " " << std::setprecision(10) << energy << std::endl;
        std::cout << "### Run done ###\n" << std::endl;
    }
    Mbs::integrator_cleanup();
    file.close();

    return 0;
}
