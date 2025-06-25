#include <iostream>
#include <fstream>
#include "fock/Bogoliubov.h"
#include "fock/MeanFieldSolver.h"
#include "fock/Numeric.h"

int main()
{
    std::vector<double> gns = Mbs::linspace(4.96552, 18.0, 22);
    Mbs::MeanFieldSolver mfs0(0);
    Mbs::MeanFieldSolver mfs1(1);

    std::ofstream file0("energies_gpe_bdg_0.txt");
    std::ofstream file1("energies_gpe_bdg_1.txt");

    file0 << "# gn mf mf+bdg dN" << std::endl;
    file1 << "# gn mf mf+bdg dN" << std::endl;
    for (double gn : gns) {
        std::cout << "### Running with gn " << gn << " ###" << std::endl;
        mfs0.solve(gn/6, 6, 10);
        double mf_ener = mfs0.energy();
        auto bdg = Mbs::bdg_correction(6, gn/6, 0);
        file0 << gn << " " << mf_ener << " " << mf_ener + bdg.dE << " " << bdg.dN << std::endl;

        mfs1.solve(gn/6, 6, 10);
        mf_ener = mfs1.energy();
        bdg = Mbs::bdg_correction(6, gn/6, 1);
        file1 << gn << " " << mf_ener << " " << mf_ener + bdg.dE << " " << bdg.dN << std::endl;
    }
    file0.close();
    file1.close();

    return 0;
}
