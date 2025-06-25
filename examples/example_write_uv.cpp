#include <iostream>
#include <fstream>
#include "fock/Bogoliubov.h"
#include <iomanip>

int main()
{
    Mbs::SolutionBdG bdg = Mbs::solve_BdG(6, 18.0/6.0, 1, 0, 5.5);
    uint idx = 0;
    double dN = 0.0;
    for (uint i = 0; i < bdg.size(); ++i) {
        if (bdg.family(i) != Mbs::SolutionBdG::Family::Plus)
            continue;
        auto v = bdg.v(i);
        const double curr_dN = v.squaredNorm();
        std::cout << "dN = " << std::scientific << std::setprecision(3) << curr_dN << ", e = " << bdg.energy(i).real() << std::endl;
        if (curr_dN > dN) {
            dN = curr_dN;
            idx = i;
        }
    }
    /*std::cout << "dN = " << dN << std::endl;*/

    std::ofstream file("uv_0.txt");
    file << "# r u v" << std::endl;
    for (double r = -1.0; r < 1.001; r += 0.01) {
        double u_val = bdg.compute_u(idx, r);
        double v_val = bdg.compute_v(idx, r);
        file << r << " " << u_val << " " << v_val << std::endl;
    }

    return 0;
}
