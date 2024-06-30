#include "src/examples.h"

int main(int argc, char **argv) {
        /*
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(6, 130.0, 6);
        FockBasis b2 = fock_basis_alloc(6, 130.0);
        uint max = 0;
        for (uint i = 0; i < basis.size; ++i)
                for (uint j = 0; j < 6; ++j)
                        max = MAX(max, basis.states[i].states[j]);
        ++max;
        printf("# size = %u\n", basis.size);
        printf("# k    = %u\n", max);
        printf("# b2   = %u\n", b2.size);
        */

        //example_random_sampling(100.0, 6);
        //example_gpe_energy(10);
        //example_rho_evec(10, 6.0, 3000);
        //example_bdg(10, 18.0, 100);
        //example_rho_diag(10, 0.0, 3.0, 20, 1500);
        //example_rho_diag(10, 4.0, 18.0, 15, 3000);
        //example_bdg_write(0, 12.0, false);
        example_bdg_uv(10, 12.0);
        
        /*
        double gamma_vals[6] = {
                0.0,
                6.0,
                12.0,
                18.0
        };
        for (uint i = 0; i < sizeof(gamma_vals)/sizeof(gamma_vals[0]); ++i) {
                double gamma = gamma_vals[i];
                example_rho_gpe_compariton(6, gamma, 3000);
        }
        */

	return 0;
}

