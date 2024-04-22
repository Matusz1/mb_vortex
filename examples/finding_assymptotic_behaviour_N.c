#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>

#include <rgc_memory.h>
#include <rgc_io.h>

int main(int argc, char **argv) {

        // Energy cutofs for states with 0 angular momentum
        const double cutoff_energies[21] = {
                000.0, // Unused
                000.0, // Unused
                900.0, // 2  Particles
                300.0,
                190.0,
                150.0, // 5  Particles
                130.0,
                125.0,
                124.0,
                123.0,
                125.0, // 10 Particles
                125.0,
                126.0,
                130.0,
                133.0,
                136.0, // 15 Particles
                138.0,
                140.0,
                143.0,
                146.0,
                149.0, // 20 Particles
        };

        // Energy cutofs for states with 0 angular momentum
        /*
        const double cutoff_energies[21] = {
                000.0, // Unused
                000.0, // Unused
                500.0, // 2  Particles
                200.0,
                130.0,
                110.0, // 5  Particles
                100.0,
                100.0,
                100.0,
                100.0,
                100.0, // 10 Particles
                105.0,
                110.0,
                110.0,
                113.0,
                116.0, // 15 Particles
                118.0,
                120.0,
                123.0,
                126.0,
                129.0, // 20 Particles
        };
        */

        // Energy cutofs for states with N angular momentum
        /*
        const double cutoff_energies[21] = {
                000.0, // Unused
                000.0, // Unused
                900.0, // 2  Particles
                300.0,
                190.0,
                150.0, // 5  Particles
                140.0,
                135.0,
                134.0,
                138.0,
                140.0, // 10 Particles
                145.0,
                150.0,
                155.0,
                160.0,
                165.0, // 15 Particles
                172.0,
                180.0,
                186.0,
                193.0,
                198.0, // 20 Particles
        };
        */

        printf("# npart   0.0    0.5    1.0    1.5     2.0\n");
        for (uint npart = 2; npart < 21; ++npart) {
                FockStateSet set = fock_stateset_alloc_same_angular_momentum(0, cutoff_energies[npart], 500, npart);
                printf("%2u\t%u\n", npart, set.size);
                fock_stateset_free(&set);
                continue;
                double gs[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
                double *T = rgc_malloc(sizeof(*T)*set.size*set.size);
                double *V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);
                printf("%2u", npart);
                for (uint i = 0; i < sizeof(gs)/sizeof(gs[0]); ++i) {
                        double g = gs[i] * 5 / (npart-1);
                        EigenSolutionSet sol = fock_eigen_vt_solve(g, set.size, T, V);
                        double ground_energy_pparticle = sol.values[0] / npart;
                        printf("   %lf", ground_energy_pparticle);
                        fock_eigen_free(&sol);
                }
                putchar('\n');

                rgc_free(T);
                rgc_free(V);
                fock_stateset_free(&set);
        }

	return 0;
}

