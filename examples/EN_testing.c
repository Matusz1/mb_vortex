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
        const double cutoff_energies[16] = {
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
        };

        printf("# npart   0.0    0.5    1.0    1.5     2.0\n");
        for (uint npart = 13; npart < 16; ++npart) {
                FockStateSet set = fock_stateset_alloc_same_angular_momentum(0, cutoff_energies[npart], 500, npart);
                //printf("%2u\t%u\n", npart, set.size);
                //fock_stateset_free(&set);
                //continue;
                double gs[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
                double *T = rgc_malloc(sizeof(*T)*set.size*set.size);
                double *V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);
                printf("%2u", npart);
                for (uint i = 0; i < sizeof(gs)/sizeof(gs[0]); ++i) {
                        double g = gs[i] * 6 / npart;
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

