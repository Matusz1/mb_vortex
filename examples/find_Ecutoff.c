#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>

#include <rgc_memory.h>
#include <rgc_io.h>

int main(int argc, char **argv) {
	double g = 0.5;
        uint idx = 0;
        double omega = 0.0;
        uint seed = 0;
        for (uint i = 1; i < argc; ++i) {
                if (strcmp(argv[i], "-g") == 0)
                        sscanf(argv[++i], "%lf", &g);
                else if (strcmp(argv[i], "-idx") == 0)
                        sscanf(argv[++i], "%u", &idx);
                else if (strcmp(argv[i], "-w") == 0)
                        sscanf(argv[++i], "%lf", &omega);
                else if (strcmp(argv[i], "-seed") == 0)
                        sscanf(argv[++i], "%u", &seed);
                else {
                        fprintf(stderr, "# ERROR: unknown argument %s.\n", argv[i]);
                        exit(EXIT_FAILURE);
                }
        }

        uint n_part = 6;
        uint tot_momentum = 6;
        //FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 80.0, 400, n_part);

        // 6:0 - 60:135
        // 6:6 - 80:145
        char fname[256];
        sprintf(fname, "./results/data/cutoffs/cutoff_n%u_m%u.txt", n_part, tot_momentum);
        FILE* file = rgc_fopen(fname, "w");
        fprintf(file, "Ecutoff size g0.5      g1.0      g1.5      g2.0      g3.0\n");
        for (double Ecut = 80.0; Ecut <= 145.5; Ecut += 5.0) {
                fprintf(stderr, "# Ecut = %.2lf\n", Ecut);
                FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, Ecut, 400, n_part);
                double* T = rgc_malloc(sizeof(*T)*set.size*set.size);
                double* V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);

                fprintf(file, "%7.1lf %4u", Ecut, set.size);
                double gs[] = { 0.5, 1.0, 1.5, 2.0, 3.0 };
                for (uint i = 0; i < sizeof(gs)/sizeof(*gs); ++i) {
                        EigenSolutionSet sol = fock_eigen_vt_solve(gs[i], set.size, T, V);
                        fprintf(file, " %lf", sol.values[0]);
                        fock_eigen_free(&sol);
                }
                fputc('\n', file);
                fflush(file);

                /* Free memory */
                rgc_free(T);
                rgc_free(V);
                fock_stateset_free(&set);
        }
        rgc_fclose(file);

	return 0;
}

