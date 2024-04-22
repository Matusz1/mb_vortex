#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>

#include <rgc_memory.h>
#include <rgc_io.h>

int main(int argc, char **argv) {
	double g = 0.5;
        uint idx = 0;
        bool save_to_file = false;
        double omega = 0.0;
        uint seed = 0;
        for (uint i = 1; i < argc; ++i) {
                if (strcmp(argv[i], "-g") == 0)
                        sscanf(argv[++i], "%lf", &g);
                else if (strcmp(argv[i], "-idx") == 0)
                        sscanf(argv[++i], "%u", &idx);
                else if (strcmp(argv[i], "-w") == 0)
                        sscanf(argv[++i], "%lf", &omega);
                else if (strcmp(argv[i], "-s") == 0)
                        save_to_file = true;
                else if (strcmp(argv[i], "-seed") == 0)
                        sscanf(argv[++i], "%u", &seed);
                else {
                        fprintf(stderr, "# ERROR: unknown argument %s.\n", argv[i]);
                        exit(EXIT_FAILURE);
                }
        }

        uint n_part = 6;
        FockStateSet set = fock_stateset_alloc(STATE_BESSEL, 200, n_part, 80.0);

        double *V, *T;
        if (save_to_file) {
                T = rgc_malloc(sizeof(*T)*set.size*set.size);
                V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);
                util_bin_tofile(T, sizeof(*T)*set.size*set.size, "T.dat");
                util_bin_tofile(V, sizeof(*V)*set.size*set.size, "V.dat");
                return 0;
        } else {
                V = util_bin_alloc_fromfile(sizeof(*V)*set.size*set.size, "V.dat");
                T = util_bin_alloc_fromfile(sizeof(*T)*set.size*set.size, "T.dat");
        }

        for (double w = 5.0; w <= 7.0; w += 0.1) {
                fprintf(stderr, "%lf\n", w);
                EigenSolutionSet sol = fock_eigen_omega_vt_solve(g, w, set, T, V);
                uint i = 0;
                for (; i < sol.size; ++i)
                        if (sol.vectors[i] != 0.0)
                                break;
                int m = fock_bessel_angular_momentum(&set.states[i]);
                printf("%.3lf   %d\n", w, m);
                fflush(stdout);
                fock_eigen_free(&sol);
        }

        fock_stateset_free(&set);
	return 0;
}

