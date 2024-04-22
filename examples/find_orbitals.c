#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_heapsort.h>

#include <rgc_memory.h>
#include <rgc_io.h>

int compare_doubles(const double * a, const double * b) {
        if (*a > *b)
                return 1;
        else if (*a < *b)
                return -1;
        else
                return 0;
}

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
        uint tot_momentum = 3;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 130.0, 400, n_part);
        //FockStateSet set = fock_stateset_alloc(STATE_BESSEL, 200, n_part, 75.0);

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

        EigenSolutionSet sol = fock_eigen_vt_solve(g, set.size, T, V);

        double orbital_weights[400] = { 0.0 };
        for (uint i = 0; i < sol.size; ++i) {
                const double w = sol.vectors[i] * sol.vectors[i];
                for (uint j = 0; j < n_part; ++j) {
                        orbital_weights[set.states[i].states[j]] += w / 6.0;
                }
        }
        size_t s_index[400];
        gsl_heapsort_index(s_index, orbital_weights, 400, sizeof(*orbital_weights), (gsl_comparison_fn_t) compare_doubles);
        //fock_stateset_print_weighted(set, sol.vectors, 20);
        printf("%lf\n", g);
        for (uint i = 0; i < 10; ++i) {
                printf("%lf ", orbital_weights[s_index[399-i]]);
                fock_single_particle_print(s_index[399-i], STATE_BESSEL);
                putchar('\n');
        }
        putchar('\n');

        fock_stateset_free(&set);
	return 0;
}

