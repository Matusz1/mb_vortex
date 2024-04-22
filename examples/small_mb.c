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
        uint tot_momentum = 0;
        //FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 120.0, 400, n_part);
        FockStateSet set = { 0 };
        FockState focks[5];

        focks[0] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
        });

        focks[1] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 2)),
        });

        focks[2] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct(-1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
        });

        focks[3] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 3)),
        });

        focks[4] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct(-1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 2)),
        });
        /*
        focks[0] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
        });

        focks[1] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 3, 1)),
        });

        focks[2] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
        });

        focks[3] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 0, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct(-1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 3, 1)),
        });

        focks[4] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct(-1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 2, 1)),
        });

        focks[5] = fock_init_direct(6, STATE_BESSEL, (uint[]){
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
                fock_bessel_find_index(fock_bessel_state_direct( 1, 1)),
        });
        */


        for (uint i = 0; i < sizeof(focks)/sizeof(focks[0]); ++i) {
                fock_stateset_append(&set, focks + i);
        }

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

        for (double gg = 0.0; gg < 3.0; gg += 0.1) {
                EigenSolutionSet sol = fock_eigen_solve(gg, set);
                printf("%lf\t%lf\n", gg, sol.values[0]);
                fock_eigen_free(&sol);
        }


        fock_stateset_free(&set);
	return 0;
}

