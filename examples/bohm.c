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
        uint tot_momentum = 3;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 120.0, 500, n_part);

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

        FockCoeffStateSet cset = fock_coeffstateset_alloc_dcoeff(set, sol.vectors, 0.00001);
        seed = 8748;
        srand(seed);
        double* pos = fock_coeffstateset_alloc_random(cset, 1);
        double dt = 0.01;
        uint niter = 1000;
        double path[2*n_part*niter];
        for (uint i = 0; i < niter; ++i) {
                printf("%5u / %5u\n", i+1, niter);
                memcpy(path+(2*n_part)*i, pos, 2*n_part*sizeof(*pos));
                fock_stateset_bohm_next(set, pos, pos+6, sol.vectors, dt);
        }
        char output_fname[256];
        sprintf(output_fname, "./final/data/paths/path_m%u.dat", tot_momentum);
        util_bin_tofile(path, sizeof(*path)*2*n_part*niter, output_fname);

        fock_stateset_free(&set);
	return 0;
}

