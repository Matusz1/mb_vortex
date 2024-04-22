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
        uint tot_momentum = 6;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 130.0, 400, n_part);

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

        Lattice2D lattice = fock_lattice_alloc((Vec2d){-1.0, -1.0}, (Vec2d){1.0, 1.0}, 50, 50);
        char fname[256];
        sprintf(fname, "./results/data/lambda_n%um%u.txt", n_part, tot_momentum);
        FILE* f = rgc_fopen(fname, "a");
        //FILE* f = stdout;
        for (double gg = 0.0; gg <= 0.201; gg += 0.01) {
                printf("g = %.1lf\n", gg);
                EigenSolutionSet sol = fock_eigen_vt_solve(gg, set.size, T, V);
                G1State gfunc = fock_g1_alloc(set, sol.vectors);

                complex double* matrix = fock_g1_matrix_alloc(lattice, &gfunc);
                double lambda = 0.0;
                double* diag = fock_g1_matrix_calc_top_eigenvector(matrix, lattice.size, 1000, &lambda);
                fprintf(f, "%.2lf %lf\n", gg, lambda);
                fflush(f);

                rgc_free(matrix);
                rgc_free(diag);
                fock_eigen_free(&sol);
                fock_g1_free(&gfunc);
        }
        rgc_fclose(f);

        fock_stateset_free(&set);
	return 0;
}

